#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

#include "miniaudio.h"

#include <SDL.h>

#ifdef BUILD_LINUX
#define GL_GLEXT_PROTOTYPES
#include <GL/gl.h>
#elif BUILD_MACOS
#define GL_SILENCE_DEPRECATION
#include <OpenGL/gl3.h>
#else
#error "missing define for gl.h"
#endif

#ifdef BUILD_LINUX
#define NANOVG_GLES3_IMPLEMENTATION
#elif BUILD_MACOS
#define NANOVG_GL3_IMPLEMENTATION
#else
#error "missing BUILD_* define"
#endif

#include <nanovg.h>
#include <nanovg_gl.h>
#include <nanovg_gl_utils.h>

#include <FLAC/metadata.h>
#include <FLAC/stream_decoder.h>
#include <FLAC/stream_encoder.h>

#define STB_SPRINTF_IMPLEMENTATION
#include "stb_sprintf.h"

#define PP_IMPLEMENTATION
#include "pitch_prefect.h"

#define OPT_IMPLEMENTATION
#include "opt.h"


//#define DEFAULT_NOISE_GATE_RMS_THRESHOLD (1e-3)
#define DEFAULT_NOISE_GATE_RMS_THRESHOLD (0)

static inline int clampi(int v, int min, int max)
{
	return v > max ? max : (v < min ? min : v);
}

struct PP pp;

struct RecSimple {
	float f0;
	float expected_f0;
	float rms;
	float nsdf_best_frequency;
	float spectrogram_corrected_best_frequency;
	// PP_state if I ever get one?
};

struct Rec {
	int cap;
	int n;
	SDL_mutex* n_lock;
	struct RecSimple* simple;
	float* float_data; // window, spectrogram, nsdf?
};

#define REC_N_SECTIONS (3)
#define REC_STRIDE (REC_N_SECTIONS * PP_WINDOW_SIZE)

static void Rec_push_pp(struct Rec* rec, struct PP* pp)
{
	int n = rec->n;
	assert(n < rec->cap);

	struct RecSimple* simple = &rec->simple[n] ;
	simple->f0 = pp->f0;
	simple->rms = pp->rms;
	simple->nsdf_best_frequency = pp->nsdf_best_frequency;
	simple->spectrogram_corrected_best_frequency = pp->spectrogram_corrected_best_frequency;

	float* float_data = &rec->float_data[n * REC_STRIDE];
	size_t sz = PP_WINDOW_SIZE * sizeof(*float_data);
	memcpy(&float_data[0*PP_WINDOW_SIZE], pp->window, sz);
	memcpy(&float_data[1*PP_WINDOW_SIZE], pp->spectrogram, sz); // spectrogram is only half-size, but whatever...
	memcpy(&float_data[2*PP_WINDOW_SIZE], pp->nsdf, sz);

	SDL_LockMutex(rec->n_lock);
	rec->n++;
	SDL_UnlockMutex(rec->n_lock);
}

static int Rec_n(struct Rec* rec)
{
	SDL_LockMutex(rec->n_lock);
	int n = rec->n;
	SDL_UnlockMutex(rec->n_lock);
	return n;
}

static void Rec_crop(struct Rec* rec, int i0, int n)
{
	assert(i0 >= 0);
	assert(n > 0);
	assert(i0+n <= rec->n);
	rec->n = n;
	rec->cap -= i0;
	rec->simple = &rec->simple[i0];
	rec->float_data = &rec->float_data[i0 * REC_STRIDE];
}

static void Rec_fill_expected_f0(struct Rec* rec)
{
	for (int i = 0; i < rec->n; i++) {
		struct RecSimple* s = &rec->simple[i];
		s->expected_f0 = s->f0;
	}
}

static float Rec_usage(struct Rec* rec)
{
	return (float)Rec_n(rec) / (float)rec->cap;
}

#define REC_WINDOW 0
#define REC_SPECTROGRAM 1
#define REC_NSDF 2
static float* Rec_get_float_data(struct Rec* rec, int type, int index)
{
	return &rec->float_data[REC_STRIDE * index + PP_WINDOW_SIZE * type];
}

struct Rec rec;

static void Rec_init_with_cap(struct Rec* rec, int cap)
{
	memset(rec, 0, sizeof *rec);

	rec->cap = cap;

	assert((rec->n_lock = SDL_CreateMutex()) != NULL);

	size_t sz = 0;

	size_t simple_sz = cap * sizeof(*rec->simple);
	assert((rec->simple = malloc(simple_sz)) != NULL);
	sz += simple_sz;

	size_t float_data_sz = cap * REC_STRIDE * sizeof(*rec->float_data);
	assert((rec->float_data = malloc(float_data_sz)) != NULL);
	sz += float_data_sz;

	printf("Rec_init_with_cap: cap=%d alloc'd=%zd\n", cap, sz);
}

// XXX can't be bothered to use frexpf/ldexpf to make a portable
// encoder/decoder :D greetings if you're reading this comment on a big-endian
// machine!

static inline uint32_t f32_encode(float v)
{
	union { uint32_t u; float f; } x;
	x.f = v;
	return x.u;
}

static inline float f32_decode(uint32_t v)
{
	union { uint32_t u; float f; } x;
	x.u = v;
	return x.f;
}

static inline void readn(FILE* f, void* buf, size_t n_bytes)
{
	if (fread(buf, 1, n_bytes, f) != n_bytes) {
		fprintf(stderr, "failed to read %zd bytes\n", n_bytes);
		exit(EXIT_FAILURE);
	}
}

char* ptc_fourcc = "PTC0";

char buf_fourcc[5];
static inline char* read_fourcc(FILE* f)
{
	readn(f, buf_fourcc, 4);
	buf_fourcc[4] = 0;
	return buf_fourcc;
}

static inline uint32_t read_u32(FILE* f)
{
	uint32_t v;
	readn(f, &v, sizeof v);
	// XXX fix endianness (lol)
	return v;
}

static inline float read_f32(FILE *f)
{
	return f32_decode(read_u32(f));
}

struct FlacReadData {
	int expected_sample_rate;
};

static FLAC__StreamDecoderWriteStatus flac_read_callback(const FLAC__StreamDecoder *decoder, const FLAC__Frame *frame, const FLAC__int32 * const buffer[], void *client_data)
{
	struct FlacReadData* data = client_data;
	assert(frame->header.sample_rate == data->expected_sample_rate);
	assert(frame->header.channels == 1);
	assert(frame->header.bits_per_sample == 16);

	int sz = frame->header.blocksize;
	const FLAC__int32* ch0 = buffer[0];
	for (int i = 0; i < sz; i++) {
		FLAC__int32 sample = ch0[i];
		float x = (float)sample / 32768.0f;
		if (pp_write_one_decimated(&pp, x)) {
			Rec_push_pp(&rec, &pp);
		}
	}

	//printf("FRAME sz=%d sample=%ld\n", frame->header.blocksize, frame->header.number.sample_number);
	return FLAC__STREAM_DECODER_WRITE_STATUS_CONTINUE;

}

void flac_read_error_callback(const FLAC__StreamDecoder *decoder, FLAC__StreamDecoderErrorStatus status, void *client_data)
{
	fprintf(stderr, "Got error callback: %s\n", FLAC__StreamDecoderErrorStatusString[status]);
}


#define PTC_FLOAT_SANITY_CHECK (420.024f)

static void Rec_load(struct Rec* rec, const char* filename)
{
	FILE* f = fopen(filename, "rb");
	if (f == NULL) {
		fprintf(stderr, "%s: no such file\n", filename);
		exit(EXIT_FAILURE);
	}

	char* fourcc = read_fourcc(f);
	if (strcmp(fourcc, ptc_fourcc) != 0) {
		fprintf(stderr, "%s: not a .ptc file (fourcc mismatch)\n", filename);
		exit(EXIT_FAILURE);
	}

	{
		float ptc_float_sanity_check = read_f32(f);
		if (ptc_float_sanity_check != PTC_FLOAT_SANITY_CHECK) {
			fprintf(stderr, "ERROR: float sanity check failed; expected %f; got %f\n", PTC_FLOAT_SANITY_CHECK, ptc_float_sanity_check);
			abort();
		}
	}

	// read header
	int n = read_u32(f);
	uint32_t expected_pp_bucket_size = read_u32(f);
	/* reading, but ignoring -> */ read_u32(f); /* ... used to be expected
	PP "n buckets" value, but we only care about values that change
	appearance outwards; changing PP_BUCKET_SIZE changes the PP output
	frequency, but changing "n buckets" didn't */
	uint32_t pp_decimation_factor = read_u32(f);
	float pp_decimated_sample_rate_hz = read_f32(f);
	float pp_noise_gate_rms_threshold = read_f32(f);

	// init rec
	Rec_init_with_cap(rec, n);

	// read expected f0 array
	for (int i = 0; i < n; i++) {
		rec->simple[i].expected_f0 = read_f32(f);
	}

	// check that file was written by a build having the same "master configuration"
	assert(expected_pp_bucket_size == PP_BUCKET_SIZE);

	pp_init(&pp, pp_decimated_sample_rate_hz * pp_decimation_factor);
	pp.noise_gate_rms_threshold = pp_noise_gate_rms_threshold;
	assert(pp.decimated_sample_rate_hz == pp_decimated_sample_rate_hz);
	assert(pp.decimation_factor == pp_decimation_factor);

	// flac read
	{
		FLAC__bool ok = true;
		FLAC__StreamDecoder *decoder = 0;
		FLAC__StreamDecoderInitStatus init_status;

		if ((decoder = FLAC__stream_decoder_new()) == NULL) {
			fprintf(stderr, "ERROR: flac decoder allocation failed\n");
			abort();
		}

		FLAC__stream_decoder_set_md5_checking(decoder, true);

		struct FlacReadData data;
		data.expected_sample_rate = pp_decimated_sample_rate_hz;

		init_status = FLAC__stream_decoder_init_FILE(decoder, f, flac_read_callback, NULL, flac_read_error_callback, &data);
		if (init_status != FLAC__STREAM_DECODER_INIT_STATUS_OK) {
			fprintf(stderr, "ERROR: flac decoder initialization failed: %s\n", FLAC__StreamDecoderInitStatusString[init_status]);
			abort();
		}

		ok = FLAC__stream_decoder_process_until_end_of_stream(decoder);
		if (!ok) {
			fprintf(stderr, "ERROR: flac decode error\n");
			abort();
		}

		FLAC__stream_decoder_delete(decoder);
	}

	assert(rec->n == n);
}

static void Rec_reload(struct Rec* rec)
{
	int n = rec->n;
	pp_init(&pp, pp.decimated_sample_rate_hz * pp.decimation_factor);
	rec->n = 0;
	for (int i = 0; i < n; i++) {
		float* window = Rec_get_float_data(rec, REC_WINDOW, i);
		for (int j = 0; j < PP_BUCKET_SIZE; j++) {
			float x = window[PP_WINDOW_SIZE - PP_BUCKET_SIZE + j];
			if (pp_write_one_decimated(&pp, x)) {
				Rec_push_pp(rec, &pp);
			}
		}
	}
	assert(rec->n == n);
}

static inline void writen(FILE* f, void* bytes, size_t sz)
{
	if (fwrite(bytes, 1, sz, f) != sz) {
		fprintf(stderr, "failed to write %zd bytes\n", sz);
		exit(EXIT_FAILURE);
	}
}

static inline void write_u32(FILE* f, uint32_t v)
{
	writen(f, &v, sizeof v);
}

static inline void write_f32(FILE* f, float v)
{
	write_u32(f, f32_encode(v));
}

static void Rec_save(struct Rec* rec, const char* filename)
{
	FILE* f = fopen(filename, "wb");
	if (f == NULL) {
		fprintf(stderr, "%s: cannot open for writing\n", filename);
		exit(EXIT_FAILURE);
	}

	writen(f, ptc_fourcc, 4);

	write_f32(f, PTC_FLOAT_SANITY_CHECK);

	// write header
	write_u32(f, rec->n);
	write_u32(f, PP_BUCKET_SIZE);
	write_u32(f, 0); // deprecated; see Rec_load()
	write_u32(f, pp.decimation_factor);
	write_f32(f, pp.decimated_sample_rate_hz);
	write_f32(f, pp.noise_gate_rms_threshold);

	for (int i = 0; i < rec->n; i++) write_f32(f, rec->simple[i].expected_f0);

	{
		FLAC__bool ok = true;
		FLAC__StreamEncoder *encoder = 0;
		FLAC__StreamEncoderInitStatus init_status;

		if((encoder = FLAC__stream_encoder_new()) == NULL) {
			fprintf(stderr, "ERROR: flac encoder allocation failed\n");
			abort();
		}

		const int channels = 1;
		const int bits_per_sample = 16;
		const int sample_rate = pp.decimated_sample_rate_hz;

		ok &= FLAC__stream_encoder_set_verify(encoder, true);
		ok &= FLAC__stream_encoder_set_compression_level(encoder, 5);
		ok &= FLAC__stream_encoder_set_channels(encoder, channels);
		ok &= FLAC__stream_encoder_set_bits_per_sample(encoder, bits_per_sample);
		ok &= FLAC__stream_encoder_set_sample_rate(encoder, sample_rate);

		if (!ok) {
			fprintf(stderr, "ERROR: flac write setup failed\n");
			abort();
		}

		init_status = FLAC__stream_encoder_init_FILE(encoder, f, NULL, NULL);
		if (init_status != FLAC__STREAM_ENCODER_INIT_STATUS_OK) {
			fprintf(stderr, "ERROR: flac encoder initialization failed: %s\n", FLAC__StreamEncoderInitStatusString[init_status]);
			abort();
		}

		for (int i = 0; i < rec->n; i++) {
			float* window = Rec_get_float_data(rec, REC_WINDOW, i);
			FLAC__int32 bucket[PP_BUCKET_SIZE];
			int bucket_index = 0;
			for (int i = PP_WINDOW_SIZE-PP_BUCKET_SIZE; i < PP_WINDOW_SIZE; i++, bucket_index++) {
				float x = window[i];
				bucket[bucket_index] = clampi(x * 32768, -32768, 32767);
			}
			ok = FLAC__stream_encoder_process_interleaved(encoder, bucket, PP_BUCKET_SIZE);
			if (!ok) {
				fprintf(stderr, "ERROR: flac encoding error\n");
				abort();
			}
		}

		ok &= FLAC__stream_encoder_finish(encoder);
		if (!ok) {
			fprintf(stderr, "ERROR: flac encoder finish\n");
			abort();
		}

		FLAC__stream_encoder_delete(encoder);
	}

	printf("wrote %s\n", filename);
}

SDL_mutex* is_recording_lock;
static int is_recording;

SDL_mutex* playback_lock;
int playback_start = 0;
int playback_position = 0;
int playback_i0 = 0;
int playback_i1 = 0;
int playback_position_frame = 0;

static void audio_callback(ma_device* device, void* v_output, const void* v_input, ma_uint32 n_frames)
{
	// input
	{
		const float* input = v_input;

		SDL_LockMutex(is_recording_lock);
		int copy_is_recording = is_recording;
		SDL_UnlockMutex(is_recording_lock);
		if (copy_is_recording) {
			const float* ip = input;
			int n_channels = device->capture.channels;
			for (int frame = 0; frame < (int)n_frames; frame++) {
				float x = 0;
				for (int i = 0; i < n_channels; i++) {
					x += *(ip++);
				}
				if (pp_write_one(&pp, x)) {
					Rec_push_pp(&rec, &pp);
				}
			}
		}
	}

	// output
	{
		float* output = v_output;

		SDL_LockMutex(playback_lock);
		int copy_playback_start = playback_start;
		int copy_playback_position = playback_position;
		int copy_playback_i0 = playback_i0;
		int copy_playback_i1 = playback_i1;
		SDL_UnlockMutex(playback_lock);

		int stop = 0;
		int dpp = 0;

		float* ip = output;
		int n_channels = device->playback.channels;

		for (int frame = 0; frame < (int)n_frames; frame++) {
			float sample = 0.0;
			if (copy_playback_start) {
				int ppos = (copy_playback_position + dpp) / pp.decimation_factor;
				int bucket_index = copy_playback_i0 + (ppos / PP_BUCKET_SIZE);
				int bucket_offset = ppos % PP_BUCKET_SIZE;
				if (bucket_index < copy_playback_i1) {
					sample = Rec_get_float_data(&rec, REC_WINDOW, bucket_index)[PP_WINDOW_SIZE - PP_BUCKET_SIZE + bucket_offset];
				} else {
					stop = 1;
				}
				dpp++;
			}

			for (int i = 0; i < n_channels; i++) {
				*(ip++) = sample;
			}
		}

		SDL_LockMutex(playback_lock);
		if (stop) playback_start = 0;
		playback_position += dpp;
		playback_position_frame = copy_playback_i0 + playback_position / (pp.decimation_factor * PP_BUCKET_SIZE);
		SDL_UnlockMutex(playback_lock);
	}
}

ma_device audio_device;

static void audio_init()
{
	assert((is_recording_lock = SDL_CreateMutex()) != NULL);
	assert((playback_lock = SDL_CreateMutex()) != NULL);

	const int sample_rate = 48000;

	ma_device_config config = ma_device_config_init(ma_device_type_duplex);
	config.capture.format = ma_format_f32;
	config.capture.channels = 1;
	config.playback.format = ma_format_f32;
	config.playback.channels = 1;
	config.sampleRate = sample_rate;
	config.dataCallback = audio_callback;

	if (ma_device_init(NULL, &config, &audio_device) != MA_SUCCESS) {
		fprintf(stderr, "ma_device_init() failed\n");
		exit(EXIT_FAILURE);
	}

	pp_init(&pp, sample_rate);
	pp.noise_gate_rms_threshold = DEFAULT_NOISE_GATE_RMS_THRESHOLD ;

	ma_device_start(&audio_device);
}

SDL_Window* window;
SDL_GLContext glctx;
int screen_width;
int screen_height;
float pixel_ratio;
NVGcontext* vg;
int monofont;

static void window_size(int* width, int* height, float* pixel_ratio)
{
	int prev_width = *width;
	int prev_height = *height;
	SDL_GL_GetDrawableSize(window, width, height);
	if ((*width != prev_width || *height != prev_height)) {
		printf("%d×%d -> %d×%d\n", prev_width, prev_height, *width, *height);
	}

	int w, h;
	SDL_GetWindowSize(window, &w, &h);
	*pixel_ratio = *width / w;
}

static void video_init()
{
	#ifdef BUILD_LINUX
	SDL_GL_SetAttribute(SDL_GL_CONTEXT_MAJOR_VERSION, 3);
	SDL_GL_SetAttribute(SDL_GL_CONTEXT_MINOR_VERSION, 0);
	SDL_GL_SetAttribute(SDL_GL_CONTEXT_PROFILE_MASK, SDL_GL_CONTEXT_PROFILE_ES);
	#elif BUILD_MACOS
	SDL_GL_SetAttribute(SDL_GL_CONTEXT_MAJOR_VERSION, 3);
	SDL_GL_SetAttribute(SDL_GL_CONTEXT_MINOR_VERSION, 0);
	SDL_GL_SetAttribute(SDL_GL_CONTEXT_PROFILE_MASK, SDL_GL_CONTEXT_PROFILE_CORE);
	#else
	#error "missing BUILD_* define"
	#endif

	SDL_GL_SetAttribute(SDL_GL_STENCIL_SIZE, 1);

	window = SDL_CreateWindow(
			"pptool",
			SDL_WINDOWPOS_UNDEFINED, SDL_WINDOWPOS_UNDEFINED,
			1920, 1080,
			SDL_WINDOW_RESIZABLE | SDL_WINDOW_OPENGL | SDL_WINDOW_ALLOW_HIGHDPI);
	if (window == NULL) {
		fprintf(stderr, "SDL_CreateWindow failed: %s\n", SDL_GetError());
		abort();
	}
	glctx = SDL_GL_CreateContext(window);
	if (!glctx) {
		fprintf(stderr, "SDL_GL_CreateContextfailed: %s\n", SDL_GetError());
		abort();
	}

	//const int flags = NVG_ANTIALIAS | NVG_STENCIL_STROKES | NVG_DEBUG;
	const int flags = NVG_ANTIALIAS | NVG_STENCIL_STROKES;
	#ifdef BUILD_LINUX
		vg = nvgCreateGLES3(flags);
	#elif BUILD_MACOS
		vg = nvgCreateGL3(flags);
	#else
	#error "missing BUILD_* define"
	#endif

	assert(vg != NULL);

	monofont = nvgCreateFont(vg, "mono", "../3rd/VeraMono.ttf");
	assert(monofont != -1);

	window_size(&screen_width, &screen_height, &pixel_ratio);

	SDL_GL_SetSwapInterval(1);
}

static void video_exit()
{
	SDL_GL_DeleteContext(glctx);
	SDL_DestroyWindow(window);
}

static inline float f2p(float f)
{
	return log2f(f / 440.0) * 12.0;
}

char* prg;

static void usage()
{
	fprintf(stderr, "usage: %s [-h] [infile.ptc]\n", prg);
	fprintf(stderr, " -h / --help           help\n");
	exit(EXIT_FAILURE);
}

int main(int argc, char** argv)
{
	prg = argv[0];

	struct Opt opt;
	struct OptDef defs[] = {
		{OPT_FLAG, 'h', "help"},
		{0},
	};
	opt_init(&opt, defs, argc-1, argv+1);

	char* load_filename = NULL;

	while (opt_next(&opt)) {
		if (opt.is_invalid) {
			fprintf(stderr, "%s: %s\n\n", opt.arg, opt.errmsg);
			usage();
		} else if (opt.is_switch) {
			switch (opt.short_opt) {
			case 'h':
				usage();
				break;
			default:
				OPT_UNREACHABLE;
			}
		} else if (opt.is_npos) {
			if (load_filename == NULL) {
				load_filename = opt.value;
			} else {
				fprintf(stderr, "%s: unexpected non-positional argument\n", opt.value);
				usage();
				break;
			}
		} else {
			OPT_UNREACHABLE;
		}
	}

	int mode = 0;

	char* save_filename;
	if (load_filename == NULL) {
		const int n_rec_minutes = 10; // TEN MINUTES OR BUST!
		const float approx_frames_per_second = (float)PP_IDEAL_SAMPLE_RATE_HZ / (float)PP_BUCKET_SIZE;
		const float approx_frames_per_minute = approx_frames_per_second * 60.0;
		Rec_init_with_cap(&rec, n_rec_minutes * approx_frames_per_minute);
		save_filename = "RENAME-ME-pptool-save-file.ptc";
	} else {
		Rec_load(&rec, load_filename);
		save_filename = load_filename;
		mode = 2;
	}

	is_recording = (mode == 0);

	audio_init();

	assert(SDL_Init(SDL_INIT_TIMER | SDL_INIT_VIDEO) == 0);
	atexit(SDL_Quit);

	video_init();

	// scrolling spectrogram: ss/SS (python joke)
	#define SS_WIDTH (1920)
	#define SS_HEIGHT PP_SPECTROGRAM_SZ
	//const int ss_image_flags = NVG_IMAGE_REPEATX | NVG_IMAGE_NEAREST;
	const int ss_image_flags = NVG_IMAGE_NEAREST;
	int ss_handle = nvgCreateImageRGBA(vg, SS_WIDTH, SS_HEIGHT, ss_image_flags, NULL);
	uint32_t* ss_data = calloc(SS_WIDTH * SS_HEIGHT, sizeof *ss_data);
	assert(ss_data != NULL);
	NVGpaint ss_paint = nvgImagePattern(vg, 0, 0, SS_WIDTH, SS_HEIGHT, 0, ss_handle, 1.0);

	const float preset_noise_gate_rms_threshold = pp.noise_gate_rms_threshold;

	int exiting = 0;
	int fullscreen = 0;
	int mx = 0;
	//int my = 0;
	int frame = 0;
	int crop_i0 = 0;
	int crop_i1 = 0;
	int crop_i0_norm = 0;
	int crop_i1_norm = 0;
	int drag = 0;
	int case_frame = 0;
	int show_f0 = 1;
	int show_expected_f0 = 1;
	while (!exiting) {
		SDL_Event e;
		int set_crop_i0 = 0;
		int set_crop_i1 = 0;
		int no_nudge_this_frame = 1;
		while (SDL_PollEvent(&e)) {
			if (e.type == SDL_QUIT) {
				exiting = 1;
			} else if (e.type == SDL_KEYDOWN) {
				int sym = e.key.keysym.sym;
				if (sym == SDLK_ESCAPE || sym == 'q') {
					exiting = 1;
				} else if (sym == 'f') {
					fullscreen = !fullscreen;
					//SDL_SetWindowFullscreen(window, fullscreen ? SDL_WINDOW_FULLSCREEN : 0);
					SDL_SetWindowFullscreen(window, fullscreen ? SDL_WINDOW_FULLSCREEN_DESKTOP : 0);
				}

				if (mode == 0) {
					if (sym == SDLK_RETURN) {
						SDL_LockMutex(is_recording_lock);
						is_recording = 0;
						SDL_UnlockMutex(is_recording_lock);
						mode = 1;
					}
				} else if (mode == 1) {
					if (sym == SDLK_RETURN) {
						Rec_crop(&rec, crop_i0_norm, crop_i1_norm-crop_i0_norm);
						Rec_fill_expected_f0(&rec);
						Rec_reload(&rec);
						// TODO sæt selection til hele Rec?
						// XXX ehh... ved hop til mode=2 ville det nok være bedst hvis vi reloadede alt data...
						mode = 2;
					} else if (sym == '[') {
						set_crop_i0 = 1;
					} else if (sym == ']') {
						set_crop_i1 = 1;
					}
				} else if (mode == 2) {
					if (sym == 's') {
						Rec_save(&rec, save_filename);
					} else if (sym == 'r') {
						Rec_reload(&rec);
						printf("reloaded...\n");
					} else if (sym == 'z') {
						pp.noise_gate_rms_threshold = preset_noise_gate_rms_threshold;
						Rec_reload(&rec);
					} else if (sym == 'n') {
						Rec_fill_expected_f0(&rec);
					} else if (sym == SDLK_LEFT) {
						case_frame--;
					} else if (sym == SDLK_RIGHT) {
						case_frame++;
					} else if (sym == '1') {
						show_f0 = !show_f0;
					} else if (sym == '2') {
						show_expected_f0 = !show_expected_f0;
					}
				}

				{
					const float nudge = 2e-4;
					int nudged = 0;
					if (no_nudge_this_frame) {
						if (sym == SDLK_UP) {
							pp.noise_gate_rms_threshold += nudge;
							nudged = 1;
						} else if (sym == SDLK_DOWN) {
							pp.noise_gate_rms_threshold -= nudge;
							if (pp.noise_gate_rms_threshold < 0) pp.noise_gate_rms_threshold = 0;
							nudged = 1;
						}
					}
					if (nudged && mode > 0) {
						Rec_reload(&rec);
						no_nudge_this_frame = 0;
					}
				}

				if ((mode == 1 || mode == 2) && sym == SDLK_SPACE) {
					SDL_LockMutex(playback_lock);
					playback_start = !playback_start;
					playback_position = 0;
					playback_i0 = playback_i1 = 0;
					if (mode == 1) {
						playback_i0 = crop_i0_norm;
						playback_i1 = crop_i1_norm;
					} else if (mode == 2) {
						playback_i0 = case_frame;
						playback_i1 = rec.n;
					}
					//printf("%s [%d; %d]\n", playback_start ? "PLAY" : "STOP", playback_i0, playback_i1);
					SDL_UnlockMutex(playback_lock);
				}
			} else if (e.type == SDL_MOUSEMOTION) {
				mx = e.motion.x;
				//my = e.motion.y;
			} else if (e.type == SDL_MOUSEBUTTONDOWN || e.type == SDL_MOUSEBUTTONUP) {
				if (mode == 1 || mode == 2) drag++;
				//e.button;
			} else if (e.type == SDL_WINDOWEVENT) {
				if (e.window.event == SDL_WINDOWEVENT_RESIZED) {
					window_size(&screen_width, &screen_height, &pixel_ratio);
				}
			}
		}

		glViewport(0, 0, screen_width, screen_height);
		glClearColor(0, 0.05, 0.15, 0);
		glClear(GL_COLOR_BUFFER_BIT | GL_STENCIL_BUFFER_BIT);

		glEnable(GL_BLEND);
		glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
		glEnable(GL_CULL_FACE);
		glDisable(GL_DEPTH_TEST);

		nvgBeginFrame(vg, screen_width / pixel_ratio, screen_height / pixel_ratio, pixel_ratio);
		nvgGlobalCompositeOperation(vg, NVG_LIGHTER); // additive color blending

		// plot usage: we crash when the Rec buffer is filled up :D
		if (mode == 0) {
			nvgSave(vg);
			float usage = Rec_usage(&rec);

			const float h = 4;
			nvgBeginPath(vg);
			nvgRect(vg, 0, 0, screen_width, h);
			nvgFillColor(vg, nvgRGB(0, 100, 200));
			nvgFill(vg);

			nvgBeginPath(vg);
			nvgRect(vg, 0, 0, screen_width * usage, h);
			nvgFillColor(vg, nvgRGB(100, 200, 255));
			nvgFill(vg);

			nvgRestore(vg);
		}

		const int window_width = screen_width;
		float xstep = (float)screen_width / (float)window_width;

		int rec_n = Rec_n(&rec);
		int i1 = rec_n - 1;
		int i0 = i1 - window_width;
		int mouse_index = i0 + ((float)mx * xstep);
		if (i0 < 0) i0 = 0;
		const float x0 = screen_width - (i1-i0)*xstep;

		if (drag == 1) {
			if (mode == 1) set_crop_i0 = 1;
			drag++;
		}
		if (drag == 2 && mode == 1) set_crop_i1 = 1;
		if (drag == 3) drag = 0;
		if (set_crop_i0) crop_i0 = mouse_index;
		if (set_crop_i1) crop_i1 = mouse_index;

		if (mode == 2 && drag > 0) {
			case_frame = mouse_index;
		}

		crop_i0_norm = clampi(crop_i0, 0, i1);
		crop_i1_norm = clampi(crop_i1, 0, i1) + 1;
		if (crop_i0_norm > crop_i1_norm) {
			int tmp = crop_i0_norm;
			crop_i0_norm = crop_i1_norm;
			crop_i1_norm = tmp;
		}

		if (case_frame < 0) case_frame = 0;
		if (case_frame > i1) case_frame = i1;

		float stack_y = 25;

		int windows_index = mode < 2 ? i1 : case_frame;
		SDL_LockMutex(playback_lock);
		if (playback_start) windows_index = playback_position_frame;
		SDL_UnlockMutex(playback_lock);

		const float scope_height = 200;
		const float scope_width = screen_width*0.5;
		// plot NSDF
		{
			const float nsdf_height = scope_height;
			float* nsdf = Rec_get_float_data(&rec, REC_NSDF, windows_index);
			nvgSave(vg);

			nvgTranslate(vg, 0, stack_y + nsdf_height/2);
			nvgBeginPath(vg);
			nvgMoveTo(vg, 0, 0);
			for (int i = 0; i < PP_WINDOW_SIZE; i++) {
				float x = ((float)i / (float)PP_WINDOW_SIZE) * scope_width;
				float y = nsdf[i] * -nsdf_height/2;
				nvgLineTo(vg, x, y);
			}
			nvgLineTo(vg, scope_width, 0);
			nvgClosePath(vg);

			nvgFillColor(vg, nvgRGBA(0, 255, 0, 50));
			nvgFill(vg);

			nvgStrokeColor(vg, nvgRGBA(0, 255, 0, 220));
			nvgStrokeWidth(vg, 2.0);
			nvgStroke(vg);

			float nsdf_best_frequency = rec.simple[windows_index].nsdf_best_frequency;
			if (nsdf_best_frequency > 0) {
				float x = ((((float)pp.decimated_sample_rate_hz / (float)nsdf_best_frequency)) / (float)PP_WINDOW_SIZE) * scope_width;
				nvgBeginPath(vg);
				const float w = 1;
				nvgRect(vg, x-w*0.5, -scope_height*0.5, w, scope_height);
				nvgFillColor(vg, nvgRGBA(100, 255, 100, 100));
				nvgFill(vg);
			}

			nvgRestore(vg);
		}

		// plot spectrogram
		{
			const float spectrogram_height = scope_height;
			float* spectrogram = Rec_get_float_data(&rec, REC_SPECTROGRAM, windows_index);
			nvgSave(vg);

			nvgTranslate(vg, scope_width, stack_y + spectrogram_height);
			nvgBeginPath(vg);
			nvgMoveTo(vg, 0, 0);
			for (int i = 0; i < PP_SPECTROGRAM_SZ; i++) {
				float x = ((float)i / (float)PP_SPECTROGRAM_SZ) * scope_width;
				float y = spectrogram[i] * 6e-2 * -spectrogram_height;
				nvgLineTo(vg, x, y);
			}
			nvgLineTo(vg, scope_width, 0);
			nvgClosePath(vg);

			nvgFillColor(vg, nvgRGBA(255, 50, 0, 50));
			nvgFill(vg);

			nvgStrokeColor(vg, nvgRGBA(255, 50, 0, 220));
			nvgStrokeWidth(vg, 2.0);
			nvgStroke(vg);

			float spectrogram_corrected_best_frequency = rec.simple[windows_index].spectrogram_corrected_best_frequency;
			if (spectrogram_corrected_best_frequency > 0) {
				nvgBeginPath(vg);
				for (int i = 0; i < PP_SPECTROGRAM_SZ; i++) {
					float x = ((float)(i) / (float)PP_SPECTROGRAM_SZ) * scope_width;
					const float fft_frequency_step = (((float)pp.decimated_sample_rate_hz * 0.5f) / (float)PP_SPECTROGRAM_SZ);
					float y = -pp__spectrogram_scoring_function((float)(i+1) * (fft_frequency_step / spectrogram_corrected_best_frequency)) * 20;
					if (i == 0) {
						nvgMoveTo(vg, x, y);
					} else {
						nvgLineTo(vg, x, y);
					}
				}
				nvgStrokeColor(vg, nvgRGBA(255, 255, 255, 180));
				nvgStrokeWidth(vg, 1.0);
				nvgStroke(vg);
			}

			nvgRestore(vg);
			//stack_y += spectrogram_height;
		}

		stack_y += scope_height;
		int scroller_begin_y = stack_y;

		// plot RMS
		{
			nvgSave(vg);
			const float rms_height = 50;
			nvgTranslate(vg, 0, stack_y);
			nvgBeginPath(vg);
			nvgMoveTo(vg, x0, rms_height);
			float x = x0;
			const float rms_scale = 1.0f;
			for (int i = i0; i < i1; i++) {
				struct RecSimple* s = &rec.simple[i];
				float rms = s->rms;
				float y = rms_height - rms * rms_height * rms_scale;
				nvgLineTo(vg, x, y);
				x += xstep;
			}
			nvgLineTo(vg, screen_width, rms_height);
			nvgClosePath(vg);

			nvgFillColor(vg, nvgRGB(255, 255, 255));
			nvgFill(vg);

			// show noise gate
			{
				// TODO show "above rms" in different color, hmm?
				float y = rms_height - pp.noise_gate_rms_threshold * rms_height * rms_scale;
				nvgBeginPath(vg);
				nvgMoveTo(vg, 0, y);
				nvgLineTo(vg, screen_width, y);
				nvgStrokeColor(vg, nvgRGBA(255, 0, 0, 200));
				nvgStrokeWidth(vg, 1.0);
				nvgStroke(vg);
			}

			nvgRestore(vg);
			stack_y += rms_height;
		}

		// plot scrolling spectrogram
		{
			// update image...
			int ss_index = 0;
			for (int y = 0; y < SS_HEIGHT; y++) {
				int si = SS_HEIGHT - 1 - y;
				for (int x = 0; x < SS_WIDTH; x++) {
					uint32_t color = 0xff000000;
					int rec_index = i1 - SS_WIDTH + x;
					if (rec_index >= 0) {
						float* spectrogram = Rec_get_float_data(&rec, REC_SPECTROGRAM, rec_index);
						float v = spectrogram[si] * 0.2;
						float r = v;
						float g = v*v*v*v;
						float b = v*v*v*v*v*v;
						//assert(isfinite(r));
						//assert(isfinite(g));
						//assert(isfinite(b));
						int ri = clampi(r * 255.0f, 0x0, 0xff);
						int gi = clampi(g * 255.0f, 0x0, 0xff);
						int bi = clampi(b * 255.0f, 0x0, 0xff);
						int ai = 0xff;
						color = (ri << 0) | (gi << 8) | (bi << 16) | (ai << 24);
					}
					ss_data[ss_index++] = color;
				}
			}
			nvgUpdateImage(vg, ss_handle, (const unsigned char*)ss_data);

			nvgSave(vg);
			nvgTranslate(vg, 0, stack_y);
			nvgBeginPath(vg);
			nvgRect(vg, 0, 0, SS_WIDTH, SS_HEIGHT);
			nvgFillPaint(vg, ss_paint);
			nvgFill(vg);
			nvgRestore(vg);
			stack_y += SS_HEIGHT;
		}

		// plot pitch track
		nvgSave(vg);
		nvgTranslate(vg, 0, screen_height);
		for (int track = 0; track < 2; track++) {
			if (track == 1 && mode != 2) break;
			if (track == 0 && !show_f0) continue;
			if (track == 1 && !show_expected_f0) continue;
			for (int pass = 0; pass < 2; pass++) {
				float x = x0;
				int prev_voiced = 0;
				nvgBeginPath(vg);
				for (int i = i0; i < i1; i++) {
					struct RecSimple* s = &rec.simple[i];
					float f0 = track == 0 ? s->f0 : s->expected_f0;
					int voiced = f0 > 0.0;

					if (voiced) {
						int attack = voiced && !prev_voiced;
						float pitch = f2p(f0) + 12.0*3.0;
						float y = -pitch * 8.0f;

						const int y_max = -10;
						if (y > y_max) y = y_max;

						if (attack) {
							if (pass == 1) nvgCircle(vg, x, y, 5);
							if (pass == 0) nvgMoveTo(vg, x, y);
						} else if (voiced) {
							if (pass == 0) nvgLineTo(vg, x, y);
						}
					}

					prev_voiced = voiced;
					x += xstep;
				}
				if (pass == 0) {
					const int a = 255;
					nvgStrokeColor(vg, track == 0 ? nvgRGBA(0, 255, 0, a) : nvgRGBA(255, 0, 0, a));
					nvgStrokeWidth(vg, 1.5);
					nvgStroke(vg);
				} else {
					const int a = 50;
					nvgFillColor(vg, track == 0 ? nvgRGBA(0, 255, 0, a) : nvgRGBA(255, 0, 0, a));
					nvgFill(vg);
				}
			}

			nvgBeginPath(vg);
			int prev_voiced = 0;
			int voiced_start = 0;
			for (int i = i0; i < i1; i++) {
				struct RecSimple* s = &rec.simple[i];
				float f0 = track == 0 ? s->f0 : s->expected_f0;
				int voiced = f0 > 0.0;

				int voiced_next = 0;
				if (i+1 < i1) {
					struct RecSimple* s_next = &rec.simple[i+1];
					float f0_next = track == 0 ? s_next->f0 : s_next->expected_f0;
					voiced_next = f0_next > 0.0;
				}

				if (voiced && !prev_voiced) voiced_start = i;
				if (voiced && !voiced_next) {
					const float h = screen_height - stack_y;
					nvgRect(vg, x0 + voiced_start * xstep, -h, (i-voiced_start) * xstep, h);
				}

				prev_voiced = voiced;
			}
			const int a = 25;
			nvgFillColor(vg, track == 0 ? nvgRGBA(0, 255, 0, a) : nvgRGBA(255, 0, 0, a));
			nvgFill(vg);
		}
		nvgRestore(vg);

		// info text
		{
			const float font_size = 18;
			nvgFontSize(vg, font_size);
			nvgTextAlign(vg, NVG_ALIGN_LEFT);
			nvgFontFaceId(vg, monofont);
			char buf[1000];
			char* mode_str;
			switch (mode) {
			case 0: mode_str = "REC"; break;
			case 1: mode_str = "CROP"; break;
			case 2: mode_str = "CASE"; break;
			default:
				mode_str = "???";
				break;
			}

			struct RecSimple* simple = &rec.simple[windows_index];
			float rms = simple->rms;
			if (mode < 2) {
				stbsp_snprintf(buf, sizeof buf, "MODE=%s rms=%.6f floor=%.6f", mode_str, rms, pp.noise_gate_rms_threshold);
			} else {
				float f0 = simple->f0;
				//float expected_f0 = simple->expected_f0;
				stbsp_snprintf(
					buf,
					sizeof buf,
					"MODE=%s frame=%d rms=%.6f floor=%.6f f0=%.1f (FFT:NSDF=%.2f)",
					mode_str, case_frame, rms, pp.noise_gate_rms_threshold, f0, simple->spectrogram_corrected_best_frequency/simple->nsdf_best_frequency);

			}
			nvgFillColor(vg, nvgRGBA(200,255,255,255));
			nvgText(vg, 5, 22, buf, NULL);
		}

		// tuner
		{
			float f0 = rec.simple[windows_index].f0;
			const float font_size = 70;
			nvgFontSize(vg, font_size);
			nvgTextAlign(vg, NVG_ALIGN_CENTER);
			nvgFontFaceId(vg, monofont);
			char buf[1000];
			float dc;
			if (f0 > 0) {
				float pitch = f2p(f0);
				float note = pitch + 0.5 - 3.0 + 12.0 * 4;
				dc = (fmodf(note, 1.0) - 0.5) * 2.0f;
				int octave = (int)note/12;
				if (octave < 0) octave = 0;
				if (octave > 9) octave = 9;
				int rnote = ((int)note) % 12;
				if (rnote < 0) rnote = 0; // mod be weird yo
				const char* notes = "C-C#D-D#E-F-F#G-G#A-A#B-";
				stbsp_snprintf(buf, sizeof buf, "%c%c%c", notes[rnote*2], notes[rnote*2+1], '0'+octave);
			} else {
				stbsp_snprintf(buf, sizeof buf, "---");
			}

			nvgFillColor(vg, nvgRGBA(255,255,50,255));
			nvgText(vg, screen_width - 90, 80, buf, NULL);

			if (f0 > 0) {
				const float other_font_size = 20;
				nvgFontSize(vg, other_font_size);
				stbsp_snprintf(buf, sizeof buf, "%.1fhz", f0);
				nvgFillColor(vg, nvgRGBA(255,50,50,255));
				nvgText(vg, screen_width - 90, 130, buf, NULL);
			}

			nvgBeginPath(vg);
			float x = screen_width - 90;
			float y = 100;
			float w = 75;
			float h = 10;
			nvgMoveTo(vg, x-w, y);
			nvgLineTo(vg, x+w, y);
			nvgMoveTo(vg, x, y-h);
			nvgLineTo(vg, x, y+h);
			nvgStrokeColor(vg, nvgRGBA(255, 255, 255, f0 > 0 ? 150 : 50));
			nvgStrokeWidth(vg, 1.0);
			nvgStroke(vg);

			if (f0 > 0) {
				nvgBeginPath(vg);
				nvgCircle(vg, x+dc*w, y, 5);
				nvgFillColor(vg, nvgRGBA(255, 255, 0, 255));
				nvgFill(vg);
			}
		}

		// cursor/selection stuff
		float sy = scroller_begin_y;
		float sh = screen_height - scroller_begin_y;
		if (mode == 1) {

			if (mouse_index >= 0 && mouse_index < i1) {
				float w = xstep;
				float x = x0 + mouse_index * xstep;
				nvgBeginPath(vg);
				nvgRect(vg, x, sy, w, sh);
				nvgFillColor(vg, (frame & 16) ? nvgRGBA(255,255,0,200) : nvgRGBA(255,0,0,200));
				nvgFill(vg);
			}

			nvgBeginPath(vg);
			nvgRect(vg, x0 + crop_i0_norm * xstep, sy, (crop_i1_norm-crop_i0_norm) * xstep, sh);
			nvgFillColor(vg, (frame & 8) ? nvgRGBA(50,100,255,70) : nvgRGBA(50,100,255,50));
			nvgFill(vg);
		} else if (mode == 2) {
			float w = xstep;
			float x = x0 + case_frame * xstep;
			nvgBeginPath(vg);
			nvgRect(vg, x, sy, w, sh);
			nvgFillColor(vg, (frame & 8) ? nvgRGBA(255,255,0,150) : nvgRGBA(255,0,0,200));
			nvgFill(vg);
		}

		SDL_LockMutex(playback_lock);
		int copy_playback_start = playback_start;
		int copy_playback_position_frame = playback_position_frame;
		SDL_UnlockMutex(playback_lock);
		if (copy_playback_start) {
			float x = x0 + copy_playback_position_frame * xstep;
			float w = xstep;
			nvgBeginPath(vg);
			nvgRect(vg, x, sy, w, sh);
			nvgFillColor(vg, nvgRGBA(255,255,255,200));
			nvgFill(vg);
		}

		nvgEndFrame(vg);
		SDL_GL_SwapWindow(window);
		frame++;
	}

	video_exit();

	return EXIT_SUCCESS;
}
