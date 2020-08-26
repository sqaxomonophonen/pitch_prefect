#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <time.h>

#include <SDL.h>

#include <soundio/soundio.h>

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

#define STB_SPRINTF_IMPLEMENTATION
#include "stb_sprintf.h"

#define PP_IMPLEMENTATION
#include "pitch_prefect.h"

#include "sound.h"

struct PP pp;


float pp_output_f0;
SDL_mutex* pp_output_lock;

static inline float f2p(float f)
{
	return log2f(f / 440.0f) * 12.0f;
}

static inline float p2f(float p)
{
	const float ti = 1.0f / 12.0f;
	return exp2f(p * ti) * 440.0f;
}

static inline void pitch_name(char* s, float p)
{
	int pi = roundf(p) + 12*4;
	int note = pi%12;
	int octave = pi/12;
	if (octave < 0) octave = 0;
	if (octave > 9) octave = 9;
	const char* nn = "A-A#B-C-C#D-D#E-F-F#G-G#";
	s[0] = nn[note*2];
	s[1] = nn[note*2+1];
	s[2] = '0' + octave;
}

static void audio_in_callback(struct SoundIoInStream *instream, int frame_count_min, int frame_count_max)
{
	#if 0
	double latency;
	int e = soundio_instream_get_latency(instream, &latency);
	printf("instream latency: %f (%d)\n", latency, e);
	#endif

	int frames_left = frame_count_max;

	int n_read = 0;
	while (frames_left > 0) {
		struct SoundIoChannelArea *areas;
		int frame_count = frames_left;
		int err;
		if ((err = soundio_instream_begin_read(instream, &areas, &frame_count))) {
			fprintf(stderr, "AUDIO ERROR begin read: %s", soundio_strerror(err));
			exit(EXIT_FAILURE);
		}
		if (frame_count == 0) break;

		if (areas) {
			for (int frame = 0; frame < frame_count; frame++) {
				const int ch = 0;
				float x = *((float*)areas[ch].ptr);
				if (pp_write_one(&pp, x)) {
					SDL_LockMutex(pp_output_lock);
					pp_output_f0 = pp.f0;
					SDL_UnlockMutex(pp_output_lock);
				}
				areas[ch].ptr += areas[ch].step;
			}
		}

		if ((err = soundio_instream_end_read(instream))) {
			fprintf(stderr, "AUDIO ERROR end read: %s", soundio_strerror(err));
			exit(EXIT_FAILURE);
		}

		n_read += frame_count;
		frames_left -= frame_count;
		if (frames_left <= 0) break;
	}
}

#define N_VOICES (4)

struct Voice {
	struct Op1 op1;
	int age;
};

struct Voice voices[N_VOICES];
SDL_mutex* voices_lock;

static void voices_init()
{
	for (int i = 0; i < N_VOICES; i++) {
		struct Voice* voice = &voices[i];
		Op1_init(&voice->op1);
		voice->op1.numer = 4;
		voice->op1.denom = 1;
	}
}

static inline void note_on(float pitch)
{
	SDL_LockMutex(voices_lock);
	int best_voice = 0;
	int best_voice_age = 0;
	for (int i = 0; i < N_VOICES; i++) {
		int age = voices[i].age;
		if (i == 0 || age > best_voice_age) {
			best_voice = i;
			best_voice_age = age;
		}
	}

	struct Voice* voice = &voices[best_voice];
	voice->age = 0;
	Op1_pitch(&voice->op1, pitch);
	SDL_UnlockMutex(voices_lock);
}

#define PI (3.141592653589793)
#define PI2 (PI*2.0)

static void audio_out_callback(struct SoundIoOutStream *outstream, int frame_count_min, int frame_count_max)
{
	#if 0
	double latency;
	int e = soundio_outstream_get_latency(outstream, &latency);
	printf("outstream latency: %f (%d)\n", latency, e);
	#endif
	// TODO

	const struct SoundIoChannelLayout *layout = &outstream->layout;
	struct SoundIoChannelArea *areas;
	int frames_left = frame_count_max;
	int err;

	while (frames_left > 0) {
		int frame_count = frames_left;
		if ((err = soundio_outstream_begin_write(outstream, &areas, &frame_count))) {
			fprintf(stderr, "%s\n", soundio_strerror(err));
			abort();
		}

		SDL_LockMutex(voices_lock);
		for (int frame = 0; frame < frame_count; frame++) {
			float sample = 0.0;
			for (int i = 0; i < N_VOICES; i++) {
				struct Voice* voice = &voices[i];
				int mod = 80 - (voice->age / 100);
				if (mod < 0) mod = 0;
				voice->op1.mod = mod;
				float ampl = 1.0f - (voice->age * 3e-5);
				if (ampl > 0) {
					sample += Op1_sample(&voice->op1) * ampl;
				}
				voice->age++;
			}
			sample *= 0.1f;
			for (int channel = 0; channel < layout->channel_count; channel += 1) {
				float *ptr = (float*)(areas[channel].ptr + areas[channel].step * frame);
				*ptr = sample;
			}
		}
		SDL_UnlockMutex(voices_lock);

		if ((err = soundio_outstream_end_write(outstream))) {
			fprintf(stderr, "%s\n", soundio_strerror(err));
			abort();
		}

		frames_left -= frame_count;
	}
}

static void audio_in_overflow_callback(struct SoundIoInStream *instream) {
	fprintf(stderr, "AUDIO-IN OVERRUN\n");
}

static void audio_in_error_callback(struct SoundIoInStream *instream, int err) {
	fprintf(stderr, "AUDIO-IN ERROR %d\n", err);
}

static void audio_out_error_callback(struct SoundIoOutStream *outstream, int err) {
	fprintf(stderr, "AUDIO-OUT ERROR %d\n", err);
}

struct SoundIo *sio;
static void audio_init()
{
	assert((sio = soundio_create()));

	int err;
	err = soundio_connect(sio);
	if (err) {
		fprintf(stderr, "error connecting: %s", soundio_strerror(err));
		exit(EXIT_FAILURE);
	}

	soundio_flush_events(sio);

	const int sample_rate = 48000;
	const float software_latency = 1.0 / 100.0;

	struct SoundIoInStream *instream;
	{
		int device_index = soundio_default_input_device_index(sio);
		struct SoundIoDevice *selected_device = soundio_get_input_device(sio, device_index);

		fprintf(stderr, "Input device: %s\n", selected_device->name);

		if (selected_device->probe_error) {
			fprintf(stderr, "Unable to probe device: %s\n", soundio_strerror(selected_device->probe_error));
			exit(EXIT_FAILURE);
		}

		soundio_device_sort_channel_layouts(selected_device);

		if (!soundio_device_supports_sample_rate(selected_device, sample_rate)) {
			fprintf(stderr, "Could not set sample rate to %d\n", sample_rate);
			exit(EXIT_FAILURE);
		}

		enum SoundIoFormat fmt = SoundIoFormatFloat32NE;

		instream = soundio_instream_create(selected_device);
		assert(instream);
		instream->format = fmt;
		instream->sample_rate = sample_rate;
		instream->software_latency = software_latency;
		instream->read_callback = audio_in_callback;
		instream->overflow_callback = audio_in_overflow_callback;
		instream->error_callback = audio_in_error_callback;
		instream->userdata = NULL;

		if ((err = soundio_instream_open(instream))) {
			fprintf(stderr, "unable to open input stream: %s", soundio_strerror(err));
			exit(EXIT_FAILURE);
		}

		//pp_init(&pp, instream->sample_rate);

	}

	struct SoundIoOutStream *outstream;
	{
		int device_index = soundio_default_output_device_index(sio);
		struct SoundIoDevice *selected_device = soundio_get_output_device(sio, device_index);

		fprintf(stderr, "Output device: %s\n", selected_device->name);

		if (selected_device->probe_error) {
			fprintf(stderr, "Unable to probe device: %s\n", soundio_strerror(selected_device->probe_error));
			exit(EXIT_FAILURE);
		}

		soundio_device_sort_channel_layouts(selected_device);

		const int sample_rate = 48000;
		if (!soundio_device_supports_sample_rate(selected_device, sample_rate)) {
			fprintf(stderr, "Could not set sample rate to %d\n", sample_rate);
			exit(EXIT_FAILURE);
		}

		enum SoundIoFormat fmt = SoundIoFormatFloat32NE;

		outstream = soundio_outstream_create(selected_device);
		assert(outstream);
		outstream->format = fmt;
		outstream->sample_rate = sample_rate;
		outstream->software_latency = software_latency * 10;
		outstream->write_callback = audio_out_callback;
		//outstream->overflow_callback = overflow_callback;
		outstream->error_callback = audio_out_error_callback;
		outstream->userdata = NULL;

		if ((err = soundio_outstream_open(outstream))) {
			fprintf(stderr, "unable to open input stream: %s", soundio_strerror(err));
			exit(EXIT_FAILURE);
		}

	}

	printf("instream->sample_rate = %d\n", instream->sample_rate);
	printf("instream->software_latency = %f\n", instream->software_latency);
	printf("outstream->sample_rate = %d\n", outstream->sample_rate);
	printf("outstream->software_latency = %f\n", outstream->software_latency);
	//outstream_sample_rate = outstream->sample_rate;

	set_sample_rate(outstream->sample_rate);
	voices_init();

	pp_init(&pp, instream->sample_rate);
	pp.noise_gate_rms_threshold = 0;

	if ((err = soundio_instream_start(instream))) {
		fprintf(stderr, "unable to start input device: %s", soundio_strerror(err));
		exit(EXIT_FAILURE);
	}

	if ((err = soundio_outstream_start(outstream))) {
		fprintf(stderr, "unable to start output device: %s", soundio_strerror(err));
		exit(EXIT_FAILURE);
	}
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
			"fun game",
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

char* prg;

static inline float lerp(float t, float a, float b)
{
	return a + (b-a) * t;
}

#define G1_N_NOTES (3)

struct note {
	float pitch;
	int done;
	float score;
};

static int note_compar(const void* va, const void* vb)
{
	const struct note* a = va;
	const struct note* b = vb;
	float pa = a->pitch;
	float pb = b->pitch;
	if (pa < pb) return -1;
	if (pa > pb) return 1;
	return 0;
}

struct {
	struct note notes[G1_N_NOTES];
} g1;

static void game1_frame(int frame, float f0)
{
	if (frame == 0) {
		memset(&g1, 0, sizeof g1);
		for (int i = 0; i < G1_N_NOTES; i++) {
			g1.notes[i].pitch = -24 + (rand() % 25);
		}
		qsort(g1.notes, G1_N_NOTES, sizeof g1.notes[0], note_compar);
	}

	int voiced = f0 > 0;
	float pitch = voiced ? f2p(f0) : 0;

	nvgGlobalCompositeOperation(vg, NVG_LIGHTER); // additive color blending
	nvgSave(vg);

	const int padding = 100;
	const float w = (float)(screen_width - (G1_N_NOTES+1) * padding) / (float)G1_N_NOTES;
	const float y = padding;
	const float h = screen_height - padding*2;

	const float font_size = 40;
	nvgFontSize(vg, font_size);
	nvgTextAlign(vg, NVG_ALIGN_CENTER);
	nvgFontFaceId(vg, monofont);
	char buf[1000];

	for (int i = 0; i < G1_N_NOTES; i++) {
		struct note* n = &g1.notes[i];
		float x = (padding + w) * i + padding;
		nvgBeginPath(vg);
		nvgRect(vg, x, y, w, h);
		int a = 200;
		nvgFillColor(vg, n->done ? nvgRGBA(0, 255, 0, a) : nvgRGBA(255, 0, 0, a));
		nvgFill(vg);

		n->score -= 0.01f;
		if (n->score < 0) n->score = 0;

		if (!n->done) {
			float mx = x + w*0.5;
			float my = y + h*0.5;
			float r = 40;
			nvgBeginPath(vg);
			if (!voiced) {
				nvgCircle(vg, mx, my, r);
			} else {
				if (pitch < n->pitch) {
					nvgMoveTo(vg, mx, my-r);
					nvgLineTo(vg, mx+r, my+r);
					nvgLineTo(vg, mx-r, my+r);
					nvgClosePath(vg);
				} else if (pitch > n->pitch) {
					nvgMoveTo(vg, mx, my+r);
					nvgLineTo(vg, mx-r, my-r);
					nvgLineTo(vg, mx+r, my-r);
					nvgClosePath(vg);
				}

				float dp = fabsf(pitch - n->pitch);
				if (dp < 0.4) {
					n->score += 0.02f;
					if (n->score > 1.0f) {
						n->done = 1;
					}
				}
			}
			nvgFillColor(vg, nvgRGBA(100, 100, 100, a));
			nvgFill(vg);

			if (n->score > 0) {
				nvgBeginPath(vg);
				nvgRect(vg, x, y+h-20, w*n->score, 20);
				nvgFillColor(vg, nvgRGBA(255, 255, 255, 255));
				nvgFill(vg);
			}

			//stbsp_snprintf(buf, sizeof buf, "%.0f", n->pitch+24);
			//pitch_name(buf, n->pitch);
			//buf[3] = 0;

			char pn[3];
			pitch_name(pn, n->pitch);
			stbsp_snprintf(buf, sizeof buf, "%c%c%c (%.0f)", pn[0], pn[1], pn[2], n->pitch+24);

			nvgFillColor(vg, nvgRGBA(255,255,255,100));
			nvgText(vg, mx, my + h/3, buf, NULL);
		}
	}
	
	if (voiced) {
		nvgTextAlign(vg, NVG_ALIGN_LEFT);
		char pn[3];
		pitch_name(pn, pitch);
		stbsp_snprintf(buf, sizeof buf, "%c%c%c (%.0f)", pn[0], pn[1], pn[2], pitch+24);
		nvgFillColor(vg, nvgRGBA(255,255,255,100));
		nvgText(vg, 10, 40, buf, NULL);
	}

	nvgRestore(vg);
}

#define G2_N_NOTES 3
struct {
	struct note notes[G2_N_NOTES];
	int counter;
} g2;

static void game2_frame(int frame, float f0)
{
	int new_puzzle = 0;
	if (frame == 0) {
		// init
		new_puzzle = 1;
	}

	int all_done = 1;
	for (int i = 0; i < G2_N_NOTES; i++) {
		if (!g2.notes[i].done) {
			all_done = 0;
			break;
		}
	}

	if (new_puzzle || all_done) {
		memset(&g2, 0, sizeof g2);
		for (int i = 0; i < G2_N_NOTES; i++) {
			g2.notes[i].pitch = -28 + (rand() % 20);
		}
		qsort(g2.notes, G2_N_NOTES, sizeof g2.notes[0], note_compar);
	}

	const int frames_per_note = 40;

	int voiced = f0 > 0;
	float pitch = voiced ? f2p(f0) : 0;

	int flash_index = 0;
	int flash_add = 0;

	int c = g2.counter;
	if (c < frames_per_note * G2_N_NOTES) {
		voiced = 0; // :-) don't activate while playing...
		if ((g2.counter % frames_per_note) == 0) {
			note_on(g2.notes[c / frames_per_note].pitch);
		}
		flash_index = c / frames_per_note;
		flash_add = (frames_per_note - (c % frames_per_note) - 1);
	}

	nvgGlobalCompositeOperation(vg, NVG_LIGHTER); // additive color blending
	nvgSave(vg);

	const int padding = 100;
	const float w = (float)(screen_width - (G2_N_NOTES+1) * padding) / (float)G2_N_NOTES;
	const float y = padding;
	const float h = screen_height - padding*2;

	const float font_size = 40;
	nvgFontSize(vg, font_size);
	nvgTextAlign(vg, NVG_ALIGN_CENTER);
	nvgFontFaceId(vg, monofont);
	char buf[1000];

	for (int i = 0; i < G2_N_NOTES; i++) {
		struct note* n = &g2.notes[i];
		float x = (padding + w) * i + padding;
		nvgBeginPath(vg);
		nvgRect(vg, x, y, w, h);
		int a = 200;
		int add = (i == flash_index ? flash_add : 0);
		nvgFillColor(vg, n->done ? nvgRGBA(0 + add, 255 + add, 0 + add, a + add) : nvgRGBA(255 + add, 0 + add, 0 + add, a + add));
		nvgFill(vg);

		n->score -= 0.01f;
		if (n->score < 0) n->score = 0;

		if (!n->done) {
			float mx = x + w*0.5;
			float my = y + h*0.5;
			float r = 40;
			nvgBeginPath(vg);
			if (!voiced) {
				nvgCircle(vg, mx, my, r);
			} else {
				if (pitch < n->pitch) {
					nvgMoveTo(vg, mx, my-r);
					nvgLineTo(vg, mx+r, my+r);
					nvgLineTo(vg, mx-r, my+r);
					nvgClosePath(vg);
				} else if (pitch > n->pitch) {
					nvgMoveTo(vg, mx, my+r);
					nvgLineTo(vg, mx-r, my-r);
					nvgLineTo(vg, mx+r, my-r);
					nvgClosePath(vg);
				}

				float dp = fabsf(pitch - n->pitch);
				if (dp < 0.4) {
					n->score += 0.04f;
					if (n->score > 1.0f) {
						n->done = 1;
					}
				}
			}
			nvgFillColor(vg, nvgRGBA(100, 100, 100, a));
			nvgFill(vg);

			if (n->score > 0) {
				nvgBeginPath(vg);
				nvgRect(vg, x, y+h-20, w*n->score, 20);
				nvgFillColor(vg, nvgRGBA(255, 255, 255, 255));
				nvgFill(vg);
			}

			//stbsp_snprintf(buf, sizeof buf, "%.0f", n->pitch+24);
			//pitch_name(buf, n->pitch);
			//buf[3] = 0;

			char pn[3];
			pitch_name(pn, n->pitch);
			stbsp_snprintf(buf, sizeof buf, "%c%c%c (%.0f)", pn[0], pn[1], pn[2], n->pitch+24);

			nvgFillColor(vg, nvgRGBA(255,255,255,100));
			nvgText(vg, mx, my + h/3, buf, NULL);
		}
	}
	
	if (voiced) {
		nvgTextAlign(vg, NVG_ALIGN_LEFT);
		char pn[3];
		pitch_name(pn, pitch);
		stbsp_snprintf(buf, sizeof buf, "%c%c%c (%.0f)", pn[0], pn[1], pn[2], pitch+24);
		nvgFillColor(vg, nvgRGBA(255,255,255,100));
		nvgText(vg, 10, 40, buf, NULL);
	}

	nvgRestore(vg);

	g2.counter++;
}

int main(int argc, char** argv)
{
	srand(time(NULL)); // randomize seed

	prg = argv[0];

	assert((pp_output_lock = SDL_CreateMutex()) != NULL);

	audio_init();

	assert(SDL_Init(SDL_INIT_TIMER | SDL_INIT_VIDEO) == 0);
	atexit(SDL_Quit);

	video_init();

	void(*game_frame)(int, float) = game2_frame;
	void(*set_game_frame)(int, float) = NULL;

	int exiting = 0;
	int fullscreen = 0;
	int frame = 0;
	while (!exiting) {
		SDL_Event e;
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
				} else if (sym == '1') {
					set_game_frame = game1_frame;
				} else if (sym == '2') {
					set_game_frame = game2_frame;
				}

			} else if (e.type == SDL_MOUSEMOTION) {
				//mx = e.motion.x;
				//my = e.motion.y;
			} else if (e.type == SDL_MOUSEBUTTONDOWN || e.type == SDL_MOUSEBUTTONUP) {
				//e.button;
			} else if (e.type == SDL_WINDOWEVENT) {
				if (e.window.event == SDL_WINDOWEVENT_RESIZED) {
					window_size(&screen_width, &screen_height, &pixel_ratio);
				}
			}
		}

		if (set_game_frame != NULL) {
			frame = 0;
			game_frame = set_game_frame;
			set_game_frame = NULL;
		}

		SDL_LockMutex(pp_output_lock);
		float f0 = pp_output_f0;
		SDL_UnlockMutex(pp_output_lock);

		glViewport(0, 0, screen_width, screen_height);
		glClearColor(0, 0.05, 0.15, 0);
		glClear(GL_COLOR_BUFFER_BIT | GL_STENCIL_BUFFER_BIT);

		glEnable(GL_BLEND);
		glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
		glEnable(GL_CULL_FACE);
		glDisable(GL_DEPTH_TEST);

		nvgBeginFrame(vg, screen_width / pixel_ratio, screen_height / pixel_ratio, pixel_ratio);

		game_frame(frame, f0);

		nvgEndFrame(vg);
		SDL_GL_SwapWindow(window);
		frame++;
	}

	video_exit();

	return EXIT_SUCCESS;
}
