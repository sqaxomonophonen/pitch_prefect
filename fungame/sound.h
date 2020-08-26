#include <math.h>

#define PI (3.141592653589793)
#define PI2 (PI*2.0)

#define OP_WTLEN (256)

#define BASE_FREQ (440.0)

#define F2P(f) (log2f(f / BASE_FREQ) * 12.0f)
#define P2F(p) (exp2f(p * (1.0f/12.0f)) * BASE_FREQ)

float sample_rate; /* FLW: one global sample rate is fine */
float sample_rate_inv; // derived


static void set_sample_rate(float r)
{
	sample_rate = r;
	sample_rate_inv = 1.0f / r;
}

struct Op {
	uint8_t* wavetable;
	uint32_t phase; // 8:16 fixed point
	uint32_t inc;
};


static inline void Op_take_wavetable(struct Op* op, uint8_t* wavetables, int index)
{
	op->wavetable = &wavetables[index * OP_WTLEN];
}

static inline void Op_freq(struct Op* op, float hz)
{
	op->inc = (uint32_t)roundf((float)(1<<24) * sample_rate_inv * hz);
}

static inline void Op_pitch(struct Op* op, float pitch)
{
	Op_freq(op, P2F(pitch));
}

static inline int Op_step(struct Op* op, int mod)
{
	uint32_t pha = op->phase;
	int idx0 = (pha >> 16) & 0xff;
	int idx1 = (idx0+1) & 0xff;
	uint8_t v0 = op->wavetable[idx0];
	uint8_t v1 = op->wavetable[idx1];
	int v = ((v0<<16) + (v1-v0) * (pha & 0xffff)) - (1<<23); // lerp; v0+(v1-v0)*t but in fixed-point
	op->phase = (op->phase + op->inc + mod) & ((1<<24)-1); // TODO try both here and as first statement (determines if mod affects now or delayed)
	return v;
}

static inline void Op_sin(struct Op* op)
{
	for (int i = 0; i < OP_WTLEN; i++) {
		op->wavetable[i] = (uint8_t)roundf(sinf(((float)i / (float)OP_WTLEN) * PI2) * 127.0f + 127.0f);
	}
}

static uint8_t* Op_wavetable_alloc(int n)
{
	// TODO 64-byte / L1 align?
	uint8_t* p = malloc(n * OP_WTLEN);
	assert(p != NULL);
	return p;
}

#define OP_SAMPLE(v) (((float)(v) / (float)(1<<24))-0.5f)
#define OP_MOD(x, y) (((x)*(y)) >> 12)

struct Op1 {
	struct Op op[2];
	int mod, numer, denom;
};

static inline void Op1_init(struct Op1* o)
{
	const int n = 2;
	memset(o, 0, sizeof *o);
	uint8_t* wavetables = Op_wavetable_alloc(n);
	for (int i = 0; i < n; i++) {
		Op_take_wavetable(&o->op[i], wavetables, i);
		Op_sin(&o->op[i]);
	}
}

static inline void Op1_pitch(struct Op1* o, float pitch)
{
	Op_pitch(&o->op[1], pitch);
	o->op[0].inc = (o->op[1].inc * o->numer) / o->denom;
}

static inline float Op1_sample(struct Op1* o)
{
	int x0 = Op_step(&o->op[0], 0);
	int x1 = Op_step(&o->op[1], OP_MOD(x0, o->mod));
	return OP_SAMPLE(x1);
}
