#ifndef PITCH_PREFECT_H

#ifdef PP_STATIC
#define PP_DEF static
#else
#define PP_DEF extern
#endif

#ifdef __cplusplus
extern "C" {
#endif

//const float voiced_threshold = 0.8f;

struct PP__FFT {
	int n;
	float* wsave;
	int* ifac;
};

struct PP {
	// config:
	float noise_gate_rms_threshold; // rms < noise_gate_rms_threshold will be unvoiced
	float pitch0_f0; // set to PP_PITCH0_F0 during pp_init() but you can change it online. if zero, pitch will not be outputted

	// output:
	float f0; // fundamental frequency, or 0 if unvoiced
	float pitch; // 12.0/oct, relative to pitch0_f0
	float rms; // root mean square

	int frame_counter;

	float* window;
	float* spectrogram;
	float* nsdf;

	float nsdf_best_frequency;
	float spectrogram_corrected_best_frequency;

	int decimation_factor;
	float decimated_sample_rate_hz;
	int decimation_counter;

	struct PP__FFT fft;
	float* fft_data;
	float* fft_window_function;

	// bucketry; input data is written here
	float* bucket_cursor;
	int bucket_remaining;

	int window_scroll_is_pending;
};

PP_DEF void pp_init(struct PP* pp, int input_sample_rate_hz);

PP_DEF void pp__flush(struct PP* pp);
PP_DEF void pp__scroll_window(struct PP* pp);

static inline int pp_write_one_decimated(struct PP* pp, float x)
{
	if (pp->window_scroll_is_pending) pp__scroll_window(pp);
	*pp->bucket_cursor = x;
	pp->bucket_cursor++;
	pp->bucket_remaining--;
	if (pp->bucket_remaining == 0) {
		pp__flush(pp);
		return 1;
	} else {
		return 0;
	}
}

static inline int pp_write_one(struct PP* pp, float x)
{
	if (pp->decimation_counter == 0) {
		pp->decimation_counter = pp->decimation_factor - 1;
		return pp_write_one_decimated(pp, x);
	} else {
		pp->decimation_counter--;
		return 0;
	}
}

PP_DEF void pp_reset(struct PP* pp);
PP_DEF void pp_destroy(struct PP* pp);

#ifdef __cplusplus
}
#endif

#define PITCH_PREFECT_H
#endif


//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////


#ifdef PP_IMPLEMENTATION

/* Pitch Prefect internally operates at a sample rate Rp:
 *   Rp=Ra/D
 * Ra is the sample rate of the data fed to Pitch Prefect (typically the audio
 * sample rate), and D is the "decimation factor"; a positive integer.
 * pp_init() will pick D so that Rp comes closest to PP_IDEAL_SAMPLE_RATE_HZ.
 * */
#ifndef PP_IDEAL_SAMPLE_RATE_HZ
#define PP_IDEAL_SAMPLE_RATE_HZ (12000)
#endif

/* Pitch Prefect generates a new f0 guess whenever the next bucket is filled.
 * */
#ifndef PP_BUCKET_SIZE
#define PP_BUCKET_SIZE (128)
#endif

/* Number of buckets determines the full NSDF/spectrogram window size (i.e.
 * PP_BUCKET_SIZE * PP_N_BUCKETS); buckets are shifted left whenever a new
 * bucket is filled. */
#ifndef PP_N_BUCKETS
#define PP_N_BUCKETS (4)
#endif

/* fundamental frequency of pitch=0 */
#ifndef PP_PITCH0_F0
#define PP_PITCH0_F0 (440.0f)
#endif

#define PP_WINDOW_SIZE (PP_BUCKET_SIZE * PP_N_BUCKETS)

/* FFTPACK port-to-C stolen from monty/xiph - public domain - thanks! see:
 *    http://www.netlib.org/fftpack/fft.c
 * The code here does forward real-FFT only; using it to make a spectrogram */

static void pp__drfti1(int n, float *wsave, int *ifac)
{
	static int ntryh[4] = { 4,2,3,5 };
	static float tpi = 6.28318530717958647692528676655900577;
	float arg,argh,argld,fi;
	int ntry=0,i,j=-1;
	int k1, l1, l2, ib;
	int ld, ii, ip, is, nq, nr;
	int ido, ipm, nfm1;
	int nl=n;
	int nf=0;

L101:
	j++;
	if (j < 4)
		ntry=ntryh[j];
	else
		ntry+=2;

L104:
	nq=nl/ntry;
	nr=nl-ntry*nq;
	if (nr!=0) goto L101;

	nf++;
	ifac[nf+1]=ntry;
	nl=nq;
	if(ntry!=2)goto L107;
	if(nf==1)goto L107;

	for (i=1;i<nf;i++){
		ib=nf-i+1;
		ifac[ib+1]=ifac[ib];
	}
	ifac[2] = 2;

L107:
	if(nl!=1)goto L104;
	ifac[0]=n;
	ifac[1]=nf;
	argh=tpi/n;
	is=0;
	nfm1=nf-1;
	l1=1;

	if(nfm1==0)return;

	for (k1=0;k1<nfm1;k1++){
		ip=ifac[k1+2];
		ld=0;
		l2=l1*ip;
		ido=n/l2;
		ipm=ip-1;

		for (j=0;j<ipm;j++){
			ld+=l1;
			i=is;
			argld=(float)ld*argh;
			fi=0.;
			for (ii=2;ii<ido;ii+=2){
				fi+=1.;
				arg=fi*argld;
				wsave[i++]=cosf(arg);
				wsave[i++]=sinf(arg);
			}
			is+=ido;
		}
		l1=l2;
	}
}

static void pp__dradf4(int ido,int l1,float *cc,float *ch,float *wa1, float *wa2,float *wa3)
{
	static float hsqt2 = .70710678118654752440084436210485;
	int i,k,t0,t1,t2,t3,t4,t5,t6;
	float ci2,ci3,ci4,cr2,cr3,cr4,ti1,ti2,ti3,ti4,tr1,tr2,tr3,tr4;
	t0=l1*ido;

	t1=t0;
	t4=t1<<1;
	t2=t1+(t1<<1);
	t3=0;

	for(k=0;k<l1;k++){
		tr1=cc[t1]+cc[t2];
		tr2=cc[t3]+cc[t4];
		ch[t5=t3<<2]=tr1+tr2;
		ch[(ido<<2)+t5-1]=tr2-tr1;
		ch[(t5+=(ido<<1))-1]=cc[t3]-cc[t4];
		ch[t5]=cc[t2]-cc[t1];

		t1+=ido;
		t2+=ido;
		t3+=ido;
		t4+=ido;
	}

	if(ido<2)return;
	if(ido==2)goto L105;

	t1=0;
	for(k=0;k<l1;k++){
		t2=t1;
		t4=t1<<2;
		t5=(t6=ido<<1)+t4;
		for(i=2;i<ido;i+=2){
			t3=(t2+=2);
			t4+=2;
			t5-=2;

			t3+=t0;
			cr2=wa1[i-2]*cc[t3-1]+wa1[i-1]*cc[t3];
			ci2=wa1[i-2]*cc[t3]-wa1[i-1]*cc[t3-1];
			t3+=t0;
			cr3=wa2[i-2]*cc[t3-1]+wa2[i-1]*cc[t3];
			ci3=wa2[i-2]*cc[t3]-wa2[i-1]*cc[t3-1];
			t3+=t0;
			cr4=wa3[i-2]*cc[t3-1]+wa3[i-1]*cc[t3];
			ci4=wa3[i-2]*cc[t3]-wa3[i-1]*cc[t3-1];

			tr1=cr2+cr4;
			tr4=cr4-cr2;
			ti1=ci2+ci4;
			ti4=ci2-ci4;
			ti2=cc[t2]+ci3;
			ti3=cc[t2]-ci3;
			tr2=cc[t2-1]+cr3;
			tr3=cc[t2-1]-cr3;


			ch[t4-1]=tr1+tr2;
			ch[t4]=ti1+ti2;

			ch[t5-1]=tr3-ti4;
			ch[t5]=tr4-ti3;

			ch[t4+t6-1]=ti4+tr3;
			ch[t4+t6]=tr4+ti3;

			ch[t5+t6-1]=tr2-tr1;
			ch[t5+t6]=ti1-ti2;
		}
		t1+=ido;
	}
	if(ido%2==1)return;

L105:

	t2=(t1=t0+ido-1)+(t0<<1);
	t3=ido<<2;
	t4=ido;
	t5=ido<<1;
	t6=ido;

	for(k=0;k<l1;k++){
		ti1=-hsqt2*(cc[t1]+cc[t2]);
		tr1=hsqt2*(cc[t1]-cc[t2]);
		ch[t4-1]=tr1+cc[t6-1];
		ch[t4+t5-1]=cc[t6-1]-tr1;
		ch[t4]=ti1-cc[t1+t0];
		ch[t4+t5]=ti1+cc[t1+t0];
		t1+=ido;
		t2+=ido;
		t4+=t3;
		t6+=ido;
	}
}

static void pp__dradf2(int ido,int l1,float *cc,float *ch,float *wa1)
{
	int i,k;
	float ti2,tr2;
	int t0,t1,t2,t3,t4,t5,t6;

	t1=0;
	t0=(t2=l1*ido);
	t3=ido<<1;
	for(k=0;k<l1;k++){
		ch[t1<<1]=cc[t1]+cc[t2];
		ch[(t1<<1)+t3-1]=cc[t1]-cc[t2];
		t1+=ido;
		t2+=ido;
	}

	if(ido<2)return;
	if(ido==2)goto L105;

	t1=0;
	t2=t0;
	for(k=0;k<l1;k++){
		t3=t2;
		t4=(t1<<1)+(ido<<1);
		t5=t1;
		t6=t1+t1;
		for(i=2;i<ido;i+=2){
			t3+=2;
			t4-=2;
			t5+=2;
			t6+=2;
			tr2=wa1[i-2]*cc[t3-1]+wa1[i-1]*cc[t3];
			ti2=wa1[i-2]*cc[t3]-wa1[i-1]*cc[t3-1];
			ch[t6]=cc[t5]+ti2;
			ch[t4]=ti2-cc[t5];
			ch[t6-1]=cc[t5-1]+tr2;
			ch[t4-1]=cc[t5-1]-tr2;
		}
		t1+=ido;
		t2+=ido;
	}

	if(ido%2==1)return;

L105:
	t3=(t2=(t1=ido)-1);
	t2+=t0;
	for(k=0;k<l1;k++){
		ch[t1]=-cc[t2];
		ch[t1-1]=cc[t3];
		t1+=ido<<1;
		t2+=ido;
		t3+=ido;
	}
}

static void pp__dradfg(int ido,int ip,int l1,int idl1,float *cc,float *c1, float *c2,float *ch,float *ch2,float *wa)
{
	static float tpi=6.28318530717958647692528676655900577;
	int idij,ipph,i,j,k,l,ic,ik,is;
	int t0,t1,t2,t3,t4,t5,t6,t7,t8,t9,t10;
	float dc2,ai1,ai2,ar1,ar2,ds2;
	int nbd;
	float dcp,arg,dsp,ar1h,ar2h;
	int idp2,ipp2;

	arg=tpi/(float)ip;
	dcp=cosf(arg);
	dsp=sinf(arg);
	ipph=(ip+1)>>1;
	ipp2=ip;
	idp2=ido;
	nbd=(ido-1)>>1;
	t0=l1*ido;
	t10=ip*ido;

	if(ido==1)goto L119;
	for(ik=0;ik<idl1;ik++)ch2[ik]=c2[ik];

	t1=0;
	for(j=1;j<ip;j++){
		t1+=t0;
		t2=t1;
		for(k=0;k<l1;k++){
			ch[t2]=c1[t2];
			t2+=ido;
		}
	}

	is=-ido;
	t1=0;
	if(nbd>l1){
		for(j=1;j<ip;j++){
			t1+=t0;
			is+=ido;
			t2= -ido+t1;
			for(k=0;k<l1;k++){
				idij=is-1;
				t2+=ido;
				t3=t2;
				for(i=2;i<ido;i+=2){
					idij+=2;
					t3+=2;
					ch[t3-1]=wa[idij-1]*c1[t3-1]+wa[idij]*c1[t3];
					ch[t3]=wa[idij-1]*c1[t3]-wa[idij]*c1[t3-1];
				}
			}
		}
	}else{

		for(j=1;j<ip;j++){
			is+=ido;
			idij=is-1;
			t1+=t0;
			t2=t1;
			for(i=2;i<ido;i+=2){
				idij+=2;
				t2+=2;
				t3=t2;
				for(k=0;k<l1;k++){
					ch[t3-1]=wa[idij-1]*c1[t3-1]+wa[idij]*c1[t3];
					ch[t3]=wa[idij-1]*c1[t3]-wa[idij]*c1[t3-1];
					t3+=ido;
				}
			}
		}
	}

	t1=0;
	t2=ipp2*t0;
	if(nbd<l1){
		for(j=1;j<ipph;j++){
			t1+=t0;
			t2-=t0;
			t3=t1;
			t4=t2;
			for(i=2;i<ido;i+=2){
				t3+=2;
				t4+=2;
				t5=t3-ido;
				t6=t4-ido;
				for(k=0;k<l1;k++){
					t5+=ido;
					t6+=ido;
					c1[t5-1]=ch[t5-1]+ch[t6-1];
					c1[t6-1]=ch[t5]-ch[t6];
					c1[t5]=ch[t5]+ch[t6];
					c1[t6]=ch[t6-1]-ch[t5-1];
				}
			}
		}
	}else{
		for(j=1;j<ipph;j++){
			t1+=t0;
			t2-=t0;
			t3=t1;
			t4=t2;
			for(k=0;k<l1;k++){
				t5=t3;
				t6=t4;
				for(i=2;i<ido;i+=2){
					t5+=2;
					t6+=2;
					c1[t5-1]=ch[t5-1]+ch[t6-1];
					c1[t6-1]=ch[t5]-ch[t6];
					c1[t5]=ch[t5]+ch[t6];
					c1[t6]=ch[t6-1]-ch[t5-1];
				}
				t3+=ido;
				t4+=ido;
			}
		}
	}

L119:
	for(ik=0;ik<idl1;ik++)c2[ik]=ch2[ik];

	t1=0;
	t2=ipp2*idl1;
	for(j=1;j<ipph;j++){
		t1+=t0;
		t2-=t0;
		t3=t1-ido;
		t4=t2-ido;
		for(k=0;k<l1;k++){
			t3+=ido;
			t4+=ido;
			c1[t3]=ch[t3]+ch[t4];
			c1[t4]=ch[t4]-ch[t3];
		}
	}

	ar1=1.;
	ai1=0.;
	t1=0;
	t2=ipp2*idl1;
	t3=(ip-1)*idl1;
	for(l=1;l<ipph;l++){
		t1+=idl1;
		t2-=idl1;
		ar1h=dcp*ar1-dsp*ai1;
		ai1=dcp*ai1+dsp*ar1;
		ar1=ar1h;
		t4=t1;
		t5=t2;
		t6=t3;
		t7=idl1;

		for(ik=0;ik<idl1;ik++){
			ch2[t4++]=c2[ik]+ar1*c2[t7++];
			ch2[t5++]=ai1*c2[t6++];
		}

		dc2=ar1;
		ds2=ai1;
		ar2=ar1;
		ai2=ai1;

		t4=idl1;
		t5=(ipp2-1)*idl1;
		for(j=2;j<ipph;j++){
			t4+=idl1;
			t5-=idl1;

			ar2h=dc2*ar2-ds2*ai2;
			ai2=dc2*ai2+ds2*ar2;
			ar2=ar2h;

			t6=t1;
			t7=t2;
			t8=t4;
			t9=t5;
			for(ik=0;ik<idl1;ik++){
				ch2[t6++]+=ar2*c2[t8++];
				ch2[t7++]+=ai2*c2[t9++];
			}
		}
	}

	t1=0;
	for(j=1;j<ipph;j++){
		t1+=idl1;
		t2=t1;
		for(ik=0;ik<idl1;ik++)ch2[ik]+=c2[t2++];
	}

	if(ido<l1)goto L132;

	t1=0;
	t2=0;
	for(k=0;k<l1;k++){
		t3=t1;
		t4=t2;
		for(i=0;i<ido;i++)cc[t4++]=ch[t3++];
		t1+=ido;
		t2+=t10;
	}

	goto L135;

L132:
	for(i=0;i<ido;i++){
		t1=i;
		t2=i;
		for(k=0;k<l1;k++){
			cc[t2]=ch[t1];
			t1+=ido;
			t2+=t10;
		}
	}

L135:
	t1=0;
	t2=ido<<1;
	t3=0;
	t4=ipp2*t0;
	for(j=1;j<ipph;j++){

		t1+=t2;
		t3+=t0;
		t4-=t0;

		t5=t1;
		t6=t3;
		t7=t4;

		for(k=0;k<l1;k++){
			cc[t5-1]=ch[t6];
			cc[t5]=ch[t7];
			t5+=t10;
			t6+=ido;
			t7+=ido;
		}
	}

	if(ido==1)return;
	if(nbd<l1)goto L141;

	t1=-ido;
	t3=0;
	t4=0;
	t5=ipp2*t0;
	for(j=1;j<ipph;j++){
		t1+=t2;
		t3+=t2;
		t4+=t0;
		t5-=t0;
		t6=t1;
		t7=t3;
		t8=t4;
		t9=t5;
		for(k=0;k<l1;k++){
			for(i=2;i<ido;i+=2){
				ic=idp2-i;
				cc[i+t7-1]=ch[i+t8-1]+ch[i+t9-1];
				cc[ic+t6-1]=ch[i+t8-1]-ch[i+t9-1];
				cc[i+t7]=ch[i+t8]+ch[i+t9];
				cc[ic+t6]=ch[i+t9]-ch[i+t8];
			}
			t6+=t10;
			t7+=t10;
			t8+=ido;
			t9+=ido;
		}
	}
	return;

L141:

	t1=-ido;
	t3=0;
	t4=0;
	t5=ipp2*t0;
	for(j=1;j<ipph;j++){
		t1+=t2;
		t3+=t2;
		t4+=t0;
		t5-=t0;
		for(i=2;i<ido;i+=2){
			t6=idp2+t1-i;
			t7=i+t3;
			t8=i+t4;
			t9=i+t5;
			for(k=0;k<l1;k++){
				cc[t7-1]=ch[t8-1]+ch[t9-1];
				cc[t6-1]=ch[t8-1]-ch[t9-1];
				cc[t7]=ch[t8]+ch[t9];
				cc[t6]=ch[t9]-ch[t8];
				t6+=t10;
				t7+=t10;
				t8+=ido;
				t9+=ido;
			}
		}
	}
}

static void pp__drftf1(int n,float *c,float *ch,float *wa,int *ifac)
{
	int i,k1,l1,l2;
	int na,kh,nf;
	int ip,iw,ido,idl1,ix2,ix3;

	nf=ifac[1];
	na=1;
	l2=n;
	iw=n;

	for (k1=0;k1<nf;k1++) {
		kh=nf-k1;
		ip=ifac[kh+1];
		l1=l2/ip;
		ido=n/l2;
		idl1=ido*l1;
		iw-=(ip-1)*ido;
		na=1-na;

		if(ip!=4)goto L102;

		ix2=iw+ido;
		ix3=ix2+ido;
		if(na!=0)
			pp__dradf4(ido,l1,ch,c,wa+iw-1,wa+ix2-1,wa+ix3-1);
		else
			pp__dradf4(ido,l1,c,ch,wa+iw-1,wa+ix2-1,wa+ix3-1);
		goto L110;

L102:
		if(ip!=2)goto L104;
		if(na!=0)goto L103;

		pp__dradf2(ido,l1,c,ch,wa+iw-1);
		goto L110;

L103:
		pp__dradf2(ido,l1,ch,c,wa+iw-1);
		goto L110;

L104:
		if(ido==1)na=1-na;
		if(na!=0)goto L109;

		pp__dradfg(ido,ip,l1,idl1,c,c,c,ch,ch,wa+iw-1);
		na=1;
		goto L110;

L109:
		pp__dradfg(ido,ip,l1,idl1,ch,ch,ch,c,c,wa+iw-1);
		na=0;

L110:
		l2=l1;
	}

	if(na==1)return;

	for(i=0;i<n;i++)c[i]=ch[i];
}

//////////////////////////////////////////////////////////////////////////////
// END OF FFTPACK CODE
//////////////////////////////////////////////////////////////////////////////

void pp__fft_init(struct PP__FFT* fft, int n)
{
	memset(fft, 0, sizeof *fft);
	fft->n = n;
	assert((fft->wsave = calloc(3*n, sizeof *fft->wsave)) != NULL);
	assert((fft->ifac = calloc(32, sizeof *fft->ifac)) != NULL);
	if (fft->n > 1) pp__drfti1(fft->n, fft->wsave+fft->n, fft->ifac);
}

float* pp__fft_alloc(struct PP__FFT* fft)
{
	float* r = calloc(fft->n, sizeof *r);
	assert(r != NULL);
	return r;
}

void pp__fft_forward(struct PP__FFT* fft, float* data)
{
	if (fft->n==1) return;
	pp__drftf1(fft->n, data, fft->wsave, fft->wsave + fft->n, fft->ifac);
}

static inline void pp__bucket_reset(struct PP* pp)
{
	pp->bucket_cursor = pp->window + PP_WINDOW_SIZE - PP_BUCKET_SIZE;
	pp->bucket_remaining = PP_BUCKET_SIZE;
}

static float* pp__window_alloc()
{
	float* p = calloc(PP_WINDOW_SIZE, sizeof(float));
	assert(p != NULL);
	return p;
}

static inline float pp__hann(float x)
{
	float s = sinf(x * 3.1415);
	return s * s;
}

static inline float pp__fft_window_function(float x)
{
	return pp__hann(x);
}

PP_DEF void pp_init(struct PP* pp, int input_sample_rate_hz)
{
	memset(pp, 0, sizeof *pp);

	pp->pitch0_f0 = PP_PITCH0_F0;

	const int desired_ds_freq = PP_IDEAL_SAMPLE_RATE_HZ;
	int best_decimation_factor = -1;
	int best_diff = 0;
	for (int k = 1; k < 20; k++) {
		int ds_freq = input_sample_rate_hz / k;
		int diff  = ds_freq - desired_ds_freq;
		if (diff < 0) diff = -diff;
		if (k == 1 || diff < best_diff) {
			best_decimation_factor = k;
			best_diff = diff;
		}
	}
	assert(best_decimation_factor >= 1);
	pp->decimation_factor = best_decimation_factor;
	pp->decimated_sample_rate_hz = input_sample_rate_hz / best_decimation_factor;

	pp__fft_init(&pp->fft, PP_WINDOW_SIZE*2); // XXX is this correct?
	pp->fft_data = pp__fft_alloc(&pp->fft);
	assert((pp->fft_window_function = calloc(PP_WINDOW_SIZE, sizeof *pp->fft_window_function)) != NULL);
	for (int i = 0; i < PP_WINDOW_SIZE; i++) {
		pp->fft_window_function[i] = pp__fft_window_function((float)i / (float)PP_WINDOW_SIZE);
		//printf("%f\n", pp->fft_window_function[i]);
	}

	pp->window = pp__window_alloc();
	pp->nsdf = pp__window_alloc();
	pp->spectrogram = pp__window_alloc();

	pp->frame_counter= 0;

	pp__bucket_reset(pp);
}

#if 0
static inline float pp__periodic_parabola(float x, float max, float min)
{
	float x0 = fmodf(x, 1.0);
	float h4 = (max-min)*4;
	return max + h4*x0*x0 - h4*x0;
}
#endif

#if 0
/* Produces various shapes with hih peaks at x=[0, 1, 2,...] and low peaks at
 * x=[0.5, 1.5, 2.5,...]. order=1 interpolates linearly between peaks, while
 * order=2+ flattens between peaks. min/max determines the height of low/high
 * peaks respectively (but be aware that this is for order=1; power is applied
 * after, so correct min/max if necessary; pass e.g. min^(1/3) for order=2,
 * min^(1/5) for order=3, etc) */
static inline float pp__periodic_power_triangle(float x, float min, float max, int order)
{
	const float xp = fmodf(x, 1.0f);
	const float unit_triangle = fabsf(1.0f - xp*2.0f);
	const float minmax_triangle = unit_triangle * (max-min) + min;
	float result = minmax_triangle;
	// return result^1 for order 1, result^3 for order 2, result^5 for order 3, etc
	for (int i = 1; i < order; i++) {
		result *= minmax_triangle*minmax_triangle;
	}
	return result;
}

static inline float pp__spectrogram_scoring_function(float x)
{
	const float min = -1.0f;
	const float max = 1.7f;
	const int order = 1;
	const float offset = 0.0f;
	return pp__periodic_power_triangle(x, min, max, order) + offset;
}
#else
static inline float pp__spectrogram_scoring_function(float x)
{
	if (x < 0.5f) return 0.0f;
	const float min = -1.0f;
	const float max = 1.7f;
	const float xp = fmodf(x, 1.0f);
	const float unit_triangle = fabsf(1.0f - xp*2.0f);
	return unit_triangle * (max-min) + min;
}
#endif

static inline float pp__n_center(int n)
{
	return (float)n * 0.5f - 0.5f;
}

static void pp__parabolic_regression(int n, float* ys, float* b_out, float* c_out)
{
	// consider the parabola:
	//
	//   y = a + b*x + c*x^2
	//
	// we want to find [a b c] so that the parabola best fits our x,y input
	// values. x values are given implicitly; they're sequential and
	// centered around x=0, so x[n-1]+1=x[n], and the sum of x values is 0.
	// to find [a b c] we need to solve the normal equations:
	//
	//   sum[y]     = a*n        + b*sum[x]   + c*sum[x^2]
	//   sum[x*y]   = a*sum[x]   + b*sum[x^2] + c*sum[x^3]
	//   sum[x^2*y] = a*sum[x^2] + b*sum[x^3] + c*sum[x^4]
	//
	// where sum[expression] is the sum of the expression for all input
	// values.


	// start by calculating sums ...
	//float sum_x = 0.0f;
	float sum_y = 0.0f;    // sum[y]
	float sum_xx = 0.0f;   // sum[x^2]
	float sum_xxx = 0.0f;  // sum[x^3]
	float sum_xxxx = 0.0f; // sum[x^4]
	float sum_xy = 0.0f;   // sum[x*y]
	float sum_xxy = 0.0f;  // sum[x^2*y]
	float center = pp__n_center(n);
	for (int i = 0; i < n; i++) {
		float x   = (float)i - center;
		float y   = ys[i];
		float xx  = x*x;
		//sum_x    += x;
		sum_y    += y;
		sum_xx   += xx;
		sum_xxx  += xx*x;
		sum_xxxx += xx*xx;
		sum_xy   += x*y;
		sum_xxy  += xx*y;
	}


	// now arrange the normal equations into matrix form:
	//
	//  B = A*X
	//
	// where B is:
	//    sum[y]
	//    sum[x*y]
	//    sum[x^2*y]
	//
	// and A is:
	//    n        sum[x]     sum[x^2]
	//    sum[x]   sum[x^2]   sum[x^3]
	//    sum[x^2] sum[x^3]   sum[x^4]
	//
	// and we need to find X:
	//    a
	//    b
	//    c
	//
	// the solution is X = inverse(A) * B


	// matrix A
	//    a11    a12    a13
	//    a21    a22    a23
	//    a31    a32    a33
	float a11 = n;
	//float a12 = sum_x; // sum_x=0
	float a13 = sum_xx;
	//float a21 = sum_x; // sum_x=0
	float a22 = sum_xx;
	float a23 = sum_xxx;
	float a31 = sum_xx;
	float a32 = sum_xxx;
	float a33 = sum_xxxx;


	// calculate the inverse of matrix of A
	//
	//   inverse(A) = 1/determinant(A) * cofactor_matrix(A)
	//
	float determinant_a =
		     a11 * (  a22*a33   - a23*a32)
		// - a12 * (  a21*a33   - a23*a31)
		   + a13 * (/*a21*a32*/ - a22*a31);
	float one_over_determinant_a = 1.0f / determinant_a;
	//float ia11 = (   a22*a33     - a23*a32   ) * one_over_determinant_a;
	//float ia12 = ( /*a21*a33*/   - a23*a31   ) * one_over_determinant_a;
	//float ia13 = ( /*a21*a32*/   - a22*a31   ) * one_over_determinant_a;
	float ia21 = ( /*a12*a33*/   - a13*a32   ) * one_over_determinant_a;
	float ia22 = (   a11*a33     - a13*a31   ) * one_over_determinant_a;
	float ia23 = (   a11*a32   /*- a12*a31*/ ) * one_over_determinant_a;
	float ia31 = ( /*a12*a23*/   - a13*a22   ) * one_over_determinant_a;
	float ia32 = (   a11*a23   /*- a13*a21*/ ) * one_over_determinant_a;
	float ia33 = (   a11*a22   /*- a12*a21*/ ) * one_over_determinant_a;


	// matrix B
	float b1 = sum_y;
	float b2 = sum_xy;
	float b3 = sum_xxy;


	// now calculate X = inverse(A) * B
	//float a = ia11*b1 + ia12*b2 + ia13*b3;
	float b = ia21*b1 + ia22*b2 + ia23*b3;
	float c = ia31*b1 + ia32*b2 + ia33*b3;


	// finally return requested values
	//if (a_out) *a_out = a;
	if (b_out) *b_out = b;
	if (c_out) *c_out = c;
}

static inline float pp__parabola_pivot(int n, float* ys)
{
	float b, c;
	pp__parabolic_regression(n, ys, &b, &c);
	return -b / (2*c) + pp__n_center(n);
}

static inline float pp__spectrogram_function(float x)
{
	return log2f(x * 4.0 + 1.0);
}

PP_DEF void pp__flush(struct PP* pp)
{
	pp->frame_counter++;

	float* window = pp->window;
	float* nsdf = pp->nsdf;
	float* spectrogram = pp->spectrogram;

	/* calculate RMS */
	{
		float sum = 0.0f;
		for (int i = 0; i < PP_WINDOW_SIZE; i++) {
			float x = window[i];
			sum += x*x;
		}
		pp->rms = sqrtf(sum / (float)PP_WINDOW_SIZE);
	}

	/* calculate NSDF the slooow way; TODO vectorization? (NOTE:
	 * the rsum part used to be calculated using dual FFTs, because
	 * it's the autocorrelation part and is faster to calculate in
	 * the frequency domain, but since I had to calculate the
	 * normalizing value (msum) as well, why not do it all here and
	 * skip an FFT?) */
	for (int i = 0; i < PP_WINDOW_SIZE; i++) {
		float rsum = 0;
		float msum = 0;
		for (int j = 0; j < (PP_WINDOW_SIZE - i); j++) {
			float a = window[j];
			float b = window[j+i];
			rsum += a*b;
			msum += a*a + b*b;
		}
		nsdf[i] = (2*rsum) / msum;
	}

	/* Use the MPM peak picking algorithm from "A Smarter Way to Find
	 * Pitch" by McLeod/Wyvill; this is fairly robust, but sometimes
	 * produces "octave errors"; these are f0 guesses that are a integer
	 * multiple (larger than 1) of the correct f0 */
	float nsdf_best_score = 0.0;
	float nsdf_best_frequency = 0.0;
	{
		int i = 0;

		// advance to first positive-slope crossing
		for (; i < PP_WINDOW_SIZE && nsdf[i] > 0; i++) {}
		for (; i < PP_WINDOW_SIZE && nsdf[i] <= 0; i++) {}

		const int i0 = i;

		/* find overall max and use it to calculate threshold
		 * */
		float x_max = 0;
		for (i = i0; i < PP_WINDOW_SIZE; i++) {
			float x = nsdf[i];
			if (x > x_max) x_max = x;
		}
		const float k = 0.80; // magical constant woooOOoo
		float threshold = x_max * k;

		/* first highest point within threshold is our f0 guess
		 * */
		for (i = i0; i < PP_WINDOW_SIZE && nsdf[i] < threshold; i++) {}
		x_max = 0;
		int i_max = 0;
		int i_end = (PP_WINDOW_SIZE*5)/6; // skip garbage at the end?
		for (; i < i_end && nsdf[i] >= threshold; i++) {
			float x = nsdf[i];
			if (x > x_max) {
				x_max = x;
				i_max = i;
			}
		}

		const int parabola_window_size = 5;
		if (x_max > 0 && i_max >= (parabola_window_size/2) && i_max < (PP_WINDOW_SIZE-(parabola_window_size/2))) {
			float better_i = pp__parabola_pivot(5, &nsdf[i_max-(parabola_window_size/2)]) + (float)i_max - (float)(parabola_window_size/2);
			nsdf_best_score = x_max;
			nsdf_best_frequency = (float)pp->decimated_sample_rate_hz / better_i;
		}
	}
	pp->nsdf_best_frequency = nsdf_best_frequency;

	// TODO consider stopping here if nsdf_best_score (or
	// nsdf_best_frequency) is bad? I just don't really know what
	// the threshold should be... sometimes the best peaks are
	// actually pretty low... so the best guess for frequency might
	// still be valid...

	/* do forward FFT and calculate spectrogram */
	float* fft_data = pp->fft_data;
	memcpy(fft_data, window, PP_WINDOW_SIZE * sizeof *fft_data);
	float* fft_window_function = pp->fft_window_function;
	for (int i = 0; i < PP_WINDOW_SIZE; i++) fft_data[i] *= fft_window_function[i];
	memset(fft_data + PP_WINDOW_SIZE, 0, PP_WINDOW_SIZE * sizeof *fft_data);
	pp__fft_forward(&pp->fft, fft_data);

	int fi = 0;
	for (int i = 0; i < PP_WINDOW_SIZE; i++) {
		float re = fft_data[fi++];
		float im = fft_data[fi++];
		float mag = re*re + im*im;
		spectrogram[i] = mag;
	}

	/* use spectrogram to possibly correct f0 */
	const float fft_bin_size = ((float)pp->decimated_sample_rate_hz / (float)PP_WINDOW_SIZE) * 0.5;
	const float f0_min = fft_bin_size * 3;
	float best_ratio = 1;
	float best_ratio_score = 0.0f;
	//int dbg = pp->frame_counter == 315;
	//if (dbg) printf("---- nsdf best guess: %f\n", nsdf_best_frequency);
	{
		float try_f0_ratios[] = {
			// NSDF peak finder sometimes triggers on
			// sub-harmonics, so try a few multiples of
			// nsdf_best_frequency
			4, 3, 2,

			// try nsdf_best_frequency as-is; NSDF peak finder
			// isn't _always_ wrong :-)
			1,

			// and NSDF peak finder sometimes triggers on
			// harmonics, so try a bunch of them, assuming 2nd
			// harmonic, 3rd harmonic, and so on
			1.f/2.f, 1.f/3.f, 1.f/4.f, 1.f/5.f, 1.f/6.f, 1.f/7.f, 1.f/8.f,
		};
		const int n_try_f0_ratios = sizeof(try_f0_ratios) / sizeof(try_f0_ratios[0]);
		for (int i = 0; i < n_try_f0_ratios; i++) {
			float ratio = try_f0_ratios[i];
			float f0_guess = nsdf_best_frequency * ratio;
			if (f0_guess < f0_min) continue;
			float score = 0.0;
			float x = 0.0;
			const float xstep = fft_bin_size / f0_guess;
			//if (dbg) printf("trying %f ...\n", f0_guess);
			for (int j = 0; j < PP_WINDOW_SIZE; j++) {
				score += pp__spectrogram_function(spectrogram[j]) * pp__spectrogram_scoring_function(x);
				//if (dbg) printf("   %d (x=%.4f) score += %f * %f\n", j, x, spectrogram[j], pp__spectrogram_scoring_function(x));
				x += xstep;
			}
			//if (dbg) printf("guess %f - score %f\n\n", f0_guess, score);
			if (i == 0 || score > best_ratio_score) {
				best_ratio_score = score;
				best_ratio = ratio;
			}
		}
	}

	float spectrogram_corrected_best_frequency = nsdf_best_frequency * best_ratio;
	//if (dbg) printf("...aaaand the winner is: %f\n\n", spectrogram_corrected_best_frequency);
	pp->spectrogram_corrected_best_frequency = spectrogram_corrected_best_frequency;

	float f0 = spectrogram_corrected_best_frequency; // TODO peak surfing!
	float confidence = nsdf_best_score; // XXX?
	if (confidence < 0.5 || pp->rms < pp->noise_gate_rms_threshold || f0 < f0_min) {
		f0 = 0;
	}
	pp->f0 = f0;

	if (pp->pitch0_f0 > 0) {
		if (f0 > 0) {
			pp->pitch = (log2f(f0) - log2f(pp->pitch0_f0)) * 12.0f;
		} else {
			pp->pitch = 0;
		}
	}

	pp__bucket_reset(pp);
	pp->window_scroll_is_pending = 1;
}

PP_DEF void pp__scroll_window(struct PP* pp)
{
	memmove(
		&pp->window[0],
		&pp->window[PP_BUCKET_SIZE],
		(PP_WINDOW_SIZE - PP_BUCKET_SIZE) * sizeof(*pp->window)
	);
	pp->window_scroll_is_pending = 0;
}

# endif

// This is free and unencumbered software released into the public domain.
//
// Anyone is free to copy, modify, publish, use, compile, sell, or
// distribute this software, either in source code form or as a compiled
// binary, for any purpose, commercial or non-commercial, and by any
// means.
//
// In jurisdictions that recognize copyright laws, the author or authors
// of this software dedicate any and all copyright interest in the
// software to the public domain. We make this dedication for the benefit
// of the public at large and to the detriment of our heirs and
// successors. We intend this dedication to be an overt act of
// relinquishment in perpetuity of all present and future rights to this
// software under copyright law.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
// EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
// MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
// IN NO EVENT SHALL THE AUTHORS BE LIABLE FOR ANY CLAIM, DAMAGES OR
// OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE,
// ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
// OTHER DEALINGS IN THE SOFTWARE.
//
// For more information, please refer to <http://unlicense.org/>
