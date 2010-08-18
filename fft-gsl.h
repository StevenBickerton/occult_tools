#ifndef FFTGSL_LIBRARY
#define FFTGSL_LIBRARY 1 

#if defined(USE_FLOAT)
#define gsl_fft_real_wavetable_FTYPE_alloc gsl_fft_real_wavetable_float_alloc
#define gsl_fft_halfcomplex_wavetable_FTYPE_alloc gsl_fft_halfcomplex_wavetable_float_alloc
#define gsl_fft_real_workspace_FTYPE_alloc gsl_fft_real_workspace_float_alloc
#define gsl_fft_real_FTYPE_transform gsl_fft_real_float_transform
#define gsl_fft_halfcomplex_FTYPE_inverse gsl_fft_halfcomplex_float_inverse
#else
#define gsl_fft_real_wavetable_FTYPE_alloc gsl_fft_real_wavetable_alloc
#define gsl_fft_halfcomplex_wavetable_FTYPE_alloc gsl_fft_halfcomplex_wavetable_alloc
#define gsl_fft_real_workspace_FTYPE_alloc gsl_fft_real_workspace_alloc
#define gsl_fft_real_FTYPE_transform gsl_fft_real_transform
#define gsl_fft_halfcomplex_FTYPE_inverse gsl_fft_halfcomplex_inverse
#endif

void FFTtimeseries(TS *pts);
int fftconvolve(FLOAT *corr, TS *pts, SHADOW *pfresPatt, int i_vel, int *Neff);
void fft_xcorr_norm(FLOAT *corr, FLOAT *corrNorm, TS *ptsa, SHADOW *pfresPatt, FLOAT rms);

#endif
