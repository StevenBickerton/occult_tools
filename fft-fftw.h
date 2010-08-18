#ifndef FFTW_LIBRARY
#define FFTW_LIBRARY 1 

#if defined(USE_FLOAT)
#define fftwFTYPE_plan_r2r_1d fftwf_plan_r2r_1d
#define fftwFTYPE_execute fftwf_execute
#define fftwFTYPE_destroy_plan fftwf_destroy_plan
#else
#define fftwFTYPE_plan_r2r_1d fftw_plan_r2r_1d
#define fftwFTYPE_execute fftw_execute
#define fftwFTYPE_destroy_plan fftw_destroy_plan
#endif

void FFTWtimeseries(TS *pts);
int fftwconvolve(FLOAT *corr, TS *pts, SHADOW *pfresPatt, int i_vel, int *Neff);
void fftw_xcorr_norm(FLOAT *corr, FLOAT *corrNorm, TS *ptsa, SHADOW *pfresPatt, FLOAT rms);

#endif
