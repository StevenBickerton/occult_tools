#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_fft_real.h>
#include <gsl/gsl_fft_halfcomplex.h>
#include <fftw3.h>

#include "detection.h"
#include "statistics.h"
#include "fft-fftw.h"

void FFTWtimeseries(TS *pts){
  
  assert(pts != NULL);
  pts->fI = det_malloc(pts->n2 * sizeof(FLOAT));

  const int N=pts->n;
  for(int i=0;i<N;i++) {
    pts->fI[i] = pts->tI[i].I - 1.0;
  }

  const int N2 = pts->n2;
  for(int i=N;i<N2;i++) {
    pts->fI[i] = 0.0;
  }

  pts->fftwplan = fftwFTYPE_plan_r2r_1d((int) N2, pts->fI, pts->fI, FFTW_R2HC, FFTW_ESTIMATE);
  fftwFTYPE_execute( pts->fftwplan );
  fftwFTYPE_destroy_plan( pts->fftwplan );

  pts->FFTisdone = 1;
}



int fftwconvolve(FLOAT *corr, TS *pts, SHADOW *pfresPatt, int i_vel, int *Neff) {

  int i;
  const int N=pts->n;
  const int n=pfresPatt->nI[i_vel];
  FLOAT real, imag;
  
  /* FFT the time series */
  if (! pts->FFTisdone) {
    FFTWtimeseries(pts);
  }
  
  const int N2 = pts->n2;
  
  /* FFT the Kernel */
  for(i=0;i<n;i++) {
    corr[i] = pfresPatt->I[i_vel][i] - 1.0;
  }
  for(i=n;i<N2;i++) {
    corr[i] = 0.0;
  }
  
  pts->fftwplan = fftwFTYPE_plan_r2r_1d((int)N2, corr, corr, FFTW_R2HC, 
				       FFTW_ESTIMATE);
  fftwFTYPE_execute(pts->fftwplan);
  fftwFTYPE_destroy_plan(pts->fftwplan);
  
  /* do the product */
  FLOAT gauss;
  /* event width is 3.4 Fsu ... filter is 8x that */
  FLOAT gaussFreqCoeff = (SMOOTH_WIDTH_FSU*
			   1.0*pfresPatt->fsu/pfresPatt->vRet[i_vel]) /
    ( ((FLOAT) N2/N) * pts->tI[N-1].t - pts->tI[0].t );
  gaussFreqCoeff = 0.5 * M_PI * M_PI * gaussFreqCoeff * gaussFreqCoeff;
  FLOAT jclip = 5.0*sqrt(1.0/(2.0*gaussFreqCoeff));

  FLOAT jd;
  for(i=1;i<N2/2;i++) {
    jd = (FLOAT) i;
    gauss = (i<jclip) ? exp(-jd*jd * gaussFreqCoeff) : 0.0;
    real = pts->fI[i] * corr[i] + pts->fI[N2-i] * corr[N2-i];
    imag = pts->fI[N2-i] * corr[i] - pts->fI[i] * corr[N2-i];
    corr[i] = real * (1.0 - gauss);
    corr[N2-i] = imag * (1.0 - gauss);
  }
  corr[0] = 0;

  /* inverse transform back */
  pts->fftwplan = fftwFTYPE_plan_r2r_1d((int)N2, corr, corr, FFTW_HC2R,
				       FFTW_ESTIMATE);
  fftwFTYPE_execute(pts->fftwplan);
  fftwFTYPE_destroy_plan(pts->fftwplan);

  const FLOAT inv_N2 = 1.0 / N2;
  for (i=0;i<N;i++) {
    corr[i] *= inv_N2;
  }

  *Neff = N;
  return(N);
}

void fftw_xcorr_norm(FLOAT *corr, FLOAT *corrNorm, TS *ptsa, SHADOW *pfresPatt, FLOAT rms) {

  int i;
  FLOAT clip = 16.0*rms*rms;
  int const N = ptsa->n;
  int const N2 = ptsa->n2;
  FLOAT gauss;
  FLOAT gaussFreqCoeff = NORM_WIDTH_SEC/
    ( ((FLOAT) N2/N) * ptsa->tI[N-1].t - ptsa->tI[0].t );
  gaussFreqCoeff = 0.5 * M_PI * M_PI * gaussFreqCoeff * gaussFreqCoeff;
  FLOAT jclip = 5.0*sqrt(1.0/(2.0*gaussFreqCoeff));
  
  FLOAT jd;


  /* square the time series */
  for (i=0;i<N;i++) {
      corrNorm[i] = corr[i] * corr[i];
      if (corrNorm[i] > clip)
	corrNorm[i] = rms*rms; 	/* put in a dummy for outliers */
  }
  for (i=N;i<N2;i++)
    corrNorm[i] = rms*rms; 	/* dummy padding */
  
  /* fft the squared-ts to the Freq dom. */
  ptsa->fftwplan = fftwFTYPE_plan_r2r_1d(N2, corrNorm, corrNorm, FFTW_R2HC,
				    FFTW_ESTIMATE);
  fftwFTYPE_execute(ptsa->fftwplan);
  fftwFTYPE_destroy_plan(ptsa->fftwplan);

  /* product with fft-gaussian in Freq dom. */
  for (i=1;i<N2/2;i++) {
    jd = (FLOAT) i; // must cast to FLOAT or j*j=inf for N>92682
    gauss = (i<jclip) ? exp(-jd*jd * gaussFreqCoeff) : 0.0;
    corrNorm[i]            *= gauss; // real
    corrNorm[N2-i]       *= gauss; // cplx
  }


  /* invert back to time domain */
  ptsa->fftwplan = fftwFTYPE_plan_r2r_1d(N2, corrNorm, corrNorm, FFTW_HC2R,
				    FFTW_ESTIMATE);
  fftwFTYPE_execute(ptsa->fftwplan);
  fftwFTYPE_destroy_plan(ptsa->fftwplan);

  /* square root is local rms */
  /* divide corr by local rms to get corrNorm */
  const FLOAT inv_N2 = 1.0/( (FLOAT) N2 );
  for (i=0;i<N;i++) {
    corrNorm[i] = corr[i] / sqrt(inv_N2*corrNorm[i]);
  }

}
