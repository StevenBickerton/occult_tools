#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <gsl/gsl_errno.h>

#include "detection.h"
#include "statistics.h"
#include "fft-gsl.h"

void FFTtimeseries(TS *pts){
  
  assert(pts != NULL);
  pts->fI = malloc(pts->n2 * sizeof(FLOAT));

  const int N = pts->n;
  for(int i=0;i<N;i++) {
    pts->fI[i] = pts->tI[i].I - 1.0;
  }

  const int N2 = pts->n2;
  for(int i=N;i<N2;i++) {
    pts->fI[i] = 0.0;
  }

  pts->Rwavetab  = gsl_fft_real_wavetable_FTYPE_alloc(N2);
  pts->HCwavetab = gsl_fft_halfcomplex_wavetable_FTYPE_alloc(N2);
  pts->workspace = gsl_fft_real_workspace_FTYPE_alloc(N2);
  gsl_fft_real_FTYPE_transform(pts->fI, 1, N2, pts->Rwavetab, pts->workspace);
  pts->FFTisdone = 1;
}



int fftconvolve(FLOAT *corr, TS *pts, SHADOW *pfresPatt, int i_vel, int *Neff) {

  int i;
  int const N=pts->n;
  int const n=pfresPatt->nI[i_vel];
  int const N2 = pts->n2;
  FLOAT real, imag;
  
  /* FFT the time series */
  if (! pts->FFTisdone)
    FFTtimeseries(pts);

  
  /* FFT the Kernel */
  for(i=0;i<n;i++)
    corr[i] = pfresPatt->I[i_vel][i] - 1.0;
  for(i=n;i<N2;i++)
    corr[i] = 0.0;

  gsl_fft_real_FTYPE_transform(corr, 1, N2, pts->Rwavetab, pts->workspace);

  FLOAT gauss;
  /* event width is 3.4 Fsu ... filter is SMTH_WIDTHx that */
  FLOAT gaussFreqCoeff = (SMOOTH_WIDTH_FSU * 
			   3.4*pfresPatt->fsu/pfresPatt->vRet[i_vel]) /
    ( ((FLOAT) N2/N) * pts->tI[N-1].t - pts->tI[0].t );
  gaussFreqCoeff = 0.5 * M_PI * M_PI * gaussFreqCoeff * gaussFreqCoeff;
  FLOAT jclip = 5.0*sqrt(1.0/(2.0*gaussFreqCoeff));
  
  FLOAT jd;

  /* do the product */
  corr[0] = pts->fI[0] * corr[0];
  for(i=1;i<N2/2;i++) {
    jd = (FLOAT) i;
    gauss = (i<jclip) ? exp(-jd*jd * gaussFreqCoeff) : 0.0;
    real = pts->fI[2*i-1] * corr[2*i-1] + pts->fI[2*i] * corr[2*i];
    imag = pts->fI[2*i] * corr[2*i-1] - pts->fI[2*i-1] * corr[2*i];
    corr[2*i-1] = real * (1.0 - gauss);
    corr[2*i] = imag * (1.0 - gauss);
  }
  
  /* inverse transform back */
  gsl_fft_halfcomplex_FTYPE_inverse(corr, 1, N2, pts->HCwavetab, pts->workspace);

  *Neff = N;
  return(N);
}

void fft_xcorr_norm(FLOAT *corr, FLOAT *corrNorm, TS *ptsa, SHADOW *pfresPatt, FLOAT rms) {

  int i;
  FLOAT clip = 16.0*rms*rms;
  int const N = ptsa->n;
  int const N2 = ptsa->n2;
  FLOAT gauss;
  FLOAT gaussFreqCoeff = NORM_WIDTH_SEC/(ptsa->tI[N-1].t - ptsa->tI[0].t);
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
  gsl_fft_real_FTYPE_transform(corrNorm, 1, N2, ptsa->Rwavetab, ptsa->workspace);
  
  /* product with fft-gaussian in Freq dom. */
  for (i=1;i<N2/2;i++) {
    jd = (FLOAT) i; // must cast to FLOAT or j*j=inf for N>92682
    gauss = (i<jclip) ? exp(-jd*jd * gaussFreqCoeff) : 0.0;
    corrNorm[2*i-1]           *= gauss; // real
    //corrNorm[2*(N-i-1)]     *= gauss; // real -ve
    corrNorm[2*i]         *= gauss; // cplx
    //corrNorm[2*(N-i-1)+1]   *= gauss; // cplx -ve
  }


  /* invert back to time domain */
  gsl_fft_halfcomplex_FTYPE_inverse(corrNorm, 1, N2, ptsa->HCwavetab, ptsa->workspace);

  /* square root is local rms */
  /* divide corr by local rms to get corrNorm */
  for (i=0;i<N;i++) {
    corrNorm[i] = corr[i] / sqrt(corrNorm[i]);
  }

}
