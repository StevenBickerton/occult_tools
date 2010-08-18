/*            */
/* Made with makeScript, Thu Jul 21, 2005  12:13:58 DST */
/* Host: kuiper */
/* Working Directory: /1/home/bickersj/sandbox/cdiffracSim  */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <stdarg.h>
#include <assert.h>
// #include <gsl/gsl_integration.h>

#include "detection.h"
#include "statistics.h"

//#define M_SQRTPI 1.77245385090552
#define M_SQRT2PI 2.506628274631


/////////////////////////////////////////////////////////////////
FLOAT min_d(FLOAT *values, int n, int *index) {

  int i;
  FLOAT min=values[0];

  for (i=0; i<n; i++ ) {
    min = (values[i]<min) ? (values[i]) : (min);
  }

  *index = i;
  return min;
}


/////////////////////////////////////////////////////////////////
FLOAT max_d(FLOAT *values, int n, int *index) {

  int i;
  FLOAT max=values[0];

  for (i=0; i<n; i++ ) {
    max = (values[i]>max) ? (values[i]) : (max);
  }

  *index = i;
  return max;
}



/////////////////////////////////////////////////////////////////
int below_d (FLOAT *values, int n, FLOAT thresh, FLOAT *below, int *index, int nMax) {

  int i,found=0;
  
  for (i=0; i<n; i++) {
    if ( values[i] < thresh ) {
      below[found] = values[i];
      index[found++] = i;
      if (found >= nMax) {
	statErrQuit("hits exceeded array size");
      }
    }
  }

  return found;
}




/////////////////////////////////////////////////////////////////
int above_d (FLOAT *values, int n, FLOAT thresh, FLOAT *above, int *index, int nMax) {

  int i,found=0;
  
  for (i=0; i<n; i++) {
    if ( values[i] > thresh ) {
      above[found] = values[i];
      index[found++] = i;
      if (found >= nMax ) {
	statErrQuit("hits exceeded array size");
      }
    }
  }

  return found;
}


/////////////////////////////////////////////////////////////////
FLOAT sum_d(FLOAT *values, int n) {

  int i;
  FLOAT sum=0;
  
  for (i=0; i<n ; i++) {
    if (!isnan(values[i]) ) {
      sum += values[i];
    }
  }
  
  return (sum);
}


/////////////////////////////////////////////////////////////////
int window (FLOAT *values, int n, FLOAT lo, FLOAT hi, FLOAT *in_range, int *index, int nMax) {

  int i,found=0;
  for (i=0; i<n; i++) {
    if ( values[i] > lo  &&  values[i] < hi ) {
      in_range[found] = values[i];
      index[found] = i;
      found++;
      if (found >= nMax) {
	statErrQuit("hits exceeded array size");
      }
    }
  }
  return found;
}



/////////////////////////////////////////////////////////////////
/* FLOAT mean_d(FLOAT *values, int n, int subsamp) { */

/*   int i, step; */
/*   FLOAT sum=0; */


/*   subsamp = (subsamp<n) ? (subsamp) : (n);  // only have n values */
/*   subsamp = (subsamp==0) ? (n) : (subsamp);  // if 0 use all values */

/*   if (n==0) return 0; */

/*   step = (((int) n/subsamp)>0) ? ((int) n/subsamp) : 1; */

/*   for (i=0; i<subsamp ; i++) */
/*     sum += values[i*step]; */
  
/*   return (sum / subsamp); */
/* } */




/* //////////////////////////////////////////////////////////// */
/* FLOAT rms_d(FLOAT *values, int n, int subsamp) { */

/*   int i, step; */
/*   FLOAT sum, sumsqr; */
/*   sum = 0; */
/*   sumsqr = 0; */

/*   subsamp = (subsamp<n) ? (subsamp) : (n);  // only have n values */
/*   subsamp = (subsamp==0) ? (n) : (subsamp);  // if 0 use all values */

/*   if (n==0) return 0; */

/*   step = (int) n/subsamp; */

/*   for (i=0; i<subsamp; i++) { */
/*     sum += values[i*step]; */
/*     sumsqr += (values[i*step] * values[i*step]); */
/*   } */
   
/*   return  sqrt(   sumsqr/subsamp    -   sum*sum/(subsamp*subsamp)   )  ; */

/* } */


FLOAT moment(FLOAT *values, int n, int r) {

  int i,j;
  FLOAT *array, moment;

  if (n==0) {
    return 0;
  }
  array = det_malloc(n * sizeof(FLOAT));

  j=0;
  for (i=0; i<n; i++)  {
    if (!isnan(values[i])) {
      array[j++] = pow(values[i], r);
    }
  }
  moment = (j) ? (sum_d(array, j)/j) : 0;

  free(array);
  return moment;
}

FLOAT moment_mean(FLOAT *values, int n, int r, FLOAT mean) {

  int i, j;
  FLOAT *array, moment_mean;

  if (n==0) {
    return 0;
  }
  array = det_malloc(n * sizeof(FLOAT));

  j=0;
  for (i=0; i<n; i++) {
    if (! isnan(values[i]) ) {
      array[j++] = values[i] - mean;
    }
  }
  moment_mean = moment(array, j, r);
  free(array);

  return moment_mean;
}

FLOAT skew_d(FLOAT *values, int n, FLOAT mean, FLOAT rms) {

  FLOAT third_moment = moment_mean(values, n, 3, mean);
  FLOAT skew = (rms) ? ( third_moment / (rms*rms*rms) ) : 0;
  return skew;
}


FLOAT kurt_d (FLOAT *values, int n, FLOAT mean, FLOAT rms) {

  FLOAT fourth_moment = moment_mean(values, n, 4, mean);
  FLOAT kurt = (rms) ? ( fourth_moment / (rms*rms*rms*rms) ) : 0;
  return kurt;
}


STAT meanRms_d (FLOAT *values, int n, int subsamp) {


  // only have n values ... if set to 0, use all
  if (subsamp >= n || subsamp == 0) {
    subsamp = n;
  }

  // return zeros if there are no data.
  STAT stats;
  if (n==0) {
    stats.mean = 0;
    stats.rms = 0;
    return stats;
  }

  int step = n/subsamp;

  // loop through and get the sum and sum of squares
  int count = 0;
  FLOAT sum = 0;
  FLOAT sumsqr = 0;
  for (int i=0; i<subsamp-1; i++) {
    if ( values[i*step] != 1.000000 && values[i*step] != 0.000000 ) {
      sum += values[i*step];
      sumsqr += (values[i*step] * values[i*step]);
      count++;
    }
  }

  stats.mean = sum/count;
  stats.rms  = sqrt(  sumsqr/count  -  stats.mean*stats.mean  )  ;
  
  return (stats);
}



////////////////////////////////////////////////////////////////
STAT meanRmsClip_d(FLOAT *values, int n, FLOAT nSigma, int stride){

  int i, j, k, nGood, N;
  FLOAT mean_last, *goodvalues;
  STAT stats;
  int subsamp;
  subsamp = (n - 1) / stride;
 
  if (n == 0) {
    stats.mean = 0.0;
    stats.rms = -1.0;
    return stats;
  }
  goodvalues= det_malloc(n * sizeof(FLOAT));

  // starting values of mean and rms
  N = (int) (n / stride - stride);
  stats = meanRms_d(values, n, subsamp);
  mean_last = stats.mean;

  // iterate through and take only the values within nSigma of the mean
  for(k=0; k<CLIP_ITERATIONS; k++) {
    
    nGood = 0;
    for (i=0; i< N; i++) {
      j = i * stride;
      if ( fabs(values[j] - stats.mean) < nSigma*stats.rms ) {
	goodvalues[nGood++] = values[j];
      }
    }
    
    mean_last = stats.mean;

    // get the improved values
    subsamp = (nGood - 1) / stride;
    stats = meanRms_d(goodvalues, nGood, subsamp);

    if (fabs((mean_last - stats.mean)/mean_last) < CLIP_TOLERANCE) {
      break;
    }
  }

  free(goodvalues);
  return ( stats );
}

FLOAT factorial (int N) {

  int i;
  FLOAT Nf = (FLOAT) N;

  if ( N == 0 ) {
      return 1;
  }

  FLOAT fact = 1.0;

  // stirling's approximation
  if (N > 50 ) {
    fact = (M_SQRT2PI * sqrt(Nf) * pow(Nf, N) * exp(-Nf)) *
      (1.0 + 1.0 / (12.0*Nf) + 
       1.0   / (288.0*Nf*Nf) + 
       139.0 / (51840*Nf*Nf*Nf) );
    
    // if N is small, just compute it
  } else {
    for(i=N; i>0; i--) {
      fact *= i;
    }
  }

  return fact;
}


FLOAT gammaN (FLOAT n) {
  
  int m;
  int i;
  
  FLOAT gammaN;
  
  // if n is a positive integer just use factorial
  if ( n>0 && fabs(n-trunc(n))<1.0e-12  ) {
    gammaN = factorial( (int) trunc(n+0.5) - 1 );
    
    // if n is a negative integer
  } else if ( n<0 && fabs(n-trunc(n)) < 1.0e-12) {
    fprintf(stderr, "Warning: GammaN(x) not defined for neg. integers. Return 0\n");
    return 0;
    
    // if n is a positive half-integer
  } else if ( n>0 &&  fabs(n-0.5 - trunc(n-0.5) < 1.0e-12)  ) {
    m = (int) trunc(n);        
    gammaN = M_SQRTPI / pow(2.0, m);
    for(i=1; i<=m; i++) {
      gammaN *= 2.0 * i - 1.0;
    }
    
    // if n is a negative half-integer
  } else if (n < 0  &&  fabs(n - 0.5 - trunc(n-0.5) < 1.0e-12) ) {
    
    m = (int) -(n - 0.5);
    gammaN = M_SQRTPI * (m%2 ? (-1.0) : 1.0) * pow(2.0, m);
    for(i=1; i<=m; i++) {
      gammaN /= 2.0 * i - 1.0;
    }
    
    // if n is an arb num (warning ... only good to a few digits)
  } else {
    n -= 1.0;
    gammaN = M_SQRT2PI * sqrt(n) * pow(n, n) * exp(-n) *
      (1.0 + 1.0/(12.0*n) + 
       1.0 / (288.0*n*n) + 
       139.0 / (51840*n*n*n) );
  }
  
  return gammaN;
}


FLOAT chi2dist (FLOAT x2, int v) {

  FLOAT term1 = 1.0/( pow(2.0, (v / 2.0)) * gammaN(v / 2.0) );
  FLOAT term2 = pow(x2, ((v - 2.0) / 2.0) );
  FLOAT term3 = exp(-x2 / 2.0);
  
  return term1 * term2 * term3;
}

FLOAT chi2dist_int (FLOAT x, void *params) {
  return chi2dist(x, *(int *) params);
}

FLOAT Predchi2 (FLOAT x2, int v) {
  
/*   int const N = 2000; */
/*   FLOAT integral = 0; */
/*   FLOAT interr; */
/*   x2 *= v;  // because it's a reduced X2 */

/*   gsl_function F; */
/*   gsl_integration_workspace *wkspc = gsl_integration_workspace_alloc(N); */
/*   F.function = &chi2dist_int; */
/*   F.params = &v; */
/*   gsl_integration_qags (&F, 0, x2, 0, 1e-8, N, wkspc, &integral, &interr); */

  const int N = 2000;
  int i;
  FLOAT integral = 0, x2_prev = 0;
  FLOAT x2i, x2tmp;
  const FLOAT inv_N = 1.0/N;
  x2 *= v;  // because it's a reduced X2
  const FLOAT dx2 = x2*inv_N;
  
  for (i=1;i<N;i++) {
    x2i = i*dx2;
    x2tmp = chi2dist(x2i,v);
    integral += 0.5 * (x2tmp + x2_prev);
    x2_prev = x2tmp;
  }
  integral *= dx2;

  if (integral > 1.0) {
    integral = 1.0;
  }

  return(integral);

}


/////////////////////////////////////////////////////////////////
/* FLOAT meanclip_d(FLOAT *values, int n, int subsamp, FLOAT nSigma, FLOAT sigma, int iter) { */

/*   int i, j, k, nGood, step, hstep; */
/*   FLOAT mean, rms, goodvalues[n], val_tmp; */

/*   subsamp = (subsamp<n) ? (subsamp) : (n);  // only have n values */
/*   step = (int) n/subsamp;    // step size for subsampling */
/*   hstep = step / 4;          // half step */

/*   mean = mean_d(values, n, subsamp); */
/*   rms = (sigma) ? (sigma) : rms_d(values, n, subsamp);  // avoid calc if given */

/*   // plan to go thru once, but then start again to complete the sample */
/*   //      (short some due to noise spikes) */

  
/*   for (k=0; k<iter; k++) { */
/*     nGood = 0; */

/*     for (j=0; j<4; j++) { */
/*       for (i=0; i<subsamp; i++) { */
	
/* 	// if above threshhold and not a mask value (1.000000) */
/* 	val_tmp = values[i*step + j*hstep]; */
/* 	if (  (fabs(val_tmp - mean) < nSigma*rms)  && (val_tmp != 1.000000) ) */
/* 	  goodvalues[nGood++] = val_tmp; */
	
/* 	if (nGood >= subsamp) break; */
/*       } */
      
/*       if (nGood >= subsamp) break; */
/*     }   */
    
/*     mean = mean_d(goodvalues, nGood, subsamp); */
/*     rms = rms_d(goodvalues, nGood, subsamp); */
/*   } */

/*   return mean; */
/* } */



/* FLOAT rmsclip_d(FLOAT *values, int n, int subsamp, FLOAT nSigma, FLOAT sigma, int iter) { */

/*   int i, j, k, nGood, step, hstep, cycles; */
/*   FLOAT mean, rms, goodvalues[n], val_tmp; */

/*   subsamp = (subsamp<n) ? (subsamp) : (n);  // only have n values */
/*   step = (int) n/subsamp; */
/*   cycles = 4; */
/*   hstep = step / cycles; */

/*   mean = mean_d(values, n, subsamp); */
/*   rms = (sigma) ? (sigma) : rms_d(values, n, subsamp);  // avoid calc if given */

/*   // iterate thru to clip away noise */
/*   for (k=0; k<iter; k++) { */

/*     // because some points will be rejected as noise, make four possible passes */
/*     // thru the data to get the requested number of samples to use. */
/*     nGood = 0; */
/*     for (j=0; j<cycles; j++) { */
/*       for (i=0; i<subsamp; i++) { */

/* 	val_tmp = values[i*step + j*hstep]; */
/* 	if (  ( (fabs(val_tmp) - mean) < nSigma*rms) && (val_tmp != 0.000000) ) */
/* 	  goodvalues[nGood++] = val_tmp; */

/* 	if (nGood >= subsamp) break; */
/*       } */

/*       if (nGood >= subsamp) break; */
/*     } */
    
/*     mean = mean_d(goodvalues, nGood, subsamp); */
/*     rms = rms_d(goodvalues, nGood, subsamp); */
/*   } */

/*   return rms; */
/* } */



void statErrQuit (char *msg, ...) {

  char s[STAT_MAX_LINE_LENGTH];

  va_list args;
  va_start(args,msg);
  vsprintf(s, msg, args);
  va_end(args);

  perror(s);
  exit(EXIT_FAILURE);
}

