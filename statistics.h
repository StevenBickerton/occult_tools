#ifndef LIBSTATS
#define LIBSTATS 1

/*            */
/* Made with makeScript, Thu Jul 21, 2005  12:13:58 DST */
/* Host: kuiper */
/* Working Directory: /1/home/bickersj/sandbox/cdiffracSim  */

// constants
#define PI 3.141592654

#define STAT_MAX_LINE_LENGTH 80
#define CLIP_TOLERANCE 0.001
#define CLIP_ITERATIONS 3
#define USE_ALL        0

// macros
#ifndef DB
#define DB(...) { fprintf(stderr, __VA_ARGS__); }
//#else
//#define DB(...) { }
#endif

typedef struct {
  FLOAT mean;
  FLOAT rms;
} STAT;

// max's and mins
FLOAT max_d(FLOAT *values, int n, int *index);
FLOAT min_d(FLOAT *values, int n, int *index);
int below_d (FLOAT *values, int n, FLOAT thresh, FLOAT *below, int *index, int nMax);
int above_d (FLOAT *values, int n, FLOAT thresh, FLOAT *above, int *index, int nMax);

int window (FLOAT *values, int n, FLOAT lo, FLOAT hi, FLOAT *in_range, int *index, int nMax);

// mean and rms  ... median, mode eventually
//FLOAT mean_d(FLOAT *values, int n, int subsamp);
//FLOAT rms_d(FLOAT *values, int n, int subsamp);
FLOAT sum_d(FLOAT *values, int n);
FLOAT moment(FLOAT *values, int n, int r);
FLOAT moment_mean(FLOAT *values, int n, int r, FLOAT mean);
FLOAT skew_d(FLOAT *values, int n, FLOAT mean, FLOAT rms);
FLOAT kurt_d (FLOAT *values, int n, FLOAT mean, FLOAT rms);
STAT meanRms_d (FLOAT *values, int n, int subsamp);

// routines with sigma clipping
//FLOAT rmsclip_d(FLOAT *values, int n, int subsamp, FLOAT nSigma, FLOAT sigma, int iter);
//FLOAT meanclip_d(FLOAT *values, int n, int subsamp, FLOAT nSigma, FLOAT sigma, int iter);
STAT meanRmsClip_d(FLOAT *values, int n, FLOAT nSigma, int stride);

FLOAT factorial (int N);
FLOAT gammaN (FLOAT n);
FLOAT chi2dist(FLOAT x2, int v);
FLOAT Predchi2(FLOAT x2, int v);

void statErrQuit(char *msg, ...);

#endif // LIBSTATS
