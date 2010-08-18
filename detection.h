/*            */
/* Made with makeScript, Fri Jun 24, 2005  15:58:46 DST */
/* Host: kuiper */
/* Working Directory: /home/bickersj/sandbox/cdiffracSim  */

//#include <stdlib.h>
//#include <stdio.h>

//#include <dirent.h>       // directory calls
//#include <sys/types.h>    // directory calls

#ifndef DETECT_LIBRARY
#define DETECT_LIBRARY 1

#include<fftw3.h>
#define USE_FFTW 1
// #define USE_FLOAT 1

#if defined(USE_FLOAT)

typedef float FLOAT;
#define fftwFTYPE_plan fftwf_plan
#define fftwFTYPE_complex fftwf_complex

#include<gsl/gsl_fft_real_float.h>
#include<gsl/gsl_fft_halfcomplex_float.h> 
#define  gsl_fft_real_workspace_FTYPE gsl_fft_real_workspace_float
#define  gsl_fft_real_wavetable_FTYPE gsl_fft_real_wavetable_float
#define  gsl_fft_halfcomplex_wavetable_FTYPE gsl_fft_halfcomplex_wavetable_float

#else

typedef double FLOAT;
#define fftwFTYPE_plan fftw_plan
#define fftwFTYPE_complex fftw_complex

#include<gsl/gsl_fft_real.h>
#include<gsl/gsl_fft_halfcomplex.h> 
#define  gsl_fft_real_workspace_FTYPE gsl_fft_real_workspace
#define  gsl_fft_real_wavetable_FTYPE gsl_fft_real_wavetable
#define  gsl_fft_halfcomplex_wavetable_FTYPE gsl_fft_halfcomplex_wavetable

#endif



// define constants
#define PI 3.141592654
#define TWOPI (2.0*PI)
#define PIHALF (PI/2.0)
#define G  6.67259e-11    // 
#define Mo 1.989e30
#define Re 1.496e11
#define DEG 180.0/PI      // radians to degrees
#define RAD PI/180.0      // deg to rad
#define AU_M  1.496e11

#define MAX_PATTERN_POINTS 512
#define MAX_LINE_LENGTH 256*sizeof(char)
#define MAX_FILENAME 256*sizeof(char)
#define MAX_FRESFILES 2048
#define MAX_KERNEL_OVERSAMPLING 8
#define MAX_KNOWN_OBJECTS 256

// allow for a high precision option which keeps small signals
//  values below these limits will be ignored by the code
#ifdef PRECISION
#define OFFSET_MIN_I_LIMIT 0.00001
#define KERNEL_LO_I_CUT 0.00001        // I smaller than this in kernel ignored
#else
#define OFFSET_MIN_I_LIMIT 0.001
#define KERNEL_LO_I_CUT 0.0025
#endif

#define MAX_OFFSETS 100
#define N_PAST_MIN_I_LIMIT 10
#define MAX_KERNEL_POINTS 5000
//#define MAX_TS_POINTS 400000
#define N_VELOCITIES 1               // total number of velocites to test (odd)
#define T_STEP 0.001                 // time resolution to use for integral
#define ITMP_SIZE  16384              // size of array Itmp in makeKernels()



#define KERNEL_N_PAST_MIN_I_LIMIT 2  // n points outside 'signal' to include
#define N_SEEDS 20000                // number of seeds to keep stored
#define MAX_HITS 8192                // max size of hits() array

#define SIGMA_CLIP 3.0               // sigma clipping level to use for stats
#define SUBSAMPLES 1000              // minimum no. of samples to use for stats
#define HIT_WIDTH 20               // dist. betw. 1st and last above for 1 hit

#define FRESFILE_PREFIX "fres-"
#define FRESFILE_PREFIX_TA "fresTA-"
#define FRESFILE_LENGTH 16          // characters in fresfile pattern
#define FRESFILE_LENGTH_TA 18          // characters in fresfile pattern

#define N_CHECKAUTOCORR    100	/* number of xcorr point to use to check
				 and see if autocorr is significant */
#define XCORR_COMPUTE_THRESH 0.4    /* fraction of corrThresh to bother computing */
#define NORM_WIDTH_SEC      20	    /* width in seconds to use for xcorr-norm*/
#define NORMALIZING         0       // to normalize the xcorr values or not
#define XCORR_RMS_WINDOW    6.8     // half-width of the sliding RMS used
                                    // to normalize the xcorr output.
#define XCORR_CORRECTION    1.0    // factor to shift the ideal xcorr value
                                    //   to one consistent with data
#define XCORR_RMS_PERC      0.3   // fraction to scale the XCORR_RMS_WINDOW by
#define SMOOTH_WIDTH_FSU    6.8   /* width of smoothing kern in fsu 2x(peak)*/
#define XCHI_CORRECTION     1.33    // factor to shift the ideal xchi value
#define XCORR_RMS_RANGE     18.0     // rms range for xcor acceptance of a det.

#define PROXIMITY           4        // distance for ok detection
#define FSU_PROX_SAME       3.4      // fsu distance for a detection to be ok
				     // and the distance for a detection to
				     // be the same 'event' as another

#define STRIDE              8

#define SAME_HIT_RANGE      6  /* hits from diff kernels are same hit */

// flags to do very rare things that i might want to do
// #define KERNEL_DUMP         1

// flags
#define RANDOM              1
#define NO_RANDOM           0
#define OVERWRITE           1
#define NO_OVERWRITE        0
#define CENTRE              1
#define NO_CENTRE           0
#define DUMP                1
#define NO_DUMP             0
#define NO_KERNEL           NULL
#define GET_MAX             1
#define GET_MIN             0
#define AUTO_OFFSET         0.0

// labels
#define LAA      "aa"
#define LTEFF     "Teff"
#define LMAXX00  "maxX00"
#define LX00STEP "x00Step"
#define LRSTAR   "RStar"
#define LAU      "AU"
#define LOFFSET  "offset"
#define LANGLE   "angle"
#define LLAMBLO  "lambLo"
#define LLAMBHI  "lambHi"
#define LNSTARS  "nStars"
#define LNLAMBDA "nLambda"
#define LASPRAT  "aspectRatio"
#define LORDER   "order"
#define LELONG   "elong"

// macros

#ifdef DEBUG
#define DPRINTF(...) { fprintf(stderr, __VA_ARGS__); }
#else
#define DPRINTF(...) { }
#endif

#ifdef MK_DEBUG
#define MKDPRINTF(...) { fprintf(stderr, __VA_ARGS__); }
#else
#define MKDPRINTF(...) { }
#endif

#ifdef AKBO_DEBUG
#define AKBODPRINTF(...) { fprintf(stderr, __VA_ARGS__); }
#else
#define AKBODPRINTF(...) { }
#endif

#ifndef DB
#define DB(...) { fprintf(stderr, __VA_ARGS__); }
//#else
//#define DB(...) { }
#endif


#define SQR(x)  ( (x)*(x) )

// global variables and constants


typedef struct {
  int cycles;
  int NperCycle;
  FLOAT corrThresh;
  FLOAT chiThresh;
  char fresDir[MAX_FILENAME];
  int Noffset;
  int norm;
} PARAMLIST;


typedef struct {
  FLOAT x;
  FLOAT I;
} SHADOWCOORD;

typedef struct {
  FLOAT t;
  FLOAT I;
} TSCOORD;


typedef struct {
  char   tsfile[MAX_FILENAME];
  FLOAT elong;
  FLOAT incl;
  FLOAT hz;
  FLOAT dt;
  FLOAT rms;
  FLOAT kurt;
  int n;
  int n2;
  int KOlist[MAX_KNOWN_OBJECTS];
  int nKO;
  TSCOORD *tI;
  FLOAT *fI;
#if defined(USE_FFTW)
    fftwFTYPE_plan fftwplan;
    fftwFTYPE_complex *fftwI;
#else
    gsl_fft_real_workspace_FTYPE *workspace;
    gsl_fft_real_wavetable_FTYPE *Rwavetab;
    gsl_fft_halfcomplex_wavetable_FTYPE *HCwavetab;
#endif
  int FFTisdone;
  int fits;
} TS;

typedef struct {
  FLOAT mean;
  FLOAT rms;
  FLOAT autoval;
  int n;
  FLOAT *value;
} SERIES;


typedef struct {
  char   filename[MAX_FILENAME];
  FLOAT aa;
  FLOAT RStar;
  FLOAT AU;
  FLOAT offset;
  FLOAT lambLo;
  FLOAT lambHi;
  FLOAT x00Step;
  FLOAT maxX00;
  FLOAT maxOffset;
  int nStars;
  int nLambda;
  int order;
  FLOAT Teff;
  FLOAT angle;
  FLOAT aspRat;
  SHADOWCOORD *xI;
  int fullprofile;
  int i_vel;
  FLOAT *I[N_VELOCITIES];   // I
  FLOAT vRet[N_VELOCITIES];
  int nI[N_VELOCITIES];
  int icentre[N_VELOCITIES];
  FLOAT x_start;
  FLOAT t_centre;
  int nRecover[N_VELOCITIES];
  int nHit[N_VELOCITIES];
  FLOAT meanCor[N_VELOCITIES];
  FLOAT meanChi[N_VELOCITIES];
  FLOAT rmsCor[N_VELOCITIES];
  FLOAT rmsChi[N_VELOCITIES];
  FLOAT fsu;
  int ts_eff_length;
  int n;
  int is_all_zero;
  int used_fft;
} SHADOW;


typedef struct {
  FLOAT red_chi2null;
  FLOAT red_chi2;
  FLOAT Pred_chi2;
  FLOAT baryc;
  FLOAT cormag;
  FLOAT chimag;
  FLOAT cormag0;
  FLOAT chimag0;
  int index;
  int n;
  int dof;
  int used_fft;
} HIT; 


void *det_malloc(size_t size);
void *det_calloc(size_t nmemb, size_t size);
void *det_realloc(void *ptr, size_t size);
 
// function declarations
FILE *openfile (char *path, char *mode);
void closefile (FILE *stream);

void loadParameters (char *paramfile, PARAMLIST *pars);
void loadTimeseries (char *tsfile, TS *pts, TS *ptsa,FLOAT elong,FLOAT incl);
void copyTimeseries (TS *pts, TS *ptsa);

int countFresnelFiles (PARAMLIST *pars);
int getFresnelFileList (char *filelist[], PARAMLIST *pars);
void freeFresnelFileList( char *filelist[], int n_file );

void allocateFresPattMem (SHADOW *pfresPatt[], int Npatterns);

void loadFresnelFile (SHADOW *pfresPatt, char *filename);
int loadFresnelFiles (SHADOW *pfresPatt[], PARAMLIST *pars);

void offsetFresPatt (SHADOW *pfresPatt, SHADOW *pfresPatto, FLOAT offset);
int offsetFresPatts (SHADOW *pfresPatt[], PARAMLIST *pars, int Npfres, FLOAT offset);

FLOAT elong2v (FLOAT elong, FLOAT AU, FLOAT dv, int i_vel);
FLOAT elongi2v (FLOAT elong, FLOAT i, FLOAT AU, FLOAT dv, int i_vel);
FLOAT r2D (FLOAT rk, FLOAT elong);
void initPfresPatt (SHADOW *pfresPatt, TS *pts, int i_vel);

int makeKernel (FLOAT *kernel, int *icentre, SHADOW *pfresPatt, TS *pts, int i_vel, int random, int overwrite);

void addKBO (int *addlist, TS *pts, SHADOW *pfresPatt, int n_events, int i_vel, int centering, int random);
void addReadnoise(FLOAT R, TS *ptsa);

void fft_xcorr_norm(FLOAT *corr, FLOAT *corrNorm, TS *ptsa, SHADOW *pfrespatt, FLOAT rms);
void FFTtimeseries(TS *tsa);
//void invFFTtimeseries(FLOAT *ps, FLOAT *corr, int N);
//void FFTkernel(SHADOW *pfresPatt, FLOAT *ps, int N);
int FFTconvolve(FLOAT *corr, TS *ptsa, SHADOW *pfresPatt, int i_vel, int *Neff);
int convolve(FLOAT *corr, TS *ptsa, SHADOW *pfresPatt, int i_vel, int *Neff);
int checkAutocorr(TS *pts, SHADOW *pfresPatt, FLOAT autocorr, FLOAT corrThresh, FLOAT *rms, int n_check);
int xcorrelate (HIT hits[], TS *ptsa, SHADOW *pfresPatt, int i_vel, FLOAT detectionThresh, int dump, int norm);
int xcorrelateKO (HIT hits[], TS *ptsa, SHADOW *pfresPatt, int i_vel, FLOAT detectionThresh, int dump, int norm);
int xchi (HIT hits[], TS *ptsa, SHADOW *pfresPatt, int i_vel, FLOAT detectionThresh, int dump);
int xchi2 (HIT hits[], TS *ptsa, SHADOW *pfresPatt, int i_vel, SERIES *chi);
int xchi3 (HIT hits[], TS *ptsa, SHADOW *pfresPatt, int i_vel, SERIES *chi, HIT corrhits[], int ncorr);

int reducelist(FLOAT *values, int *index, int *nAbove, int n, int sign, int hit_width);

int merge(HIT *list1, HIT *list2, int n1, int n2);
int merge2(HIT *corlist, HIT *chilist, int n1, int n2, SERIES *chi, PARAMLIST *pars);

int countHits(HIT *hits, int nhit);

int addVerify(int *addlist, int nAdd, HIT *corhits, HIT *chihits, int *nHit, HIT *recovered, SHADOW *pfresPatt, TS *pts);

void printHitHeader(FILE *fp);
void printHit(FILE *fp, SHADOW *pfresPatt, HIT corrhits, HIT chihits);

void writeXcorrSeries(TS *ptsa, SHADOW *pfresPatt, FLOAT *corr, FLOAT *corrNorm, HIT hits[MAX_HITS], int n, char *suffix);

void printLenHeader(FILE *fp);
void printLen(FILE *fp, TS *ptsa, SHADOW *pfresPatt);

void printStatsHeader(FILE *fp);
void printStats(FILE *fp, SHADOW *pfresPatt, HIT corrhits, HIT chihits, PARAMLIST *pars, int i_vel);
void printBStatsHeader(FILE *fp);
void printBStats(FILE *fp, SHADOW *pfresPatt, SHADOW *pfresPattA, HIT corrhits, HIT chihits, PARAMLIST *pars, int i_vel);


// could also declare as __attribute__ ((no_return)) to avoid compiler error
void errQuit (char *msg, ...);

void debugOutput (SHADOW *pfresPatt[], PARAMLIST *pars, TS *pts, TS *ptsa, int Npfres);


#endif  // DETECT_LIBRARY
