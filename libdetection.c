/*            */
/* Made with makeScript, Fri Jun 24, 2005  15:58:46 DST */
/* Host: kuiper */
/* Working Directory: /home/bickersj/sandbox/cdiffracSim  */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <dirent.h>
#include <stdarg.h>
#include <errno.h>
#include <string.h>
#include <sys/types.h>
#include <fitsio.h>
#include <assert.h>
#include <time.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

// #define MK_DEBUG 1
// #define AKBO_DEBUG 1

#include "detection.h"
#include "statistics.h"
#include "fft-gsl.h"
#include "fft-fftw.h"

/* ***********************************************
 *
 *  private functions declarations
 *
 * ************************************************ */

DIR *opendirectory (char *path);
void closedirectory (DIR *dir);

/* malloc, calloc, realloc */
void *det_malloc(size_t size) {
  void *tmp = malloc( size );
  assert (tmp != NULL);
  return tmp;
}
void *det_calloc(size_t nmemb, size_t size) {
  void *tmp = calloc(nmemb,size);
  assert(tmp != NULL);
  return tmp;
}
void *det_realloc(void *ptr, size_t size) {
  void *tmp = realloc(ptr, size);
  assert(tmp != NULL);
  return tmp;
}


/* ***********************************************
 *
 *  File handling with error output
 *
 * ********************************************** */

FILE *openfile (char *path, char *mode) {

  FILE *fp;
  if ( (fp = fopen(path, mode)) != NULL) {
    return fp;
  } else {
    errQuit("open file (%s)", path);
  }
  return NULL;  // never reached ... surpresses compiler error
}

void closefile (FILE *stream) {
  if (fclose(stream) != 0) {
    errQuit("close file");
  }
}


DIR *opendirectory (char *path) {

  DIR *dir;
  if ( (dir = opendir(path)) != NULL) {
    return dir;
  } else {
    errQuit("open directory (%s)", path);
  }
  return NULL;  // never reached ... surpresses compiler error
}

void closedirectory (DIR *dir) {
  if (closedir(dir) != 0) {
    errQuit("close directory");
  }
}



/* ***********************************************
 *
 *  Input functions
 *
 * ********************************************** */

void loadParameters (char *paramfile, PARAMLIST *pars) {

  char paramLine[80];
  FILE *fp;
  
  // read in parameters from initialization file
  fp = openfile(paramfile, "r");
  while ( fgets(paramLine, MAX_LINE_LENGTH, fp) != NULL ) {
    
    // get rid of any comment lines or blank lines
    if ( (paramLine[0]=='#') || ( strlen(paramLine)<2 ) ) {
      continue;
    }
    
#if defined (USE_FLOAT)
#define FMT2 "%d %d %f %f %s %d %d"
#else
#define FMT2 "%d %d %le %le %s %d %d"
#endif

    if ( sscanf (paramLine, FMT2, &pars->cycles, &pars->NperCycle, &pars->corrThresh, &pars->chiThresh, pars->fresDir, &pars->Noffset, &pars->norm) != 7 ) {
      errQuit("reading parameter file");
    }

  }
  closefile(fp);

}



void loadTimeseries (char *tsfile, TS *pts, TS *ptsa, FLOAT elong, FLOAT incl) {

  int i=0;
  FLOAT dt=0;
  FLOAT tlast=0;
  FLOAT *I;
  STAT stats;
  int KOi = 0;
  
  pts->tI = det_malloc(1024 * sizeof(TSCOORD));
  I       = det_malloc(1024 * sizeof(FLOAT));

  /*  handle both FITS and text */

  size_t pathlen = strlen(tsfile);
  pts->fits =  ( strstr(tsfile+pathlen-4,"fits") == NULL ) ? 0 : 1;

  if ( pts->fits ) {

    fitsfile *fptr;
    long nrows = 0;
    int status = 0, ncolumns;
    int hdutype, anynul;
    LONGLONG firstrow = 1;
    LONGLONG firstelem = 1;
    FLOAT nulval = NAN;
    float dtf, tlastf;
    int j;
    
    if (! fits_open_file(&fptr, tsfile, READONLY, &status) ) {
      
      if ( fits_movabs_hdu(fptr, 2, &hdutype, &status) ) { /*hdu 2*/
	errQuit("HDU not found\n");
      }
      
      if (hdutype != BINARY_TBL) {
        errQuit("HDU does not contain a binary table.\n");
      } else {
        fits_get_num_rows(fptr, &nrows, &status);
        fits_get_num_cols(fptr, &ncolumns, &status);
	i = (int) nrows;

	if (ncolumns > 1) {
	  errQuit("FITS file contains more than 1 column.");
	}

	/* put the values in arrays */
	pts->tI = det_realloc(pts->tI,sizeof(TSCOORD)*(nrows));
	I = det_realloc(I,sizeof(FLOAT)*(nrows));
	float *tmpf = det_malloc(sizeof(float)*(nrows));


	/* get the dt, tlast values */
	fits_read_key(fptr, TFLOAT, "SAMPTIME", 
		      &dtf, NULL, &status);
	fits_read_key(fptr, TFLOAT, "DURATION", 
		      &tlastf, NULL, &status);
	dt = (FLOAT) (dtf * (i-1));
	tlast = (FLOAT) tlastf;
	
	fits_read_col(fptr, TFLOAT, 1, firstrow, firstelem,
		      nrows, &nulval, tmpf, &anynul, &status);

	/* copy */
	for(j=0;j<nrows;j++) {
	  I[j] = (FLOAT) tmpf[j];
	  pts->tI[j].I = I[j];
	  pts->tI[j].t = (FLOAT) (dtf*j);
	}

	free(tmpf);
	
      }
      fits_close_file(fptr, &status);
    }
    if (status) fits_report_error(stderr, status);
    

  } else {

    char paramLine[80];
    FILE *fp;
    int KOhead = 0;
    int dummyi;
    float dummyf;

    // read in parameters
    fp = openfile(tsfile, "r");
    while ( fgets(paramLine, MAX_LINE_LENGTH, fp) != NULL ) {
      
      //  printf("%d\n",i);
      pts->tI = det_realloc(pts->tI, sizeof(TSCOORD)*(i+10));
      I = det_realloc(I, sizeof(FLOAT)*(i+10));
      
      // deal with any comments
      if (paramLine[0]=='#')  {
	

	// get a list of added objects, if present

	// see if this is the header for added object list
	if ( strstr(paramLine, "Index") != NULL ) {
	  KOhead = 1;
	  continue;
	}
	
	// parse the Known object list
	if ( KOhead ) {
	  if (sscanf(paramLine,"# %d %d %f", &dummyi, &pts->KOlist[KOi++], &dummyf) !=3){
	    errQuit("reading KOlist entry in timeseries file");
	  }
	}
	continue;

      }

      // get rid of blank lines
      if ( strlen(paramLine)<2 ) {
	continue;
      }

#if defined(USE_FLOAT)
#define FMT3 "%f %f"
#else
#define FMT3 "%le %le"
#endif

      if (sscanf(paramLine,FMT3, &pts->tI[i].t, &pts->tI[i].I) !=2){
	errQuit("reading parameter file");
      }

      I[i] = pts->tI[i].I;    
      
      dt += (pts->tI[i].t - tlast);
      tlast = pts->tI[i].t;
      
      i++;
    }
    closefile(fp);
    
  }

  stats     = meanRmsClip_d(I, i, SIGMA_CLIP, STRIDE);
  pts->kurt = kurt_d(I, i, stats.mean, stats.rms);

  free(I);

  strncpy(pts->tsfile, tsfile, MAX_FILENAME);
  pts->dt = dt/(i-1);
  pts->hz = 1.0/pts->dt;
  pts->n = i;
  pts->rms = stats.rms;
  //pts->rms = rmsclip_d(I, i, SUBSAMPLES, SIGMA_CLIP, 0, 3);
  pts->elong = elong;
  pts->incl = incl;
  pts->FFTisdone = 0;
  
  pts->nKO = KOi;

  /* define n2 to pad up to the nearest power of 2 */
  double base2d, base2rem;
  base2rem = modf( log2(pts->n), &base2d);
  int base2 = (int) ( base2d + (base2rem?1:0) );
  pts->n2 = (int) ldexp(1,base2);

  copyTimeseries(pts, ptsa);



}


void copyTimeseries (TS *pts, TS *ptsa) {

  ptsa->tI = det_malloc(sizeof(TSCOORD)*pts->n);

  strcpy(ptsa->tsfile, pts->tsfile);
  ptsa->fits      = pts->fits;
  ptsa->dt        = pts->dt;
  ptsa->hz        = pts->hz;
  ptsa->n         = pts->n;
  ptsa->rms       = pts->rms;
  ptsa->kurt      = pts->kurt;
  ptsa->elong     = pts->elong;
  ptsa->incl      = pts->incl;
  ptsa->FFTisdone = pts->FFTisdone;
  ptsa->n2        = pts->n2;
  ptsa->nKO       = pts->nKO;

  for(int i=0;i<pts->n;i++) {
    ptsa->tI[i] = pts->tI[i];
  }

  for(int i=0;i<pts->nKO;i++) {
    ptsa->KOlist[i] = pts->KOlist[i];
  }

}




int countFresnelFiles (PARAMLIST *pars) {

  DIR *dir;
  struct dirent *dp;
  int i,filecount = 0;
  int fresfilecount = 0;
  char *filelist[MAX_FRESFILES], *filename;

  // get the contents of the fresnelFiles directory
  dir = opendirectory(pars->fresDir);
  while ( ( dp = readdir(dir)) != NULL ) {
    
    filelist[filecount] = det_malloc(MAX_FILENAME*sizeof(char));
    filelist[filecount] = strncpy(filelist[filecount], dp->d_name, 
				  MAX_FILENAME);
    filecount++;
    
  }
  closedirectory(dir);

  for (i=0; i<filecount; i++) {
    filename = filelist[i];
    
    // ignore files with the wrong filenames
    if ( (strstr(filename, FRESFILE_PREFIX) != filename)  ||  
	 (strlen(filename) != FRESFILE_LENGTH) ) {
      continue;
    }
    fresfilecount++;
  }
  return fresfilecount;
}



int getFresnelFileList (char *filelist[], PARAMLIST *pars) {

  DIR *dir;
  struct dirent *dp;
  int filecount = 0;
  char path[MAX_FILENAME];

  // get the contents of the fresnelFiles directory
  dir = opendirectory(pars->fresDir);
  while ( ( dp = readdir(dir)) != NULL ) {

    // ignore files with the wrong filenames
    if ( ( (strstr(dp->d_name, FRESFILE_PREFIX) != dp->d_name)  ||
	   (strlen(dp->d_name) != FRESFILE_LENGTH) ) &&
	 ( (strstr(dp->d_name, FRESFILE_PREFIX_TA) != dp->d_name)  ||
	   (strlen(dp->d_name) != FRESFILE_LENGTH_TA) )  ) {
      continue;
    }

    filelist[filecount] = det_malloc(MAX_FILENAME*sizeof(char));
    filelist[filecount] = strncpy(filelist[filecount], dp->d_name, 
				  MAX_FILENAME);
    snprintf (path, MAX_FILENAME,"%s/%s", pars->fresDir, dp->d_name);
    strcpy(filelist[filecount], path);

    filecount++;

  }
  closedirectory(dir);

  return filecount;
}


void freeFresnelFileList( char *filelist[], int n_file ) {
  for (int i_f=0; i_f<n_file; i_f++) {
    free(filelist[i_f]);
  }
}

void allocateFresPattMem (SHADOW *pfresPatt[], int Npatterns) {

  for (int j=0; j<Npatterns; j++) {
    pfresPatt[j]     = det_malloc(sizeof(SHADOW));
    
    // will be realloc'd when read in.
    pfresPatt[j]->xI = det_malloc(10*sizeof(SHADOWCOORD)); 
  }
}



void loadFresnelFile (SHADOW *pfresPatt, char *filename) {

  char paramLine[MAX_LINE_LENGTH]; // s[MAX_LINE_LENGTH];
  int j;
  FLOAT maxOffset;
  int offsetcount;   // go this many points past a specified I value
  FILE *fp;

  // s1='#'    s2=keyword    s3=comment   s4=value
  char s1[2], s2[32], s3[32], s4[32]; 

  FLOAT x, Ihole, Idisk;

  char *offset_min_I_limit_tmp = getenv("MIN_SIGNAL");
  FLOAT offset_min_I_limit = offset_min_I_limit_tmp != NULL ? 
    atof(offset_min_I_limit_tmp) :	OFFSET_MIN_I_LIMIT;
  
  // read in values from file
  fp = openfile(filename, "r");
  j=0;
  maxOffset = 0;
  offsetcount = 0;

  while ( fgets(paramLine, MAX_LINE_LENGTH, fp) != NULL ) {
    
    // ignore blank lines
    if ( strlen(paramLine)<2 )
      continue;
    
    // get the parameters from the header
    if ( paramLine[0]=='#') {
      
      // ignore other comments
      if ( sscanf (  paramLine, "%s %s %s %s", s1, s2, s3, s4) !=4 ) {
	continue;
      
      // match keywords to read in the header parameters
      } else {
	if (      strcmp(LAA,s2)        == 0 ) {
	  pfresPatt->aa      = atof(s4);
	}else if ( strcmp(LMAXX00,s2)    == 0 ) {
	  pfresPatt->maxX00  = atof(s4);
	}else if ( strcmp(LX00STEP,s2)   == 0 ) {
	  pfresPatt->x00Step = atof(s4);
	}else if ( strcmp(LRSTAR,s2)     == 0 ) {
	  pfresPatt->RStar   = atof(s4);
	}else if ( strcmp(LAU,s2)        == 0 ) {
	  pfresPatt->AU      = atof(s4);
	}else if ( strcmp(LOFFSET,s2)    == 0 ) {
	  pfresPatt->offset  = atof(s4);

	}else if ( strcmp(LNSTARS,s2)    == 0 ) {
	  pfresPatt->nStars  = atof(s4);
	}else if ( strcmp(LNLAMBDA,s2)   == 0 ) {
	  pfresPatt->nLambda = atof(s4);
	}else if ( strcmp(LASPRAT,s2)    == 0 ) {
	  pfresPatt->aspRat  = atof(s4);
	}else if ( strcmp(LORDER,s2)     == 0 ) {
	  pfresPatt->order   = atof(s4);
	}else if ( strcmp(LLAMBLO,s2)    == 0 ) {
	  pfresPatt->lambLo  = atof(s4);
	}else if ( strcmp(LLAMBHI,s2)    == 0 ) {
	  pfresPatt->lambHi  = atof(s4);
	}else if ( strcmp(LTEFF,s2)    == 0 ) {
	  pfresPatt->Teff  = atof(s4);
	}else if ( strcmp(LANGLE,s2)    == 0 ) {
	  pfresPatt->angle  = atof(s4);
	}
      }

      
      // get the data
    } else {  // ie. if not a comment 

#if defined(USE_FLOAT)
#define FMT "%f %f %f"
#else
#define FMT "%le %le %le"
#endif 

      if ( sscanf (  paramLine, FMT, &x, &Ihole, &Idisk ) != 3 ) {
	errQuit("reading fresnel file");
      }

      pfresPatt->xI = 
	det_realloc(pfresPatt->xI, (2*j + 1) * sizeof(SHADOWCOORD));
      pfresPatt->xI[j].x = x; 
      pfresPatt->xI[j].I = Idisk;
      
      
      // also figure out where the pattern fades to ~zero
      
      // increment the counter if Idisk is low, reset it otherwise;
      if ( fabs(Idisk-1.0) <  offset_min_I_limit )  {
	offsetcount++;

      // else tricky .. decrement to zero.
      //   I'd just set it to zero, but then a
      //   stray large value could reset
      } else  {
	offsetcount = (offsetcount > 0) ? ( offsetcount - 1 ) : 0;
      }
      
      // this will keep being set until offsetcount grows
      // larger than N_PAST..
      if ( offsetcount < N_PAST_MIN_I_LIMIT ) {
	maxOffset = x;
      }
      
      j++;  // counter for number of values read in
      
    } // end of if (comment) {} else {}
    
  } // end while (fgets paramLine) ...    

  // if there are negatives, we have the full profile, not just half
  pfresPatt->fullprofile = (pfresPatt->xI[0].x < 0) ? 1 : 0;

  snprintf (pfresPatt->filename, MAX_FILENAME,"%s", filename);
  pfresPatt->is_all_zero = (offsetcount == j) ? 1 : 0;
  pfresPatt->maxOffset = (offsetcount == j) ? x/2 : maxOffset;

  /* note: this changes the size of teh conv kernel */
  FLOAT ulamb = 0.5 * ( pfresPatt->lambLo + pfresPatt->lambHi );
  FLOAT Dist = AU_M * pfresPatt->AU;
  FLOAT fsu = sqrt( Dist * ulamb / 2.0 );
  pfresPatt->maxOffset = pfresPatt->RStar + pfresPatt->aa + 2.5*fsu;

  pfresPatt->n = j;
  pfresPatt->fsu = fsu;
  
  closefile(fp);
}


int loadFresnelFiles (SHADOW *pfresPatt[], PARAMLIST *pars) {

  int filecount = 0;
  int fresfilecount = 0;
  char *filelist[MAX_FRESFILES], *filename, path[MAX_FILENAME];
  DIR *dir;
  struct dirent *dp;

  // get the contents of the fresnelFiles directory
  dir = opendirectory(pars->fresDir);
  while ( ( dp = readdir(dir)) != NULL ) {
    filelist[filecount] = det_malloc(sizeof(char));
    filelist[filecount] = det_malloc(MAX_FILENAME*sizeof(char));
    strcpy(filelist[filecount++], dp->d_name);
  }
  closedirectory(dir);


  // Go through each file and read in the data
  for (int i=0; i<filecount; i++) {
    filename = filelist[i];
    
    // ignore files with the wrong filenames
    if ( ( (strstr(filename, FRESFILE_PREFIX) != filename)  ||
	   (strlen(filename) != FRESFILE_LENGTH) ) &&
	 ( (strstr(filename, FRESFILE_PREFIX_TA) != filename)  ||
	   (strlen(filename) != FRESFILE_LENGTH_TA) )  ) {
      continue;
    }
    snprintf (path, MAX_FILENAME,"%s/%s", pars->fresDir, filename);
    loadFresnelFile(pfresPatt[fresfilecount], path);
    fresfilecount++;
    //free(filelist[i]);
  }

  return fresfilecount;
}




/* offset the pattern pointed to by pfresPatt an amount 'offset'  *
 * and put the result in the location pointed to by pfresPatto */
void offsetFresPatt (SHADOW *pfresPatt, SHADOW *pfresPatto, FLOAT offset) {

  int i_lo, i_hi;
  double r_rem, r_int;  // fractional remainder and integer part of r
  FLOAT r, x, x00Step;
  FLOAT I, I_lo, I_hi;

  x00Step = pfresPatt->x00Step;

  pfresPatto->xI = 
    det_realloc(pfresPatto->xI,sizeof(SHADOWCOORD)*pfresPatt->n);
  
  // loop over all x coordinates  at this offset
  int const nfres = pfresPatt->n;
  for (int k=0; k<nfres; k++) {
    
    x = k*x00Step;
    r = sqrt(x*x + offset*offset);     // dist from cent of shadow
    r_rem = modf( r/x00Step, &r_int);  
    
    i_lo = (int) r_int;       // indices of array to interp betw
    i_hi = i_lo + 1;
    
    // intensity values to interpolate between
    I_lo = (i_hi<pfresPatt->n) ? pfresPatt->xI[i_lo].I : 1.0;
    I_hi = (i_hi<pfresPatt->n) ? pfresPatt->xI[i_hi].I : 1.0;
    
    // intensity at this x,offset
    I = I_lo + (r_rem)*(I_hi-I_lo)/x00Step;

    pfresPatto->xI[k].x = x;     // stuff the values in the array
    pfresPatto->xI[k].I = I;
    
  }  // end loop over x
  
  // make the entries in the structure (copies except for offset)
  pfresPatto->aa       = pfresPatt->aa;
  pfresPatto->RStar    = pfresPatt->RStar;
  pfresPatto->AU       = pfresPatt->AU;
  pfresPatto->offset   = offset;
  pfresPatto->x00Step  = pfresPatt->x00Step;
  pfresPatto->maxX00   = pfresPatt->maxX00;
  pfresPatto->n        = pfresPatt->n;
  pfresPatto->lambLo    = pfresPatt->lambLo;
  pfresPatto->lambHi    = pfresPatt->lambHi;
  pfresPatto->nLambda   = pfresPatt->nLambda;
  pfresPatto->nStars    = pfresPatt->nStars;
  pfresPatto->maxOffset = pfresPatt->maxOffset;
  pfresPatto->order     = pfresPatt->order;
  pfresPatto->aspRat    = pfresPatt->aspRat;
  pfresPatto->fullprofile = pfresPatt->fullprofile;
  pfresPatto->is_all_zero = pfresPatt->is_all_zero;
  pfresPatto->fsu       = pfresPatt->fsu;
  strcpy (pfresPatto->filename, pfresPatt->filename);
  
}




int offsetFresPatts (SHADOW *pfresPatt[], PARAMLIST *pars, int Npfres, FLOAT offset) {

  int N=Npfres;
  int Noffset;
  FLOAT ostep, offset_tmp=0;   // ostep = offset step

  /* if you specify an offset, just use that one
   otherwise, just use the default number */
  Noffset = (offset > 0) ? (1) : (pars->Noffset);

  // loop over all fresnel patterns
  /* the idea here is that the zero-offset patterns are in the first Npfres
   elements of the array, and the offset ones are placed in subsequent
   elements.  In practice, I only load one ata  time, so Npfres=1 and
   and the for loop is only executed once */
  if (Noffset > 0) {
    for (int i=0; i<Npfres; i++) {
      
      ostep = pfresPatt[i]->maxOffset / Noffset;  // offset step
      
      // loop over all the offsets
      if (!pfresPatt[i]->is_all_zero) { // don't bother for 'nosignal'
	for (int j=0; j<Noffset; j++) {
	  
	  offset_tmp = (offset) ? (offset) : ( (j+1)*ostep ) ;
	  offsetFresPatt(pfresPatt[i], pfresPatt[N], offset_tmp);
	  N++;
	  
	} // end loop over offsets
      } else {
	for (int j=1; j<Noffset; j++) {
	  pfresPatt[j]->is_all_zero = 1;
	}
      }

    }  // end loop over fresnel patterns
  }
    
  return N;
}



FLOAT elong2v (FLOAT elong, FLOAT rk, FLOAT dv, int i_vel) {

  FLOAT v;
  FLOAT Ve = sqrt(G*Mo/Re);

  FLOAT dvcoef = (i_vel%2) ? (-(i_vel+1)/2) : (i_vel/2);  // if odd

  v = Ve * ( 
	    sqrt( (1.0/rk) * (1.0-(1.0/(rk*rk))*SQR(sin(elong))) ) + 
	    cos(elong) 
	     );
  v += dvcoef*dv;  // now depricated, kept for backward compat.

  return v;
}



FLOAT elongi2v (FLOAT elong, FLOAT i, FLOAT rk, FLOAT dv, int i_vel) {

  // elong and i in radians, rk in AU

  FLOAT dvcoef = (i_vel%2) ? (-(i_vel+1)/2) : (i_vel/2);  // if odd

  FLOAT Ve = sqrt(G*Mo/Re);           // orb vel of earth (m/s)

  FLOAT rkp       = rk * cos(i);      // rk projected down
  FLOAT elong_lim = (rkp<1.0) ? asin(rkp) : 360.0; // elong limit for i and rk
  FLOAT elongD;                       // elong in degrees

  // get the dir cosines of the line connecting the earth to kbo
  
  if (elong > elong_lim) {
    elongD = DEG*elong_lim;
    fprintf (stderr, "Elongation cannot exceed %.2f deg at this inclination.  Returning 0\n", elongD);
    return 0;
  }
  
  FLOAT zk = rk * sin(i);          // z coord of the kbo
  FLOAT dp = r2D(rkp, elong);      // dist 'd' proj'd down on eclip
  FLOAT d  = sqrt(dp*dp + zk*zk);  // dist 'd' from earth to kbo


  // this whole process comes from Schaum's math. methods pp. 31.
  //  But it's exactly the same as taking the dot-produce betw the 
  //   relative velocity vector and the relative position vector.
  
  // get the direction cosines of the relative velocity
  FLOAT vk = Ve / sqrt(rk);
  FLOAT vex = 0.0;
  FLOAT vey = Ve;
  FLOAT vez = 0.0;
  
  FLOAT beta = asin ( sin(elong) / rkp );
  FLOAT alpha = PI - elong;

  // get the coords of the earth
  FLOAT xe = 1.0;
  FLOAT ye = 0.0;
  FLOAT ze = 0.0;

  // coords of the kbo (zk is above)
  FLOAT xk = rkp*cos(alpha-beta);
  FLOAT yk = rkp*sin(alpha-beta);

  // direction cosines of vectro from earth to kbo
  FLOAT l = (xk-xe)/d;
  FLOAT m = (yk-ye)/d;
  FLOAT n = (zk-ze)/d;

  // velocity vector of kbo when at highest ecliptic latitude
  FLOAT vkx = -vk*sin(alpha-beta);
  FLOAT vky = vk*cos(alpha-beta);
  FLOAT vkz = 0.0;
  
  // relative velocity wrt earth
  FLOAT dvx = vkx - vex;
  FLOAT dvy = vky - vey;
  FLOAT dvz = vkz - vez;
  FLOAT dv_rel = sqrt(dvx*dvx + dvy*dvy + dvz*dvz);

  // direction cosine of relative vel wrt earth
  FLOAT l2 = dvx/dv_rel;
  FLOAT m2 = dvy/dv_rel;
  FLOAT n2 = dvz/dv_rel;

  // angle betw kbo position vec (wrt earth), and vel vec (wrt earth)
  FLOAT theta = acos(l*l2 + m*m2 + n*n2);

  // proj of relative vel onto plane which is perp to line-of-sight
  FLOAT vPerp = dv_rel * sin(theta);

  vPerp += dvcoef * dv;  // now depricated, kept for backward compat

  // this is a kluge.  if you set the elongation to be <0.1deg
  //    then the incl is used as the velocity
  //   ie. you can set the velocity explicitly
  //   I needed to do this to test the chang06 result obtained
  //   with satellite data.  vRet values were > Vearth and elongi2v()
  //   was meaningless.
  
  if (elong < RAD*0.1) {
    vPerp = DEG*i;   // need DEG because it was *RAD when read in
  }

  return vPerp;  // in metres
}


// calculate the dist to an object given it's orb radius and elong
//  see page 182 of my notes.  ... it's just cosine law.
FLOAT r2D (FLOAT rk, FLOAT elong) {

    FLOAT a = Re;
    FLOAT b = rk*Re;
    FLOAT alpha = PI - elong;

    FLOAT beta = asin( sin(elong)/rk );
    FLOAT theta = alpha - beta;
    
    return sqrt(a*a + b*b - 2.0 * a * b * cos(theta) )/Re;
}

// original elong2v program without inclination
/*
FLOAT elong2v (FLOAT elong, FLOAT AU, FLOAT dv, int i_vel) {

  FLOAT sin_alpha;
  FLOAT Vk;        // kbo velocity (assumed circular orbit);
  FLOAT Ve;        // earth orbit velocity (again circular);
  FLOAT beta;      // angle betw kbo and lin of sight (sun at vertex)
  FLOAT alpha;     // angle betw earth and lin of sight (sun at vert)
  FLOAT V_k_perp;  // component of Vk perp to l.o.s.
  FLOAT V_e_perp;  // component of Ve perp to l.o.s.
  FLOAT dvcoef;    // coefficient to control sign of velocity diff dv
  FLOAT vRet;      // local version of v-retrograd

  // odd indices for negatives, even for positives
  dvcoef = (i_vel%2) ? (-(i_vel+1)/2) : (i_vel/2);  // if odd
  
  Ve = sqrt(G*Mo/Re);
  alpha = PI - elong;
  V_e_perp = Ve * cos(alpha);
  sin_alpha = sin(alpha);  // self-explanitory
  
  Vk = Ve / sqrt(AU);
  beta = asin ( sin_alpha/AU );
  V_k_perp = Vk * cos(beta);
  vRet = V_e_perp - V_k_perp + dvcoef*dv;
  
  return vRet;
  
}
*/


void initPfresPatt (SHADOW *pfresPatt, TS *pts, int i_vel) {

  int size;
  FLOAT vRet;

  // get and store vRet
  // vRet = elong2v(pts->elong, pfresPatt->AU, pfresPatt->dv, i_vel);
  vRet = elongi2v(pts->elong, pts->incl, pfresPatt->AU, 0, i_vel);
  pfresPatt->vRet[i_vel] = fabs(vRet);

  // allocate memory for the kernel
  size = 2 * (int) trunc(pfresPatt->maxX00/
			 (pfresPatt->vRet[i_vel]*pts->dt)) + 6;
			 
  pfresPatt->I[i_vel] = det_malloc(size * sizeof(FLOAT));
}




// new makeKernel
int makeKernel (FLOAT *kernel, int *icentre, SHADOW *pfresPatt, TS *pts, int i_vel, int random, int overwrite) {

  int i,k,l;    // loop indices

  FLOAT x1, x2, x, x_last, x_next, I1, I2, I, I_last;
  FLOAT dx, xmx1, xmx2, xmx_end;  // x intervals used in integration
  FLOAT t_start;                  // begin of kernel time
  FLOAT x_start, x_data, dx_sample, x_end;  // begin kernel position
  const FLOAT x00Step = pfresPatt->x00Step; // local version

  const FLOAT dt = pts->dt;    // sampling time
  FLOAT t_step;                // x00Step / vRet
  FLOAT t_pad;                 // the time for one padding unit
 
  int imax;   // highest index in x-position array for pattern
  FLOAT tmax; // highest time at imax position

  int icentre_tmp; // store the index of the centre of the pattern
  FLOAT t_centre;  // the time coord of the true pattern centre

  FLOAT vRet;      // retrograd velocity

  int nt, nx;      // number of values in array converted to 'times'
  int ntPatt_pad;  // number of 1's pad the start of the pattern with
                   // ... this allows the 'random' dt offset

  FLOAT Isum;          // store integrating sum

  int intcount;        // num of points in integrated series
  int cullcount;       // num of points left after zeros culled
  int i_start, i_stop; // flags to cull zeros from the integr. series

  char *kernel_lo_I_cut_tmp = getenv("MIN_SIGNAL");
  FLOAT kernel_lo_I_cut = (kernel_lo_I_cut_tmp != NULL) ?
    atof(kernel_lo_I_cut_tmp) : KERNEL_LO_I_CUT;
  
  MKDPRINTF("makeKernel(): Before initializations\n");

  // initialize values
  pfresPatt->i_vel = i_vel;

  vRet = pfresPatt->vRet[i_vel];
  t_step = x00Step / vRet;      // time between x00Steps at vRet
                                // ~ 2 ms for 50m at 26km/s
  imax = pfresPatt->n - 1;              // highest index
  tmax = pfresPatt->xI[imax].x / vRet;  // highest time
  nt = 2 * (int) ( trunc(tmax / dt) + 2 );    // num spaces nec. for data

  //pad with enough 1's for dt step
  // this is important for undersampled data
  //  12 for dt=0.025s and t_step=0.002s
  ntPatt_pad = (int) ( trunc(dt / t_step) + 1 );  
  nx = 2 * pfresPatt->n + 3*ntPatt_pad; // 1pad at start, 2 at end
  icentre_tmp = (pfresPatt->fullprofile) ? 
    pfresPatt->n/2  : ntPatt_pad + pfresPatt->n;
  t_pad = ntPatt_pad * t_step;          // time added due to padding
  t_centre = tmax + t_pad;

  // decl and init Iint (the integrals of x*dt ... the kernel!) and 
  // xItmp (the x-values from the templ which are to be integrated).
  FLOAT Iint[nt+2];
  for (i=0; i<(nt+2); i++) {
    Iint[i] = 1.0;
  }
  FLOAT xItmp[nx+10];  // add 10 extra as a buffer
  for (i=0; i<nx; i++) {
    xItmp[i] = 1.0;
  }


  //FLOAT intTest = 0;
  //FLOAT intTest1 = 0;
  //FLOAT intTest2 = 0;
  // do a straight integration  (used only to test output)
  //for (i=1; i<=imax; i++) {
  //  intTest += (pfresPatt->xI[i].I - 1.0);
  //}
  //intTest = 2*intTest + (pfresPatt->xI[0].I - 1.0);
  //intTest *= t_step;


  // copy the half of the pattern in order to be a whole pattern
  //   start in the middle and fill the array outward
  //  Modefied 2007-11-
  //  --> if ->fullprofile, then the pattern must be for
  //  an entire (asymmetric) profile ... don't need to x2 the pattern
  int const nfres = pfresPatt->n;
  if (!pfresPatt->fullprofile) {
    for (i=0; i<nfres; i++) {
      xItmp[ntPatt_pad + nfres + i] = 
	(i) ? pfresPatt->xI[i].I  : (pfresPatt->xI[0].I);
      xItmp[ntPatt_pad + nfres - i] = 
	(i) ? pfresPatt->xI[i].I  : (pfresPatt->xI[0].I);
    }
  } else {
    for (i=0; i<nfres; i++) {
      xItmp[i] = pfresPatt->xI[i].I;
    }
  }

  
  MKDPRINTF("makeKernel():  Before rand() call\n");

  // set the non-random start-time
  //  a=(tmax + t_pad + t_step) is the time from 0 to the centre
  //  b=(trunc (tmax/dt))*dt is the longest time evenly div by dt
  //   ... t_start = a-b means that the integ'n will have a point at
  //       the centre of the pattern

  t_start = (tmax + t_pad + t_step) - ( trunc(tmax/dt) )*dt; 
  //  add on a random offset from 0 to 'dt'
  if ( random )  {
    t_start += ( (FLOAT) rand()/RAND_MAX )*dt;
  }

  ///////////////////////////////////////////////////////
  // go to each time position and figure out the intensity
  //   by interpolating from the x-position array

  dx_sample = dt*vRet;                  // 
  x_data = 2 * pfresPatt->xI[imax].x;   // the highest x value
  x_start = t_start * vRet;             // the starting x value
  pfresPatt->x_start = x_start;
  x1 = x_start;                         
  x2 = x1 + dx_sample;
  x_end = x_start + ( trunc(x_data/dx_sample) + 2 )*dx_sample;


  Isum = 0;
  l = 0;

  // x, x_last   ... coords of values in array
  // x1, x2      ... coords of 1st,last points in integral range


  MKDPRINTF("makeKernel():  Before k-loop\n");
  for (k=0; k<nx; k++) {

    //intTest1 += (xItmp[k]-1.0) * t_step;

    //MKDPRINTF("makeKernel(): k=%d\n",k);
    x = k * x00Step;
    x_last = x - x00Step;
    x_next = x + x00Step;
    
    I = xItmp[k];
    I_last = (k>0) ? xItmp[k-1] : 1.0;
    //I_next = (k<(nx-1)) ? xItmp[k+1] : 1.0;
    
    xmx1 = fabs(x-x1);
    xmx2 = fabs(x-x2);
    xmx_end = fabs(x-x_end);
    
    if (k==icentre_tmp) {
      *icentre = l;
      t_centre = (x-x1)/(x2-x1);
    }


    //if (x00Step > MAX_KERNELOVERSAMPLING*vRet*dt)
    //  errQuit("Kernel too finely sampled for array size ... exiting\n(Try increasing MAX_KERNEL_OVERSAMPLING)\n");

    // if we're sampling courser than templ pattern ... integrate
    if ( x00Step< vRet*dt ) {

      // start integrating ... first interval
      if (   x>=x1   &&   xmx1 < x00Step    &&   x<x2  ) {
	dx = x-x1;
	I1 = I_last + (I - I_last)*(dx/x00Step);
	
	Isum = 0.5*(I1+I)*dx/vRet;
	
	// if ( between x1 and x2 )
      } else if ( x>=x1   &&   x<x2 ) {
	dx = x-x_last;
	Isum += 0.5*(I_last+I)*dx/vRet;
	
	// if ( past x2 but not the end ) and
      } else if ( x>=x2   &&  x<x_end) {
	
	// finish this interval
	dx = x2 - x_last;
	I2 = I_last + (I - I_last)*(dx/x00Step);
	Isum += 0.5*(I_last+I2)*dx/vRet;
	
	Iint[l++] = Isum/dt;  // the average value 

	
	// change the endpoint values
	x1 = x2;
	x2 += dx_sample;
	I1 = I2;

	// start the next interval going
	dx = x-x1;
	Isum = 0.5*(I1+I)*dx/vRet;
	
	
	// for the last point ... just tie things up.
      } else if ( x>=x_end   && xmx_end<=x00Step) {
	
	dx = x_end - x_last;
	I2 = I_last + (I - I_last)*(dx/x00Step);
	Isum += 0.5*(I_last+I2)*dx/vRet;

	Iint[l++] = Isum/dt; 

      }


    // what if we're sampling finer than the template pattern 
    //    we have to interpolate to integrate
    } else {

      // start with the partial integration step      
      dx = x2-x_last;
      I2 = I_last + (I-I_last)*(x2-x_last)/x00Step;
      if (x_last>x1) {  // don't do this on the first segment
	Isum += 0.5*(I2+I_last)*dx/vRet;
	if (l<(nt+2)) {
	  Iint[l++] = Isum/(dt);
	}
	Isum = 0;

      }
      x1 = x2;
      x2 = x1 + vRet*dt;
      

      // interpolate the full steps in the middle
      while (x>x2) {
	
	I1 = I_last + (I-I_last)*( (x1-x_last)/x00Step );
	I2 = I_last + (I-I_last)*( (x2-x_last)/x00Step );
	Isum = 0.5*(I1+I2)*dt;

	if (l<(nt+2)) {
	  Iint[l++] = Isum/dt;
	}

	x1 = x2;
	x2 = x1 + vRet*dt;
      }
	
      // finish-off with the partial step at the end.
      dx = x-x1;
      I1 = I2;
      Isum = 0.5*(I1+I)*dx/vRet;

    }
    
    
  } // end loop over all times (k)
    ////////////////////////////////////////////////////////////
  intcount = l;
  
  
  ///////////////////////////////////////////////////////////
  // now keep only the non-zero values (below kernel_lo_I_cut).
  MKDPRINTF("makeKernel(): Before culling out zeros from kernel\n");
  cullcount = 0;
  if (kernel != NULL) {
    kernel[cullcount] = 1.0;
  }
  if (overwrite) {
    pfresPatt->I[i_vel][cullcount] = 1.0;   // pad w init 'nosignal'
  }
  cullcount++;


  i_start = 0;
  i_stop = 0;

  // Flag first and last points w enough signal to bother keeping
  for (k=0; k<intcount; k++) {

    if ( fabs(Iint[k] - 1.0)  > kernel_lo_I_cut   &&  !i_start ) {
      i_start = k;
    }
    if ( fabs(Iint[intcount-k-1] - 1.0) > kernel_lo_I_cut &&  
	 !i_stop ) {
      i_stop = intcount - k - 1;
    }
  }
    
  // save the values that we've chosen to keep
  icentre_tmp = *icentre;
  for (k=i_start-KERNEL_N_PAST_MIN_I_LIMIT; k<i_stop+KERNEL_N_PAST_MIN_I_LIMIT; k++) {

    if (k<0 || k>intcount) {
      continue;
    }
    
    if (k == icentre_tmp) {
      *icentre = cullcount;
    }

    //intTest2 += (Iint[k] - 1.0)*dt;
    
    // only overwrite the existing pattern if asked to
    if (overwrite) {
      pfresPatt->I[i_vel][cullcount] = Iint[k];
    }
    if (kernel != NULL)  {
      kernel[cullcount] = Iint[k];
    }
    cullcount++;   
    
  }
  
  if (overwrite) {
    pfresPatt->I[i_vel][cullcount] = 1.0;  // tail with 'nosignal'
  }
  if (kernel != NULL)  {
    kernel[cullcount] = 1.0;
  }
  cullcount++;
  
  
  // set the new values for point number and centres
  if (overwrite) {
    pfresPatt->nI[i_vel] = cullcount;
    pfresPatt->icentre[i_vel] = *icentre;


  // if there were no points worth keeping set the all_zero flag
    if ( pfresPatt->nI[i_vel] < pfresPatt->icentre[i_vel] ) {
      *icentre = pfresPatt->nI[i_vel] / 2;
      pfresPatt->icentre[i_vel] = *icentre;
    }
  }
  ///////////////////////////////////////////////////////////////

  t_centre = (t_centre + (FLOAT) *icentre - 1.0)*dt;
  pfresPatt->t_centre = t_centre;

  MKDPRINTF("makeKernel(): Before return\n");
  return cullcount;   // the number of points in the kernel
}






void addKBO (int *addlist, TS *ptsa, SHADOW *pfresPatt, int n_events, int i_vel, int centering, int random) {

  int i, j, n, icentre;

  FLOAT vRet = pfresPatt->vRet[i_vel];
  int size = 2* (int) ( trunc(pfresPatt->maxX00/(vRet*ptsa->dt)) + 6);
  FLOAT kernel[size];

  int istart;

  for (i=0; i<n_events; i++) {
    
    // args (kernel, pfresPatt, ptsa, i_vel, random, overwrite)
    AKBODPRINTF ("addKBO(): Before makeKernel()\n");
    n = makeKernel(kernel, &icentre, pfresPatt, ptsa, i_vel, random, NO_OVERWRITE);

    // don't want points inserted right at the beginning or 
    //   end of the timeseries
    AKBODPRINTF ("addKBO(): Before rand()\n");
    
    istart = 0;
    if (centering) {
      istart =  (ptsa->n - n) / 2;
    } else {
      while (istart<n || istart>(ptsa->n - n) )
	istart = rand() % ptsa->n;
    }
    
    AKBODPRINTF("addKBO(): Before kernel insertion into data\n");
    // scale the intensity to add the event
    //   note that I=1.000000 are masked and won't be scaled
    for (j=0; j<n; j++ ) {
      if ( (istart+j) < ptsa->n   &&   
	   (istart+j) >= 0  &&  
	   ptsa->tI[istart+j].I != 1.000000 ) {
	ptsa->tI[istart+j].I *= kernel[j];
      }
    }

#if defined(KERNEL_DUMP)
    FILE *fp;
    char kernelfile[8];
    snprintf (kernelfile, 8, "kern%02d", i);
    fp = openfile(kernelfile, "w");
    for (j=0; j<n; j++) {
      fprintf (fp, "%.6f %.6f\n", ptsa->tI[istart+j].t, kernel[j]);
    }
    closefile(fp);
#endif

    addlist[i] = istart + icentre;

  }
}


void addReadnoise(FLOAT R, TS *ptsa) {

  const FLOAT counts = 1.0/(SQR(ptsa->rms));
  const FLOAT rms_rdnoise = R/counts;

  /* get normal variates and add them */
  const gsl_rng_type *T;
  gsl_rng *r;
  gsl_rng_env_setup();
  T = gsl_rng_default;
  r = gsl_rng_alloc(T);
  gsl_rng_set(r, time(NULL));
  for (int i=0;i<ptsa->n;i++) {
    FLOAT rdnoise = gsl_ran_gaussian(r, rms_rdnoise);
    ptsa->tI[i].I += rdnoise;
    //if (ptsa->tI[i].I < 0) {
    //  ptsa->tI[i].I = 0;
    //}
  }

}

void smooth(FLOAT *corr, FLOAT *corr_sm, int n, int Nsmth) {

  /* pre-compute the gaussian values */
  const FLOAT Nsig = 3.0;
  int mid = (int) Nsig * Nsmth;
  FLOAT gauss[2*mid];
  gauss[mid] = 1.0 / (sqrt(TWOPI)*Nsmth);
  FLOAT gsum = gauss[mid];
  for(int j=1; j<mid; j++) {
    gauss[mid+j] = (1.0/(sqrt(TWOPI)*Nsmth)) * 
      exp(- SQR(j)/(2.0*SQR(Nsmth)));
    gauss[mid-j] = gauss[mid+j];
    gsum += gauss[mid+j] + gauss[mid-j];
  }
  gauss[0] = gauss[2*mid-1] = 0;

  /* smooth */
  for (int i=0; i<n; i++ ) {
    FLOAT ggsum = 0;
    for (int j=-mid; j<mid; j++) {
      if (i+j >= 0 && i+j < n) {
	ggsum += gauss[mid+j];
	corr_sm[i] += corr[i+j] * gauss[mid+j];
      }
    }
    corr_sm[i] /= (ggsum) ? ggsum : 1.0;
  }
  
}

int convolve(FLOAT *corr, TS *ptsa, SHADOW *pfresPatt, int i_vel, int *Neff) {
  
  int i, j;
  int const N = ptsa->n;
  int const n = pfresPatt->nI[i_vel];
  const FLOAT width = SMOOTH_WIDTH_FSU;

  const FLOAT dt = ptsa->dt;
  FLOAT sigma = 0.5 * width * pfresPatt->fsu / pfresPatt->vRet[i_vel];
  const int Nsmth = (int) (sigma / dt);

  /*  smooth the time-series (subtract 1.0 to avoid edge effects) */
  FLOAT *corr_sm = det_malloc(N*sizeof(FLOAT));
  for(int i_x=0; i_x<N; ++i_x) {
    corr_sm[i_x] = 0;
  }
  for(i=0;i<N;i++) {
    corr[i] = ptsa->tI[i].I - 1.0;
  }
  smooth(corr, corr_sm, N, Nsmth);
  for(i=0;i<N;i++) {
    corr_sm[i] += 1.0;
  }

  /* do the cross-correlation */
  *Neff = 0;
  for (i=0; i<(N - n); i++) {
    
    FLOAT sumtmp = 0.0;
    for (j=0; j<n; j++) {

      sumtmp += (ptsa->tI[i+j].I - corr_sm[i+j]) * 
	(pfresPatt->I[i_vel][j] - 1.0);

      /* used to use 1.0 as a flagged value */
      if (ptsa->tI[i+j].I == 1.0) {
	sumtmp = 0.0;
	break;
      }
    }

    if ( sumtmp != 0) {
      (*Neff)++;
    }
    
    corr[i] = (ptsa->tI[i].I == 1.00000000) ? 0.0 : sumtmp;
    
  } 

  free(corr_sm);
  return(i);
}


STAT xcorr_norm(FLOAT *corr, FLOAT *corrNorm, TS *ptsa, SHADOW *pfresPatt, int nPts, FLOAT *meanAbsDev) {

  int i_vel = pfresPatt->i_vel;
  int i, j;
  FLOAT corrSum;
  FLOAT mean, rms;
  int nCorr;
  STAT stats;


  // get global stats to use for clipping in the normaliz'n routine
  stats = meanRmsClip_d(corr, nPts, SIGMA_CLIP, STRIDE);
  mean = stats.mean;
  rms = stats.rms;
  
  // initialize
  corrSum = 0;
  nCorr = 0;
  for (i=0; i<XCORR_RMS_WINDOW; i++) {
    if (  (corr[i] < SIGMA_CLIP*rms)   &&   (corr[i] != 0.0) ) {
      corrSum += fabs(corr[i]);
      nCorr++;
    }
  }
  corrNorm[0] = corr[0]/(corrSum/nCorr);
  
  // get the local standard deviation to normalize
  j=XCORR_RMS_WINDOW;
  *meanAbsDev = 0.0;
  
  const int n = pfresPatt->nI[i_vel];
  const int N = ptsa->n;
  for (i=1; i<(N - n); i++) {
      
    // if we're in the middle of the array
    if ( (i+j < (N - n) )  &&  (i-j >= 0)  ) {
      
      // add on the next one
      if (  (corr[i+j] < SIGMA_CLIP*rms)   &&   
	    (corr[i+j] != 0.0) ) {
	corrSum += fabs(corr[i+j]);
	nCorr++;
      }
      // subtract off the oldest
      if (  (corr[i-j] < SIGMA_CLIP*rms)   &&   
	    (corr[i-j] != 0.0) ) {
	corrSum -= fabs(corr[i-j]);
	nCorr--;
      }
      
      // if we're near the end of the array
    } else if (i+j > N) {
      
      // can't add on the next one
      // subtract off the oldest
      if (  (corr[i-j] < SIGMA_CLIP*rms)   &&   
	    (corr[i-j] != 0.0) ) {
	corrSum -= fabs(corr[i-j]);
	nCorr--;
      }
      
      
      // if we're at the beginning of the array
    } else if (i-j < 0 ) {
      
      // add on the next one
      if (  (corr[i+j] < SIGMA_CLIP*rms)   &&   
	    (corr[i+j] != 0.0) ) {
	corrSum += fabs(corr[i+j]);
	nCorr++;
      }
      // can't subtract off the oldest
      
    }
    
    // normalize the corr values by the local rms
    if (nCorr > 0) {
      corrNorm[i] = (corr[i]!=0.0) ?  (corr[i] / (corrSum/nCorr)) : 0.0;
    }
    
    *meanAbsDev += fabs(corr[i]);
    
  }
  *meanAbsDev /= (N - n);
  
  return stats;
}

int checkAutocorr(TS *pts, SHADOW *pfresPatt, FLOAT autocorr, FLOAT corrThresh, FLOAT *rms, int n_check) {

  const int N = pts->n;
  const int n = pfresPatt->nI[ pfresPatt->i_vel ];
  FLOAT xcorrs[n_check];

  /* just do it in the time domain */
  int i_last = 0;
  for (int i=0;i<n_check;i++) {
    i_last = i;
    int j0 = i * n;
    xcorrs[i] = 0;
    if (j0 + n > N) {
      break;
    }
    for (int j=0; j<n; j++) {
      xcorrs[i] +=  (pts->tI[j0+j].I - 1) * 
	(pfresPatt->I[pfresPatt->i_vel][j] - 1);
    }
  }
  STAT stats = meanRmsClip_d(xcorrs, i_last, SIGMA_CLIP, 1);

  /* normalize */
  FLOAT inv_rms = 1.0/stats.rms;
  for (int k=0; k < i_last; k++) {
    xcorrs[k] = (xcorrs[k] - stats.mean) * inv_rms;
  }
  autocorr = (autocorr - stats.mean) * inv_rms;
  
  int pass = 0;
  //printf ("ac: %g\n", autocorr);
  if ( autocorr > XCORR_COMPUTE_THRESH*corrThresh ) {
    pass = 1;
  }
  *rms = stats.rms;

  return(pass);
}


int xcorrelate (HIT hits[MAX_HITS], TS *ptsa, SHADOW *pfresPatt, int i_vel, FLOAT corrThresh, int dump, int norm) {

  int i, j, nhit, nPts;
  int index[MAX_HITS], nAbove[MAX_HITS], nMax=MAX_HITS;
  FLOAT *corr, autocorr, min, max;
  FLOAT *corrNorm, meanAbsDev;
  FLOAT mean, rms, above[MAX_HITS];
  int Ngoodpoints = 0;
  STAT stats;
  int const N = ptsa->n;
  int const n = pfresPatt->nI[i_vel];
  int const N2 = ptsa->n2;

  // autocorrelate to get the expected response for this kernel
  const FLOAT width = SMOOTH_WIDTH_FSU + pfresPatt->aa/pfresPatt->fsu;

  const FLOAT dt = ptsa->dt;
  FLOAT sigma = 0.5 * width * pfresPatt->fsu / pfresPatt->vRet[i_vel];
  const int Nsmth = (int) (sigma / dt);

  FLOAT *kern = det_calloc((n+12*Nsmth),sizeof(FLOAT));
  FLOAT *kern_sm = det_calloc((n+12*Nsmth),sizeof(FLOAT));
  for (j=0; j<n; j++) {
    kern[6*Nsmth+j] = pfresPatt->I[i_vel][j] - 1.0;
  }
  smooth(kern, kern_sm, n+12*Nsmth, Nsmth);
  autocorr=0;
  FLOAT autocorr_pure=0;
  for (j=0; j<n; j++) {
    autocorr_pure += SQR(pfresPatt->I[i_vel][j] - 1.0);
    autocorr += SQR(pfresPatt->I[i_vel][j] - (kern_sm[6*Nsmth+j]+1.0) );
  }
  free(kern);
  free(kern_sm);
  //printf ("autocorr: %g  %d\n", autocorr, pfresPatt->is_all_zero);
  autocorr_pure /=XCORR_CORRECTION;
  autocorr /= XCORR_CORRECTION;
  
  
  /* See if this kernel's autocorr will be above corrThresh */
  /* --need to correct corrThresh as checkAutocorr() doesn't smooth */
  FLOAT rms_tmp;
  if ( checkAutocorr(ptsa, pfresPatt, autocorr, corrThresh, 
		     &rms_tmp, N_CHECKAUTOCORR) == 0 ) {
    // printf ("skipping\n");
    hits[0].cormag0 = autocorr/rms_tmp;
    //fprintf (stderr, "Template r=%.1fm  d=%.1fAU ip=%.1fm has low signal, skipping.\n",
    //     pfresPatt->aa, pfresPatt->AU, pfresPatt->offset);
    return(0);
  }

  corr = det_malloc(N2 * sizeof(FLOAT));
  
  // calculate the convolution (cross correlation)
  if ( n < 10*log2(N) ) {
    //fprintf(stderr,"time domain\n");
    nPts = convolve(corr, ptsa, pfresPatt, i_vel, &Ngoodpoints);
    pfresPatt->used_fft = 0;
  } else {
    //fprintf(stderr,"frequency domain\n");
      
#if defined( USE_FFTW )
      nPts = fftwconvolve(corr, ptsa, pfresPatt, i_vel, &Ngoodpoints);
#else
      nPts = fftconvolve(corr, ptsa, pfresPatt, i_vel, &Ngoodpoints);
#endif

    pfresPatt->used_fft = 1;
  }
  pfresPatt->ts_eff_length = Ngoodpoints;

  // get global stats to use for clipping in the normaliz'n routine
  stats = meanRmsClip_d(corr, nPts, SIGMA_CLIP, STRIDE);
  mean = stats.mean;
  rms = stats.rms;
  meanAbsDev = rms;
  //printf ("corrrms: %g\n", rms);

    // normalize autocorr
  autocorr /= meanAbsDev;


  // // // NORMALIZATION // // // 
  if (norm) {

    /* make sure the wavetable and workspace are there */
    if ( ! ptsa->FFTisdone ) {
#if defined( USE_FFTW )
	FFTWtimeseries(ptsa);
#else
	FFTtimeseries(ptsa);
#endif
    }

    corrNorm = det_malloc(ptsa->n2 * sizeof(FLOAT));
#if defined( USE_FFTW )
    fftw_xcorr_norm(corr, corrNorm, ptsa, pfresPatt, rms);
#else
    fft_xcorr_norm(corr, corrNorm, ptsa, pfresPatt, rms);
#endif

    
  } else {

    corrNorm = corr;    //long const ni = ptsa->n;
    FLOAT inv_rms = 1.0/rms;
    for(i=0; i<N; i++) {
      corrNorm[i] = (corrNorm[i]==0.0) ? 
	0.0 : (corrNorm[i] - mean) * inv_rms;
    }

  }
  //fprintf(stderr, "mean: %f rms: %f\n", mean, rms);
  stats.rms = 1.0; // by definition
  stats.mean = 0.0; // by definition
  mean = stats.mean;
  rms = stats.rms;
  
  // get min and max acceptable value (but not below corrThresh
  min = ( ((1.0-XCORR_RMS_PERC)*autocorr-XCORR_RMS_RANGE*rms) > (mean+corrThresh*rms)) ? 
    ( (1.0-XCORR_RMS_PERC)*autocorr-XCORR_RMS_RANGE*rms ) : ( mean+corrThresh*rms );
  max = (1.0+XCORR_RMS_PERC)*autocorr + XCORR_RMS_RANGE*rms;
  
  nhit = window(corrNorm, nPts, min, max, above, index, nMax);

  const FLOAT hit_width = 
    (FSU_PROX_SAME * pfresPatt->fsu / pfresPatt->vRet[i_vel]) / ptsa->dt;
  
  nhit = reducelist(above, index, nAbove, nhit, GET_MAX, hit_width); 
  //fprintf("nPts:%d  min:%.3g  max:%.3g  auto:%.3g  n:%d\n", 
  // nPts, min, max, autocorr, n);
  
  FLOAT inv_rms2 = 1.0 / rms;
  hits[0].cormag0  = (autocorr - mean) * inv_rms2; // if nhit==0
  for (i=0; i<nhit; i++) {
    //fprintf (stderr, "%.2f %.2f %.2f  %.2f\n", min, autocorr, max, above[i]);  
    //exit(0);
    hits[i].cormag   = (above[i] - mean) * inv_rms2;
    // add on the cent pt (kern coord is the left-most pt, not mid)
    hits[i].index    = index[i] + pfresPatt->icentre[i_vel];
    hits[i].n        = nAbove[i];
    hits[i].cormag0  = (autocorr - mean) * inv_rms2;

    //printf ("cormag: %.2f   i: %d   n: %d  cormag0: %.2f\n", hits[i].cormag, hits[i].index, hits[i].n, hits[i].cormag0 );
  }

  // dump to a file

  if (dump) {
    writeXcorrSeries(ptsa, pfresPatt, corr, 
		     corrNorm, hits, nhit, "xcor");
  }

  free(corr);
  if (norm) {
    free(corrNorm);
  }
  return nhit;
}


/* KO = known objects */
int xcorrelateKO (HIT hits[MAX_HITS], TS *ptsa, SHADOW *pfresPatt, int i_vel, FLOAT corrThresh, int dump, int norm) {

  int i, j, nPts;
  FLOAT *corr, autocorr;
  FLOAT *corrNorm, meanAbsDev;
  FLOAT mean, rms;
  int Ngoodpoints = 0;
  STAT stats;
  int const N = ptsa->n;
  int const n = pfresPatt->nI[i_vel];
  int const N2 = ptsa->n2;
  int const nhit = ptsa->nKO;


  // autocorrelate to get the expected response for this kernel
  const FLOAT width = SMOOTH_WIDTH_FSU;
  const FLOAT dt = ptsa->dt;
  FLOAT sigma = 0.5 * width * pfresPatt->fsu / pfresPatt->vRet[i_vel];
  const int Nsmth = (int) (sigma / dt);

  FLOAT *kern = det_calloc((n+12*Nsmth),sizeof(FLOAT));
  FLOAT *kern_sm = det_calloc((n+12*Nsmth),sizeof(FLOAT));
  for (j=0; j<n; j++) {
    kern[6*Nsmth+j] = pfresPatt->I[i_vel][j] - 1.0;
  }
  smooth(kern, kern_sm, n+12*Nsmth, Nsmth);
  autocorr=0;
  FLOAT autocorr_pure=0;
  for (j=0; j<n; j++) {
    autocorr_pure += SQR(pfresPatt->I[i_vel][j] - 1.0);
    autocorr += SQR(pfresPatt->I[i_vel][j] - (kern_sm[6*Nsmth+j]+1.0) );
  }
  free(kern);
  free(kern_sm);
  //printf ("autocorr: %g  %d\n", autocorr, pfresPatt->is_all_zero);
  autocorr_pure /=XCORR_CORRECTION;
  autocorr /= XCORR_CORRECTION;
  
  
  /* See if this kernel's autocorr will be above corrThresh */
  /* --need to correct corrThresh as checkAutocorr() doesn't smooth */
  FLOAT rms_tmp;
  if ( checkAutocorr(ptsa, pfresPatt, autocorr, corrThresh, 
		     &rms_tmp, N_CHECKAUTOCORR) == 0 ) {
    // printf ("skipping\n");
    hits[0].cormag0 = autocorr/rms_tmp;
    return(0);
  }

  // autocorrelate to get the expected response for this kernel
  //autocorr=0;
  //for (j=0; j<n; j++) {
  //  autocorr += SQR(pfresPatt->I[i_vel][j] - 1.0);
  // }
  //printf ("autocorr: %g  %d\n", autocorr, pfresPatt->is_all_zero);
  //autocorr /= XCORR_CORRECTION;
  

  /* See if this kernel's autocorr will be above corrThresh */
  //FLOAT rms_tmp;
  // if ( checkAutocorr(ptsa, pfresPatt, autocorr, corrThresh, &rms_tmp) == 0 ) {
  //  // printf ("skipping\n");
    
  //  return(0);
  /// }

  corr = det_malloc(N2 * sizeof(FLOAT));
  
  // calculate the convolution (cross correlation)
  if ( n < 10*log2(N) ) {
    //fprintf(stderr,"time domain\n");
    nPts = convolve(corr, ptsa, pfresPatt, i_vel, &Ngoodpoints);
    pfresPatt->used_fft = 0;
  } else {
    //fprintf(stderr,"frequency domain\n");
#if defined( USE_FFTW )
      nPts = fftwconvolve(corr, ptsa, pfresPatt, i_vel, &Ngoodpoints);
#else
      nPts = fftconvolve(corr, ptsa, pfresPatt, i_vel, &Ngoodpoints);
#endif

    pfresPatt->used_fft = 1;
  }
  pfresPatt->ts_eff_length = Ngoodpoints;

  // get global stats to use for clipping in the normaliz'n routine
  stats = meanRmsClip_d(corr, nPts, SIGMA_CLIP, STRIDE);
  mean = stats.mean;
  rms = stats.rms;
  meanAbsDev = rms;
  //printf ("corrrms: %g\n", rms);

    // normalize autocorr
  autocorr /= meanAbsDev;


  // // // NORMALIZATION // // // 
  if (norm) {

    /* make sure the wavetable and workspace are there */
    if ( ! ptsa->FFTisdone ) {
#if defined( USE_FFTW )
	FFTWtimeseries(ptsa);
#else
	FFTtimeseries(ptsa);
#endif
    }

    corrNorm = det_malloc(ptsa->n2 * sizeof(FLOAT));
#if defined( USE_FFTW )
    fftw_xcorr_norm(corr, corrNorm, ptsa, pfresPatt, rms);
#else
    fft_xcorr_norm(corr, corrNorm, ptsa, pfresPatt, rms);
#endif

    
  } else {

    corrNorm = corr;    //long const ni = ptsa->n;
    FLOAT inv_rms = 1.0/rms;
    for(i=0; i<N; i++) {
      corrNorm[i] = (corrNorm[i]==0.0) ? 
	0.0 : (corrNorm[i] - mean) * inv_rms;
      //printf("%d %g %g\n", i, corr[i], corrNorm[i]);
    }
    //exit(0);
  }

  stats.rms = 1.0; // by definition
  stats.mean = 0.0; // by definition
  mean = stats.mean;
  rms = stats.rms;
  
  FLOAT inv_rms2 = 1.0/rms;
  int i_tmp;
  for (i=0; i<nhit; i++) {
    i_tmp = ptsa->KOlist[i] - pfresPatt->icentre[i_vel];    
    hits[i].cormag   = (corrNorm[i_tmp] - mean) * inv_rms2;
    hits[i].index    = i_tmp + pfresPatt->icentre[i_vel];
    hits[i].n        = 1;
    hits[i].cormag0  = (autocorr - mean) * inv_rms2;
  }
  
  // dump to a file
  if (dump) {
    writeXcorrSeries(ptsa, pfresPatt, corr, 
		     corrNorm, hits, nhit, "xcor");
  }

  free(corr);
  if (norm) {
    free(corrNorm);
  }
  return nhit;
}


int xchi (HIT *hits, TS *ptsa, SHADOW *pfresPatt, int i_vel, FLOAT chiThresh, int dump) {

  int i,j, nhit, nPts;
  int index[MAX_HITS], nAbove[MAX_HITS], nMax=MAX_HITS;
  FLOAT *chi, autochi;
  FLOAT diff, sumdiffsqr, mean, rms, above[MAX_HITS];
  STAT stats;
  const int N = ptsa->n;
  const int n = pfresPatt->nI[i_vel];

  chi = det_malloc(N*sizeof(FLOAT));
  for (i=0; i<( N - n ); i++) {
    sumdiffsqr = 0;
    for (j=0; j<n; j++) {
      diff = (ptsa->tI[i+j].I - pfresPatt->I[i_vel][j]);
      sumdiffsqr += diff*diff;
    }
    // chi[i] = sumdiffsqr;
    chi[i] = (ptsa->tI[i].I == 1.0) ? sumdiffsqr : 0.0 ;
  } 
  nPts = i;
  
  // flag anything above the threshhold
  //mean    = meanclip_d(chi, nPts, SUBSAMPLES, SIGMA_CLIP, 0, 2);
  //rms     = rmsclip_d(chi, nPts, SUBSAMPLES, SIGMA_CLIP, 0, 3);
  stats = meanRmsClip_d(chi, nPts, SIGMA_CLIP, STRIDE);
  mean = stats.mean;
  rms = stats.rms;
  
  autochi = n * (ptsa->rms * ptsa->rms);

  nhit = window(chi, nPts, 0, (mean+chiThresh*rms), 
		above, index, nMax);

  const FLOAT hit_width = 
    (FSU_PROX_SAME * pfresPatt->fsu / pfresPatt->vRet[i_vel]) / ptsa->dt;

  nhit = reducelist(above, index, nAbove, nhit, GET_MIN, hit_width);
 
  const FLOAT inv_rms = 1.0/rms;
  for (i=0; i<n; i++) {
    hits[i].chimag   = (above[i] - mean) * inv_rms;
    // add on the cent pt (kern coord is the left-most pt, not mid)
    hits[i].index = index[i] + pfresPatt->icentre[i_vel];
    hits[i].n     = nAbove[i];
    hits[i].chimag0  = ((autochi - mean) * inv_rms)/XCHI_CORRECTION;
  }


  // dump to a file
  if (dump) {
    
    FILE *fp;
    char outfile[MAX_FILENAME];

    snprintf (outfile, MAX_FILENAME, "%s.xchi",ptsa->tsfile);
    fp = openfile(outfile, "w");
    
    // print a header
    fprintf (fp, "# Created by xchi\n");
    fprintf (fp, "# Fresnel_Pattern_File: %s\n", pfresPatt->filename);
    fprintf (fp, "# %s \tSolar_elongation: \t%.1f\n", LELONG, 180.0*ptsa->elong/PI);
    fprintf (fp, "# %s \tKBO_radius: \t%g\n", LAA, pfresPatt->aa);
    fprintf (fp, "# %s \tdist_to_KBO_(AU): \t%g\n", LAU, pfresPatt->AU);
    fprintf (fp, "# %s \toffset_from_line_of_sight: \t%g\n", LOFFSET, pfresPatt->offset);
    fprintf (fp, "# %s \tlambda_Low: \t%g\n", LLAMBLO, pfresPatt->lambLo);
    fprintf (fp, "# %s \tlambda_High: \t%g\n", LLAMBHI, pfresPatt->lambHi);
    fprintf (fp, "# %s \tNum_lambdas: \t%d\n", LNLAMBDA, pfresPatt->nLambda);
    fprintf (fp, "# %s \tProjected_Star_Radius: \t%g\n", LRSTAR, pfresPatt->RStar);
    fprintf (fp, "# %s \tNum_point_sources_for_star: \t%d\n", LNSTARS, pfresPatt->nStars);


    for (i=0; i<nhit; i++) {
      fprintf (fp, "# %d %.4f %.3f %.3f %d\n", hits[i].index, hits[i].index*ptsa->dt, hits[i].chimag, hits[i].chimag0, hits[i].n);
    }

    for (i=0; i<(N - n); i++) {
      fprintf (fp, "%d %.4f %.6f\n", i, i*ptsa->dt, chi[i]);
    }

    closefile(fp);

  }
  free(chi);

  return nhit;
}

int xchi2 (HIT *hits, TS *ptsa, SHADOW *pfresPatt, int i_vel, SERIES *chi) {

  int i, j, nPts;
  FLOAT diff, sumdiffsqr;
  STAT stats;
  const int N = ptsa->n;
  const int n = pfresPatt->nI[i_vel];
  
  for (i=0; i<( N - n ); i++) {
    sumdiffsqr = 0;
    
    for (j=0; j<n; j++) {
      diff = (ptsa->tI[i+j].I - pfresPatt->I[i_vel][j]);
      sumdiffsqr += diff*diff;
    }
    chi->value[i+pfresPatt->icentre[i_vel]] = sumdiffsqr;

  } 
  nPts = i;

  // flag anything above the threshhold
  //chi->mean = meanclip_d(chi->value, nPts, SUBSAMPLES, SIGMA_CLIP, 0, 2);
  //chi->rms  = rmsclip_d(chi->value, nPts, SUBSAMPLES, SIGMA_CLIP, 0, 3);
  
  stats = meanRmsClip_d(chi->value, nPts, SIGMA_CLIP, STRIDE);
  chi->mean = stats.mean;
  chi->rms  = stats.rms;
  chi->autoval = n * (ptsa->rms * ptsa->rms);
  chi->n       = nPts + pfresPatt->icentre[i_vel];

  hits->chimag0 = ( (chi->autoval - chi->mean) / chi->rms )/XCHI_CORRECTION;

  return chi->n;
}


/* same as xchi2 but only do the values for which there are xcorr hits */
int xchi3 (HIT *hits, TS *ptsa, SHADOW *pfresPatt, int i_vel, SERIES *chi, HIT corrhits[], int ncorr) {

  int i, j, k;
  FLOAT diff, sumdiffsqr, sumIsqr, sumI;
  STAT stats;
  int dof;
  int Nrandom = 200;
  const int N = ptsa->n;
  const int n = pfresPatt->nI[i_vel];
  FLOAT r, baryc, m;

  SERIES chi_rand;
  chi_rand.value = det_malloc(Nrandom*sizeof(FLOAT));
  for (k=0; k<ncorr; k++) {

    i = corrhits[k].index - pfresPatt->icentre[i_vel];

    sumdiffsqr = 0.0;
    sumIsqr = 0.0;
    sumI = 0.0;
    baryc = 0.0;
    for (j=0; j<n; j++) {
      diff = (ptsa->tI[i+j].I - pfresPatt->I[i_vel][j]);
      sumIsqr += SQR( ptsa->tI[i+j].I - 1.0 );
      sumdiffsqr += diff*diff;
      r = (i + j - corrhits[k].index) * pfresPatt->vRet[i_vel] * ptsa->dt /
	pfresPatt->fsu;
      m = (1.0-ptsa->tI[i+j].I) * fabs(pfresPatt->I[i_vel][j]);
      baryc += r * m;
      sumI += m;
      //printf("%d %d %.6f %.6f %.6f %.6f\n",i, j, diff, r, fabs(r), baryc);
    }
    dof = n - 2;
    chi->value[i+pfresPatt->icentre[i_vel]] = sumdiffsqr;
    corrhits[k].red_chi2null = sumIsqr / ( dof * SQR(ptsa->rms ) );
    corrhits[k].red_chi2 = sumdiffsqr / ( dof * SQR(ptsa->rms) );
    corrhits[k].Pred_chi2 = Predchi2(corrhits[k].red_chi2,dof);
    corrhits[k].dof = dof;

    corrhits[k].baryc = baryc / sumI;
    corrhits[k].used_fft = pfresPatt->used_fft;
  } 
  

  /* check a few random locations to get the stats */
  for (k=0;k<Nrandom;k++) {
    i = (int) ((k+0.5) * (N - n) / Nrandom);
    sumdiffsqr = 0.0;
    for (j=0; j<n; j++) {
      diff = (ptsa->tI[i+j].I - pfresPatt->I[i_vel][j]);
      sumdiffsqr += diff*diff;
    }
    chi_rand.value[k] = sumdiffsqr;
  }

  // flag anything above the threshhold
  //chi->mean = meanclip_d(chi->value, nPts, SUBSAMPLES, SIGMA_CLIP, 0, 2);
  //chi->rms  = rmsclip_d(chi->value, nPts, SUBSAMPLES, SIGMA_CLIP, 0, 3);
  stats = meanRmsClip_d(chi_rand.value, Nrandom, SIGMA_CLIP, 1);
  chi->mean = stats.mean;
  chi->rms  = stats.rms;
  chi->autoval = n * (ptsa->rms * ptsa->rms);
  chi->n       = ncorr;
  
  hits->chimag0 = ( (chi->autoval - chi->mean) / chi->rms )/XCHI_CORRECTION;

  free(chi_rand.value);
  return chi->n;
}


int reducelist(FLOAT *values, int *index, int *nXtrem, int n, int sign, int hit_width) {

  int i, j, xind;
  FLOAT xval;

  
  // step through the hits
  xval = values[0];
  xind = index[0];
  nXtrem[0] = 1;
  j = 0;
  for (i=0; i<n; i++) {
    
    // if this is still the same hit
    if ( fabs(index[i] - xind) < hit_width ) {
      
      nXtrem[j]++;        // count these

      // if 'sign' is set ... get the max, otherwise get the min
      if (sign) {
	if ( values[i] > xval ) { 
	  xind = index[i];  
	  xval = values[i];
	}
      } else {
	if ( values[i] < xval ) { 
	  xind = index[i];
	  xval = values[i];
	}
      }

    // if this is a new hit
    } else {
    
      // save the old values
      values[j] = xval;
      index[j++]  = xind;
      
      nXtrem[j] = 1;        // initialize the new counter
      xind = index[i];
      xval = values[i];
    }

  }

  // need to record the stuff for the last detection

  if (n) {
    values[j] = xval;
    index[j++] = xind;
  }

  return j;

}



int merge(HIT *corlist, HIT *chilist, int n1, int n2) {

  int i,j,nMatch;
  HIT corlisttmp[MAX_HITS], chilisttmp[MAX_HITS];
  
  // check each point from xcorr hits
  nMatch=0;
  for (i=0; i<n1; i++) {
    for (j=0; j<n2; j++) {      
      if ( fabs(corlist[i].index - chilist[j].index)<PROXIMITY ) {
	corlisttmp[nMatch].index = corlist[i].index;
	chilisttmp[nMatch].index = chilist[j].index;
	corlisttmp[nMatch].cormag = corlist[i].cormag;
	chilisttmp[nMatch].chimag = chilist[j].chimag;
	corlisttmp[nMatch].n = corlist[i].n;
	chilisttmp[nMatch++].n = chilist[j].n;
      }
    }
  }

  // now give back the same lists with the matched numbers in them
  for (i=0; i<nMatch; i++) {
    corlist[i].index = corlisttmp[i].index;
    chilist[i].index = chilisttmp[i].index;
    corlist[i].cormag = corlisttmp[i].cormag;
    chilist[i].chimag = chilisttmp[i].chimag;
    corlist[i].n = corlisttmp[i].n;
    chilist[i].n = chilisttmp[i].n;
  }

  return nMatch;

}

int merge2(HIT *corlist, HIT *chilist, int nCorHit, int nChiHit, SERIES *chi, PARAMLIST *pars) {

  int i,nMatch;
  FLOAT mean     = chi->mean;
  FLOAT rms      = chi->rms;
  FLOAT autochi  = chi->autoval;

  FLOAT chi_rms;
  FLOAT cormag0;

  // check each point from xcorr hits
  nMatch=0;
  for (i=0; i<nCorHit; i++) {

    cormag0 = corlist[i].cormag0;
    chi_rms = (chi->value[corlist[i].index] - mean)/rms;

    // if it's also high in xchi, it must be noise
    if (  (chi_rms < pars->chiThresh)  &&  
	  ( cormag0 > pars->corrThresh ) ) {

      corlist[nMatch].baryc    = corlist[i].baryc;
      corlist[nMatch].red_chi2null = corlist[i].red_chi2null;
      corlist[nMatch].red_chi2 = corlist[i].red_chi2;
      corlist[nMatch].Pred_chi2 = corlist[i].Pred_chi2;
      corlist[nMatch].dof      = corlist[i].dof;
      corlist[nMatch].index    = corlist[i].index;
      corlist[nMatch].cormag   = corlist[i].cormag;
      corlist[nMatch].chimag   = ( chi->value[corlist[i].index]-mean )/rms;
      corlist[nMatch].chimag0  = ( autochi - mean )/rms/XCHI_CORRECTION;
      corlist[nMatch].n        = corlist[i].n;
      corlist[nMatch].used_fft = corlist[i].used_fft;


      chilist[nMatch].baryc    = corlist[i].baryc;
      chilist[nMatch].red_chi2null = corlist[i].red_chi2null;
      chilist[nMatch].red_chi2 = corlist[i].red_chi2;
      chilist[nMatch].Pred_chi2 = corlist[i].Pred_chi2;
      chilist[nMatch].dof      = corlist[i].dof;
      chilist[nMatch].index    = corlist[i].index;
      chilist[nMatch].cormag   = corlist[i].cormag;
      chilist[nMatch].chimag   = corlist[nMatch].chimag;
      chilist[nMatch].chimag0  = corlist[nMatch].chimag0;
      chilist[nMatch].n      = corlist[i].n;
      chilist[nMatch].used_fft = corlist[nMatch].used_fft;
      nMatch++;
    }
  }

  return nMatch;

}


int countHits(HIT *hits, int nhit) {

  int i, j, ncullhit=0;
  int seen;

  for (i=0;i<nhit;i++) {
    
    /* check and see if we've seen this one before */
    seen = 0;
    for (j=0;j<ncullhit;j++) {
      if ( fabs(hits[i].index - hits[j].index) < SAME_HIT_RANGE ) {
	seen++;
      }
    }
    
    /* if seen, ignore it -- if new count it */
    if ( ! seen ) {
      ncullhit++;
    }

  }

  return ncullhit;
}

int addVerify(int *addlist, int nAdd, HIT *corhits, HIT *chihits, int *nHit, HIT *recovered, SHADOW *pfresPatt, TS *pts) {

  int i,j,n;
  int recover;
  const int i_vel = pfresPatt->i_vel;
  int nRecover = pfresPatt->nRecover[i_vel];
  const FLOAT fsu_index = (FSU_PROX_SAME*pfresPatt->fsu/pfresPatt->vRet[i_vel])/pts->dt;
  const FLOAT prox = (fsu_index > PROXIMITY) ? (fsu_index) : PROXIMITY;
  

  n=0;
  for (i=0; i<(*nHit); i++) {

    recover = 0;
    for (j=0; j<nAdd; j++) {

      //if (pfresPatt->offset > 0) {
      //fprintf (stderr, "corind: %d  chiind: %d  addlist: %d\n", corhits[i].index, chihits[i].index, addlist[j]);
      //}

      // if we added it, log the detection
      if ( (fabs(corhits[i].index - addlist[j]) < prox) ) {
	//fprintf (stderr, "%.1f %.1f  %d %d\n", prox, fsu_index, corhits[i].index, addlist[j]);
	// &&  (fabs(chihits[i].index - addlist[j]) < PROXIMITY)    )  {

	recovered[nRecover].cormag = corhits[i].cormag;
	recovered[nRecover].chimag = chihits[i].chimag;
	recovered[nRecover].index = chihits[i].index;
	nRecover++;
	recover++;

      }
    }
    // if we didn't add it ... it's a hit candidate
    if (! recover) {
      corhits[n].index    = corhits[i].index;
      corhits[n].cormag   = corhits[i].cormag;
      corhits[n].cormag0  = corhits[i].cormag0;
      corhits[n].n        = corhits[i].n;
      
      chihits[n].index    = chihits[i].index;
      chihits[n].chimag   = chihits[i].chimag;
      chihits[n].chimag0  = chihits[i].chimag0;
      chihits[n].n        = chihits[i].n;
      n++;
    }
    
  }  // end for i<(*nHit)

  //if (pfresPatt->offset > 0) {
    //fprintf(stderr, "%d\n", n);
  //  exit(0);
  //}

  (*nHit) = n;
  pfresPatt->nHit[i_vel] = n;
  return nRecover;
}

void printHitHeader(FILE *fp) {
  fprintf (fp, "%-6s %6s %6s %6s   %7s %6s %6s %4s   %7s %7s %7s %4s  %6s %5s %6s %3s  %5s  %3s\n", 
	   "# RStar", "aa", "AU", "offset", 
	   "cor_i", "cormg", "cormg0", "co_n", 
	   "chi_i", "chimg", "chimg0", "ch_n", 
	   "rX2", "P>X2", "rXnul", "dof", "baryc", "fft" );
}

void printHit(FILE *fp, SHADOW *pfresPatt, HIT corrhits, HIT chihits) {
  fprintf (fp, "%-6.0f %6.0f %6.0f %6.0f   %7d %6.3f %6.3f %4d   %7d %7.3f %7.3f %4d  %6.3f %5.3e %6.3f %3d  %5.2f %3d\n", 
	   pfresPatt->RStar, pfresPatt->aa, 
	   pfresPatt->AU, pfresPatt->offset, 
	   corrhits.index, corrhits.cormag, corrhits.cormag0, 
	   corrhits.n, chihits.index, chihits.chimag, 
	   chihits.chimag0, chihits.n, 
	   chihits.red_chi2, 1.0-chihits.Pred_chi2, 
	   chihits.red_chi2null, chihits.dof, chihits.baryc, chihits.used_fft
	   );  
}


void writeXcorrSeries(TS *ptsa, SHADOW *pfresPatt, FLOAT *corr, FLOAT *corrNorm, HIT hits[MAX_HITS], int n, char *suffix) {

  int i_vel = pfresPatt->i_vel;
  int i;
  FILE *fp;
  char outfile[MAX_FILENAME];
  int N = ptsa->n - pfresPatt->nI[i_vel];

  if (ptsa->fits) {

    float t_total = (float) ptsa->n*ptsa->dt;
    float *fitsdata, dt = ptsa->dt;
    fitsfile *fptr;
    int status=0, tfields=1;
    char extname[] = "Flux_norm_add";
    char *ttype[] = {"I"};
    char *tform[] = {"1E"};
    char *tunit[] = {"NA"};
    LONGLONG firstrow=1, firstelem=1;

    fitsdata = det_malloc(N*sizeof(float));
    for (i=0;i<N;i++) {
      fitsdata[i] = (float) corrNorm[i];
    }

    char path[MAX_FILENAME];
    snprintf(path, MAX_FILENAME, "%s.%s.fits", ptsa->tsfile, suffix);
    
    if (! fits_create_file(&fptr, path, &status) ) {
      
      /* write dt to the header */
      fits_create_tbl(fptr, BINARY_TBL, n, tfields, 
                      ttype, tform, tunit, extname, &status);
      
      fits_write_key(fptr, TFLOAT, "SAMPTIME", &dt, 
                     "Sampling time", &status);
      fits_write_key(fptr, TFLOAT, "DURATION", &t_total, 
                     "Duration of timeseries", &status);
   
      // print a header
      FLOAT elongdeg = 180.0 * ptsa->elong/PI;
      fits_write_key(fptr, TSTRING, "FRESPATT", pfresPatt->filename,
		     "Fresnel Pattern Filename", &status);
      fits_write_key(fptr, TDOUBLE, "ELONG", &elongdeg, 
		     "Solar Elongation", &status);
      fits_write_key(fptr, TDOUBLE, "KBO_RAD", &pfresPatt->aa, 
		     "KBO Radius", &status);
      fits_write_key(fptr, TDOUBLE, "KBO_DIST", &pfresPatt->AU, 
		     "Distance to KBO (AU)", &status);
      fits_write_key(fptr, TDOUBLE, "IMPAR",  &pfresPatt->offset, 
		     "Occultation Impact Parameter", &status);
      fits_write_key(fptr, TDOUBLE, "LAMBDALO",  &pfresPatt->lambLo,
		     "Lowest wavelength (m)", &status);
      fits_write_key(fptr, TDOUBLE, "LAMBDAHI",  &pfresPatt->lambHi,
		     "Highest wavelength (m)", &status);
      fits_write_key(fptr, TINT, "NLAMBDA", &pfresPatt->nLambda, 
		     "Number of lambdas", &status);
      fits_write_key(fptr, TDOUBLE, "RSTAR",  &pfresPatt->RStar, 
		     "Projected Stellar Radius", &status);
      fits_write_key(fptr, TINT, "NSTAR", &pfresPatt->nStars, 
		     "Number of sources in stellar disk", &status);
      

      char label[9];
      char value[21];
      for (i=0; i<n; i++) {
          snprintf(label, 9, "HIT%03d",i);
          snprintf(value, 21, "%7d  %4.1f %4.1f", 
                   hits[i].index, hits[i].cormag, hits[i].cormag0);
          fits_write_key(fptr, TSTRING, label, value, "i cormag cormag0", 
                         &status);
      }
      
      /* write the data */
      fits_write_col(fptr, TFLOAT, 1, firstrow, firstelem,
                     N, fitsdata, &status);
      
    }
    fits_close_file(fptr, &status);
    if (status) {
      fits_report_error(stderr, status);
    }


  } else {
      snprintf (outfile, MAX_FILENAME, "%s.xcor",ptsa->tsfile);
    fp = openfile(outfile, "w");
    
    // print a header
    fprintf (fp, "# Created by xcorrelate\n");
    fprintf (fp, "# Fresnel_Pattern_File: %s\n", pfresPatt->filename);
    fprintf (fp, "# 1 index\n");
    fprintf (fp, "# 2 time\n");
    fprintf (fp, "# 3 corrNorm\n");
    fprintf (fp, "# 4 corr\n");
    fprintf (fp, "# %s Solar_elongation: %.1f\n", LELONG, 180.0*ptsa->elong/PI);
    fprintf (fp, "# %s \tKBO_radius: \t%g\n", LAA, pfresPatt->aa);
    fprintf (fp, "# %s \tdist_to_KBO_(AU): \t%g\n", LAU, pfresPatt->AU);
    fprintf (fp, "# %s \toffset_from_line_of_sight: \t%g\n", LOFFSET, pfresPatt->offset);
    fprintf (fp, "# %s \tlambda_Low: \t%g\n", LLAMBLO, pfresPatt->lambLo);
    fprintf (fp, "# %s \tlambda_High: \t%g\n", LLAMBHI, pfresPatt->lambHi);
    fprintf (fp, "# %s \tNum_lambdas: \t%d\n", LNLAMBDA, pfresPatt->nLambda);
    fprintf (fp, "# %s \tProjected_Star_Radius: \t%g\n", LRSTAR, pfresPatt->RStar);
    fprintf (fp, "# %s \tNum_point_sources_for_star: \t%d\n", LNSTARS, pfresPatt->nStars);
    
    for (i=0; i<n; i++) {
      fprintf (fp, "# %d %.4f %.3f %.3f %d  %d %.4f\n", hits[i].index, hits[i].index*ptsa->dt, hits[i].cormag, hits[i].cormag0, hits[i].n, pfresPatt->icentre[i_vel], pfresPatt->t_centre);
    }
    
    for (i=0; i<N; i++) {
      fprintf (fp, "%d %.4f %.8f %.8f\n", i, i*ptsa->dt, corrNorm[i], corr[i]);
    }
    
    closefile(fp);
  }
}

void printLenHeader(FILE *fp) {
  fprintf (fp, "# %-6s %-6s   %6s %6s  %8s %8s  %5s\n","rad", "AU", "total", "eff", "t_total", "t_eff", "percent");
}

void printLen(FILE *fp, TS *ptsa, SHADOW *pfresPatt) {
  fprintf (fp, "%6.1f %6.1f   %6d  %6d   %8.3f %8.3f    %5.3f\n", 
	   pfresPatt->aa, pfresPatt->AU, ptsa->n,  pfresPatt->ts_eff_length,
	   ptsa->n*ptsa->dt, pfresPatt->ts_eff_length*ptsa->dt,
	   pfresPatt->ts_eff_length * ptsa->dt / (ptsa->n*ptsa->dt)    );
}

void printStatsHeader(FILE *fp) {
  fprintf (fp, "%-6s %6s %6s %6s   %4s %4s %4s   %6s %5s %5s %5s %6s %6s %5s\n", 
	   "# RStar", "aa", "AU", "offset",  
	   "nRec", "nAdd", "nHit",
	   "vRet", "mnCor", "cor0", "rmsCo",
	   "mnChi", "chi0", "rmsCh"        );
}

void printStats(FILE *fp, SHADOW *pfresPatt, HIT corrhits, HIT chihits, PARAMLIST *pars, int i_vel) {
      fprintf (fp, "%-6.0f %6.0f %6.0f %6.0f   %4d %4d %4d   %6.0f %5.2f %5.2f %5.2f %6.2f %6.2f %5.2f\n", 
	       pfresPatt->RStar,
	       pfresPatt->aa,
	       pfresPatt->AU,
	       pfresPatt->offset,

	       pfresPatt->nRecover[i_vel],
	       pars->NperCycle * pars->cycles,
	       pfresPatt->nHit[i_vel],

	       pfresPatt->vRet[i_vel],
	       pfresPatt->meanCor[i_vel],
	       corrhits.cormag0,
	       pfresPatt->rmsCor[i_vel],

	       pfresPatt->meanChi[i_vel],
	       chihits.chimag0,
	       pfresPatt->rmsChi[i_vel]        
	       );

}

void printBStatsHeader(FILE *fp) {
  fprintf (fp, "%-6s %6s %6s %6s %6s %6s %6s   %4s %4s %4s   %6s %5s %5s %5s %6s %6s %5s\n", 
	   "# RStar", "aa", "AU", "offset", "aaA", "AUA", "offsA",  
	   "nRec", "nAdd", "nHit",
	   "vRet", "mnCor", "cor0", "rmsCo",
	   "mnChi", "chi0", "rmsCh"        );
}

void printBStats(FILE *fp, SHADOW *pfresPatt, SHADOW *pfresPattA, HIT corrhits, HIT chihits, PARAMLIST *pars, int i_vel) {
      fprintf (fp, "%-6.0f %6.0f %6.0f %6.0f %6.0f %6.0f %6.0f   %4d %4d %4d   %6.0f %5.2f %5.2f %5.2f %6.2f %6.2f %5.2f\n", 
	       pfresPatt->RStar,
	       pfresPatt->aa,
	       pfresPatt->AU,
	       pfresPatt->offset,
	       pfresPattA->aa,
	       pfresPattA->AU,
	       pfresPattA->offset,

	       pfresPatt->nRecover[i_vel],
	       pars->NperCycle * pars->cycles,
	       pfresPatt->nHit[i_vel],

	       pfresPatt->vRet[i_vel],
	       pfresPatt->meanCor[i_vel],
	       corrhits.cormag0,
	       pfresPatt->rmsCor[i_vel],

	       pfresPatt->meanChi[i_vel],
	       chihits.chimag0,
	       pfresPatt->rmsChi[i_vel]        
	       );

}

void errQuit (char *msg, ...) {

  char s[MAX_LINE_LENGTH];

  va_list args;
  va_start(args,msg);
  vsprintf(s, msg, args);
  va_end(args);

  perror(s);
  exit(EXIT_FAILURE);
}



void debugOutput (SHADOW *pfresPatt[], PARAMLIST *pars, TS *pts, TS *ptsa, int Npfres) {

  // input parameters
  fprintf (stderr,"cycles: %d  NperCycle:  %d   corrThresh: %.3f    chiThresh: %.3f   fresDir: %s   Noffset: %d    norm: %d\n", pars->cycles, pars->NperCycle, pars->corrThresh,  pars->chiThresh,  pars->fresDir, pars->Noffset, pars->norm);
  

  // example data stored in a pfresPatt (SHADOW) stucture
  fprintf (stderr, "x0: %.4f   I0: %.4f   n: %d   Npfres: %d   aa: %.1f  maxX00: %.1f  x00Step: %.1f  RStarProj: %.1f   AU: %.1f  offset: %.1f  maxOffset: %.1f\n", pfresPatt[0]->xI[0].x, pfresPatt[0]->xI[0].I, pfresPatt[0]->n, Npfres, pfresPatt[0]->aa, pfresPatt[0]->maxX00, pfresPatt[0]->x00Step, pfresPatt[0]->RStar, pfresPatt[0]->AU, pfresPatt[0]->offset, pfresPatt[0]->maxOffset);
  


  // a few points from a timeseries read in.
  fprintf (stderr, "t0: %.3f   I0: %.4f   dt: %.3f  hz: %.3f\n", pts->tI[0].t, pts->tI[0].I, pts->dt, pts->hz);



  // test the offsets  ////////////////////////////////////////////
  int i,j;
  FILE *fp;
  //  char *prefix = "tmpfres-";
  char filename[MAX_FILENAME];

  for (i=0; i<Npfres; i++) {

      snprintf (filename, MAX_FILENAME, "%s%05d","tmpfres-",i);
    fp = openfile (filename, "w");
    for (j=0; j<(pfresPatt[i]->n); j++) {
      fprintf(fp, "%.4f %.4f %.4f\n", pfresPatt[i]->xI[j].x, 
	      pfresPatt[i]->offset, pfresPatt[i]->xI[j].I);
    }
    closefile (fp);

  }
  // done testing offsets  ////////////////////////////////



  // testing makeKernels integrations  /////////////////////

  int i_vel;

  for (i=0; i<Npfres; i++) {

    for (i_vel=0; i_vel<N_VELOCITIES; i_vel++) {
        snprintf(filename, MAX_FILENAME, "%s%03d-%01d","kernel-",i,i_vel);
      
      fp = openfile (filename, "w");
      for (j=0; j<(pfresPatt[i]->nI[i_vel]); j++) {
	fprintf(fp, "%d %.4f\n", j, pfresPatt[i]->I[i_vel][j]);
      }
      closefile (fp);
    }
      
  }
  // done testing makeKernels integrations  ///////////////////////




  // testing addKBO  //////////////////////////////////////


  //for (i=0; i<ptsa->n; i++) {
  // printf ("%.4f %.4f\n", ptsa->tI[i].t, ptsa->tI[i].I);
  //}

  // done testing addKBO  /////////////////////////////////
}
