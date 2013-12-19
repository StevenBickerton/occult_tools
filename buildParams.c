/*            */
/* Steven Bickerton */
/* Dept. of Physics/Astronomy, McMaster University */
/* bick@physics.mcmaster.ca*/
/* Made with makeScript, Mon Nov  6, 2006  11:47:50 EST */
/* Host: kuiper */
/* Working Directory: /1/home/bickersj/pack_main/sandbox/cdiffracSim  */


#include <stdlib.h>
#include <stdio.h>

#include "detection.h"
#include "statistics.h"

#define MAX_PARAM 25

void usage(char execName[]) {
    printf ("Usage: %s\n", execName);
    exit(1);
}

void ddump(char param[], int value) { 
  char format[MAX_LINE_LENGTH];
  snprintf (format, MAX_LINE_LENGTH,"%%-%ds %%d\n", MAX_PARAM);
  printf (format, param, value);}
void fdump(char param[], FLOAT value) { 
  char format[MAX_LINE_LENGTH];
  snprintf (format, MAX_LINE_LENGTH,"%%-%ds %%.4f\n", MAX_PARAM);
  printf (format, param, value);}
void sdump(char param[], char value[]) { 
  char format[MAX_LINE_LENGTH];
  snprintf (format, MAX_LINE_LENGTH,"%%-%ds %%s\n", MAX_PARAM);
  printf (format, param, value);}

int main(int argc, char *argv[])
{

  printf ("\ndetection.h preprocessor constants\n\n");

  
  ddump ("MAX_PATTERN_POINTS", MAX_PATTERN_POINTS);
  ddump ("MAX_LINE_LENGTH", MAX_LINE_LENGTH);
  ddump ("MAX_FILENAME", MAX_FILENAME);
  ddump ("MAX_FRESFILES", MAX_FRESFILES);
  ddump ("MAX_KERNEL_OVERSAMPLING",MAX_KERNEL_OVERSAMPLING );

  // allow for a high precision option which keeps small signals
  //  values below these limits will be ignored by the code

#ifdef PRECISION
  sdump ("PRECISION", "enabled");
#else
  sdump ("PRECISION", "disabled");
#endif

  fdump ("OFFSET_MIN_I_LIMIT", OFFSET_MIN_I_LIMIT);
  fdump ("KERNEL_LO_I_CUT", KERNEL_LO_I_CUT);        
  // I smaller than this in kernel ignored

  ddump ("MAX_OFFSETS", MAX_OFFSETS); 
  ddump ("N_PAST_MIN_I_LIMIT", N_PAST_MIN_I_LIMIT); 
  // ddump ("MAX_KERNEL_POINTS ", MAX_KERNEL_POINTS );
  //  ddump ("MAX_TS_POINTS", MAX_TS_POINTS);
  ddump ("N_VELOCITIES ", N_VELOCITIES );  // total num velocites to test (odd)
  fdump ("T_STEP", T_STEP  );              // time resol to use for integral
  ddump ("ITMP_SIZE", ITMP_SIZE);        // size of array Itmp in makeKernels()


  ddump ("KERNEL_N_PAST_MIN_I_LIMIT", KERNEL_N_PAST_MIN_I_LIMIT);   
  // n points outside 'signal' to include

  ddump ("N_SEEDS", N_SEEDS);           // number of seeds to keep stored
  ddump ("MAX_HITS", MAX_HITS );        // max size of hits() array

  fdump ("SIGMA_CLIP", SIGMA_CLIP ); // sigma clipping level to use for stats
  ddump ("SUBSAMPLES", SUBSAMPLES ); // minimum no. of samples to use for stats
  ddump ("HIT_WIDTH", HIT_WIDTH   ); // dist. betw. 1st and last above for 1 hit

  sdump ("FRESFILE_PREFIX", FRESFILE_PREFIX);
  ddump ("FRESFILE_LENGTH", FRESFILE_LENGTH ); // characters in fresfile pattern

  ddump ("NORMALIZING", NORMALIZING);   // to normalize the xcorr values or not

  ddump ("XCORR_RMS_WINDOW", XCORR_RMS_WINDOW);  
  // half-wid of sliding RMS used
  // to normalize the xcorr output.

  fdump ("XCORR_CORRECTION", XCORR_CORRECTION);  
  // factor to shift the ideal xcorr value
  //   to one consistent with data

  fdump ("XCHI_CORRECTION", XCHI_CORRECTION);         
  // factor to shift the ideal xchi value

  fdump ("XCORR_RMS_RANGE", XCORR_RMS_RANGE );         
  // rms range for xcor acceptance of a det.

  ddump ("PROXIMITY", PROXIMITY);                  
  // distance for a detection to be ok


  printf ("\nstatistics.h preprocessor constants\n\n");

  ddump("STAT_MAX_LINE_LENGTH", STAT_MAX_LINE_LENGTH);
  ddump("CLIP_TOLERANCE", CLIP_TOLERANCE); 
  ddump("CLIP_ITERATIONS", CLIP_ITERATIONS);

  printf ("\n");
  return 0;
    
}


