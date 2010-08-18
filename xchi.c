/*            */
/* Made with makeScript, Tue Jun 28, 2005  10:57:29 DST */
/* Host: kuiper */
/* Working Directory: /1/home/bickersj/sandbox/cdiffracSim  */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>

#include "detection.h"
#include "statistics.h"

/*  Main body of code */
int main(int argc, char *argv[])
{

  int i_vel;
  int nChiHit, icentre;
  HIT chihits[MAX_HITS];
  char *thisfile, *fresfile, *tsfile;
  time_t tseed = time(NULL);

  FLOAT elong, offset, chiThresh, incl;

  TS ts, *pts, tsa, *ptsa;
  pts = &ts;        // pure time series   (ptsa has added kbo occ events)
  ptsa = &tsa;

  SHADOW *pfresPatt[2];

  srand(tseed);   // initialize random number generator with current time.


  // read in command line arguments
  thisfile = argv[0];
  if ( argc != 7 ) {
    fprintf(stderr,"Usage: %s elong incl fresnelfile timeseries offset chiThresh\n", thisfile);
    exit(EXIT_FAILURE);
  } else {
    elong = RAD * atof(argv[1]);
    incl = RAD * atof(argv[2]);
    fresfile = argv[3];
    tsfile = argv[4];
    offset = atof(argv[5]);
    chiThresh = atof(argv[6]);
  }

  
  //  Load the data
  DPRINTF("main(): calling loadTimeseries()\n");
  loadTimeseries(tsfile, pts, ptsa, elong, incl);
  DPRINTF("main(): calling loadFresnelFile()\n");
  allocateFresPattMem(pfresPatt, 2);
  loadFresnelFile(pfresPatt[0], fresfile);

  DPRINTF ("main(): tsfile: %s  dt: %.3f\n", tsfile, pts->dt);

  DPRINTF ("main(): callin offsetFresPatt().\n");
  offsetFresPatt(pfresPatt[0], pfresPatt[1], offset);
  DPRINTF ("main(): calling initPfresPatt().\n");
  i_vel = 0;
  initPfresPatt(pfresPatt[1], ptsa, i_vel);
  DPRINTF ("main(): calling makeKernel().\n");
  makeKernel(NO_KERNEL, &icentre, pfresPatt[1], ptsa, i_vel, NO_RANDOM, OVERWRITE);


  // cross-correlate the time-series
  DPRINTF("main(): calling xcorrelate()\n");
  nChiHit = xchi(chihits, ptsa, pfresPatt[1], i_vel, chiThresh, DUMP);

  // cleanup
  free(pfresPatt[1]->I[i_vel]);

  return(0);
}
