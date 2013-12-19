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

  int i, N;
  int i_vel, centering, random;
  char *tsfile, outfile[MAX_FILENAME];
  char *thisfile, *fresfile;
  time_t tseed = time(NULL);
  TS ts, *pts, tsa, *ptsa;
  pts = &ts;        // pure time series   (ptsa has added kbo occ events)
  ptsa = &tsa;
  FLOAT elong, incl, offset;
  FILE *fp;
  SHADOW *pfresPatt[2];
  unsigned int seed = (unsigned int) tseed;

  srand(tseed);   // initialize random number generator with current time.

  // read in command line arguments
  thisfile = argv[0];
  if ( argc < 9 || argc > 10) {
    fprintf(stderr,"Usage: %s elong incl fresnelfile timeseries N offset center[1|0] random[1|0] [seed]\n", thisfile);
    exit(EXIT_FAILURE);
  } else {
    elong = RAD * atof(argv[1]);
    incl  = RAD * atof(argv[2]);
    fresfile = argv[3];
    tsfile = argv[4];
    N = atoi(argv[5]);
    offset = atof(argv[6]);
    centering = atoi(argv[7]);
    random = atoi(argv[8]);
    if (argc == 10)
      seed = atoi(argv[9]);
  }

  if (seed > 0)
    srand( (unsigned int) seed);


  int addlist[N];
  
  //  Load the data
  DPRINTF("main(): calling loadTimeseries()\n");
  loadTimeseries(tsfile, pts, ptsa, elong, incl);
  DPRINTF("main(): calling loadFresnelFile()\n");
  allocateFresPattMem(pfresPatt, 2);
  loadFresnelFile(pfresPatt[0], fresfile);
  DPRINTF ("main(): callin offsetFresPatt().\n");
  offsetFresPatt(pfresPatt[0], pfresPatt[1], offset);
  DPRINTF ("main(): calling initPfresPatt().\n");
  i_vel = 0;
  initPfresPatt(pfresPatt[1], ptsa, i_vel);


  // Add the KBO to the timeseries
  DPRINTF("main(): calling addKBO()\n");
  addKBO(addlist, ptsa, pfresPatt[1], N, i_vel, centering, random);


  // output the added timeseries

  if (ptsa->fits) {

    /* quick and dirty - use writeXcorrSeries */
    HIT hits[N];
    for (i=0;i<N;i++) {
      hits[i].index   = addlist[i];
      hits[i].cormag  = 0;
      hits[i].cormag0 = 0;
    }
    FLOAT *Itmp;
    Itmp = det_malloc(ptsa->n*sizeof(FLOAT));
    for (i=0;i<ptsa->n;i++) {
      Itmp[i] = ptsa->tI[i].I;
    }
    writeXcorrSeries(ptsa, pfresPatt[1], NULL, Itmp, hits, N, "add");

  } else {
    DPRINTF("main(): output to file\n");
    snprintf (outfile, MAX_FILENAME,"%s.add", tsfile);
    fp = openfile (outfile, "w");
    
    // a header
    fprintf (fp, "# Created by addKBO\n");
    fprintf (fp, "# Fresnel_Pattern_File: %s\n", fresfile);
    fprintf (fp, "# elongation: %.1f\n", 180.0*elong/PI);
    fprintf (fp, "# Number_added: %d\n", N);
    fprintf (fp, "# %s \tKBO_radius: \t%g\n", LAA, pfresPatt[1]->aa);
    fprintf (fp, "# %s \tdist_to_KBO_(AU): \t%g\n", LAU, pfresPatt[1]->AU);
    fprintf (fp, "# %s \toffset_from_line_of_sight: \t%g\n", LOFFSET, offset);
    fprintf (fp, "# %s \tlambda_Low: \t%g\n", LLAMBLO, pfresPatt[1]->lambLo);
    fprintf (fp, "# %s \tlambda_High: \t%g\n", LLAMBHI, pfresPatt[1]->lambHi);
    fprintf (fp, "# %s \tNum_lambdas: \t%d\n", LNLAMBDA, pfresPatt[1]->nLambda);
    fprintf (fp, "# %s \tProjected_Star_Radius: \t%g\n", LRSTAR, pfresPatt[1]->RStar);
    fprintf (fp, "# %s \tNum_point_sources_for_star: \t%d\n", LNSTARS, pfresPatt[1]->nStars);
    
    fprintf (fp, "# No.  Index  time\n");
    for (i=0; i<N; i++)
      fprintf (fp, "# %02d %d %.4f\n", i, addlist[i], addlist[i]*ptsa->dt);
    
    // the added timeseries
    for (i=0;i<pts->n;i++) 
      fprintf (fp, "%.4f %.8f\n", ptsa->tI[i].t, ptsa->tI[i].I);
    
    closefile(fp);
    
  }

  // cleanup
  free(pfresPatt[1]->I[i_vel]);

  return(0);
}
