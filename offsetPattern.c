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

/*   globals    */



/*  Main body of code */

int main(int argc, char *argv[])
{

  int i, j, nOffset;
  char *fresfile, *thisfile;
  time_t tseed = time(NULL);

  SHADOW *pfresPatt[2];

  FLOAT maxOffset, offset;


  srand(tseed);   // initialize random number generator with current time.


  // read in command line arguments
  DPRINTF("main(): Reading in Command line args.\n");
  thisfile = argv[0];
  if ( argc != 4 ) {
    fprintf(stderr,"Usage: %s fresfile maxOffset nOffset\n",
	    thisfile);
    exit(EXIT_FAILURE);
  } else {
    fresfile = argv[1];
    maxOffset = atof(argv[2]);
    nOffset = atoi(argv[3]);
  }



  DPRINTF ("main(): calling loadFresnelFile().\n");
  allocateFresPattMem(pfresPatt, 2);
  loadFresnelFile(pfresPatt[0], fresfile);

  DPRINTF ("main(): calling offsetFresPatt().\n");

  DPRINTF ("main(): printing header.\n");
  fprintf (stdout, "# %s \tKBO_radius: \t%g\n", LAA, pfresPatt[0]->aa);
  fprintf (stdout, "# %s \tlambda_Low: \t%g\n", LLAMBLO, pfresPatt[0]->lambLo);
  fprintf (stdout, "# %s \tlambda_High: \t%g\n", LLAMBHI, pfresPatt[0]->lambHi);
  fprintf (stdout, "# %s \tNum_lambdas: \t%d\n", LNLAMBDA, pfresPatt[0]->nLambda);
  fprintf (stdout, "# %s \tmax_dist_from_shadow_centre(m): \t%g\n", LMAXX00, pfresPatt[0]->maxX00);
  fprintf (stdout, "# %s \tstep_size: \t%g\n", LX00STEP, pfresPatt[0]->x00Step);
  fprintf (stdout, "# %s \tProjected_Star_Radius: \t%g\n", LRSTAR, pfresPatt[0]->RStar);
  fprintf (stdout, "# %s \tNum_point_sources_for_star: \t%d\n", LNSTARS, pfresPatt[0]->nStars);
  fprintf (stdout, "# %s \tdist_to_KBO_(AU): \t%g\n", LAU, pfresPatt[0]->AU);
  fprintf (stdout, "# %s \toffset_from_line_of_sight: \t%g\n", LOFFSET, pfresPatt[0]->offset);
  fprintf (stdout, "# %s \taspectRatio_of_KBO_occulter(x/y): \t%.2g\n", LASPRAT, pfresPatt[0]->aspRat);
  fprintf (stdout, "# %s order_of_approximation: \t%d\n", LORDER, pfresPatt[0]->order);

  for (i=0; i<=nOffset; i++) {

    offset = maxOffset/((FLOAT) nOffset);

    // if only one offset was requested, only return that one (not zero-offset)
    if (nOffset==1 && i==0)
      continue;
    
    offsetFresPatt(pfresPatt[0], pfresPatt[1], ((FLOAT) i)*offset);
    
    DPRINTF ("main(): outputting to file\n");
    for (j=0; j<pfresPatt[1]->n; j++) {
      fprintf (stdout, "%.1f %.1f %.8f\n", 
	       ((FLOAT) i)*offset, pfresPatt[1]->xI[j].x, pfresPatt[1]->xI[j].I);
    }
    fprintf (stdout, "\n");

  }
  
  DPRINTF ("main(): exiting\n");
  return(0);
}
  
