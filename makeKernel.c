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

  int i, icentre;
  int centering, random;
  char *fresfile, *thisfile;
  time_t tseed = time(NULL);

  TS tsa, *ptsa;
  ptsa = &tsa;

  SHADOW *pfresPatt[2];

  FLOAT elong, dt, t, offset, incl;


  srand(tseed);   // initialize random number generator with current time.


  // read in command line arguments
  DPRINTF("main(): Reading in Command line args.\n");
  thisfile = argv[0];
  if ( argc != 8 ) {
    fprintf(stderr,"Usage: %s fresfile elong incl dt offset random[0|1] centering[0|1] \n",
	    thisfile);
    exit(EXIT_FAILURE);
  } else {
    fresfile =  argv[1];
    elong = RAD * atof(argv[2]);
    incl = RAD * atof(argv[3]);
    dt = atof(argv[4]);
    offset = atof(argv[5]);
    random = atoi(argv[6]);
    centering = atoi(argv[7]);
  }
  ptsa->elong = elong;
  ptsa->incl = incl;
  ptsa->dt = dt;


  // allocate memory for the pfresPatt structure
  DPRINTF ("main(): allocating pfresPatt memory.\n");

  DPRINTF ("main(): calling loadFresnelFile().\n");
  allocateFresPattMem(pfresPatt,2);
  loadFresnelFile(pfresPatt[0], fresfile);

  DPRINTF ("main(): calling offsetFresPatt().\n");
  if (!pfresPatt[0]->fullprofile) {
    offsetFresPatt(pfresPatt[0], pfresPatt[1], offset);
  } else {
    pfresPatt[1] = pfresPatt[0];  // no op
  }

  DPRINTF ("main(): initiallizing the fresnel pattern.\n");
  initPfresPatt(pfresPatt[1], ptsa, 0);  
  DPRINTF ("main(): making the Kernel\n");
  makeKernel(NULL, &icentre, pfresPatt[1], ptsa, 0, random, OVERWRITE); 

  DPRINTF ("main(): outputting to file\n");
  printf ("# x_start = %.2f\n", pfresPatt[1]->x_start);
  printf ("# t_centre = %.4f\n", pfresPatt[1]->t_centre);
  printf ("# t(sec) I/<I>\n");
  for (i=0; i<pfresPatt[1]->nI[0]; i++) {
    t = (centering) ? ( (i-icentre)*dt )  : (i*dt);
    printf ("%.8f %.8f\n", t, pfresPatt[1]->I[0][i]);
  }

  // free the memory
  DPRINTF ("main(): freeing memory.\n");
  free( pfresPatt[1]->I[0] );

  DPRINTF ("main(): exiting\n");
  return(0);
}
  
