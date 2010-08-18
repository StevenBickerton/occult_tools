/*            */
/* Made with makeScript, Sat Feb 9 2008 */
/* Host: apache2 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "detection.h"

int main(int argc, char *argv[])
{
  int i;
  TS ts, tsa, *pts, *ptsa;
  pts = &ts;
  ptsa = &tsa;
  char *tsfile;
  FLOAT elong = PI;
  FLOAT incl = 0;
  
  // read in command line arguments
  if ( argc != 2 ) {
    fprintf(stderr,"Usage: %s timeseries\n", argv[0]);
    exit(EXIT_FAILURE);
  } else {
    tsfile = argv[1];
  }

  /* load timeseries from files */
  loadTimeseries(tsfile, pts, ptsa, elong, incl);

  /* print it */
  printf ("# fits2ascii %s\n", tsfile);
  for (i=0;i<pts->n;i++)
    printf ("%.8f %.8f\n", pts->tI[i].t, pts->tI[i].I);

  return(0);
}
  
