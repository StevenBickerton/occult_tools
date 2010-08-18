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

  int i;
  int Nfresfiles;
  int file_i;
  char *tsfile, *thisfile;
  char outfile[MAX_FILENAME];
  time_t tseed = time(NULL);
  //FLOAT offset;
  int i_vel = 0;
  FILE *fp;

  PARAMLIST params;
  PARAMLIST *pars;
  pars = &params;
  TS ts, tsa, *pts, *ptsa;
  pts = &ts;        // pure time series   (ptsa has added kbo occ events)
  ptsa = &tsa;
  FLOAT elong, incl;


  srand(tseed);   // initialize random number generator with current time.


  // read in command line arguments
  DPRINTF("main(): Reading in Command line args.\n");
  thisfile = argv[0];
  if ( argc != 5 ) {
    fprintf(stderr,"Usage: %s elongation incl fresDir timeseries\n", thisfile);
    exit(EXIT_FAILURE);
  } else {
    elong = RAD * atof(argv[1]);
    incl  = RAD * atof(argv[2]);
    sprintf (pars->fresDir, "%s", argv[3]);
    tsfile = argv[4];
  }
  pars->NperCycle = 1;

  ///////////////////////////////////////////////////////////
  // load parameters and timeseries from files
  DPRINTF ("main(): calling loadTimeseries().\n");
  loadTimeseries(tsfile, pts, ptsa, elong, incl);



  //////////////////////////////////////////////////////////
  // prepare to load the fresnel diffraction patterns
  DPRINTF ("main(): calling getFresnelFileList().\n");
  char *fresfilelist[MAX_FRESFILES];
  Nfresfiles = getFresnelFileList(fresfilelist, pars);

  // allocate memory for the pfresPatt structure
  DPRINTF ("main(): allocating pfresPatt memory.\n");
  SHADOW *pfresPatt[1];
  allocateFresPattMem(pfresPatt,1);

  int addlist[pars->NperCycle];


  ///////////////////////////////////////////////////////////////////
  //  main loop for adding KBOs
  ///////////////////////////////////////////////////////////////////

  // open the output file and print a header
  DPRINTF("main(): output to file\n");
  sprintf (outfile, "%s.hide", tsfile);
  fp = openfile (outfile, "w");

  for (file_i=0; file_i<Nfresfiles; file_i++) {

    DPRINTF ("main(): Starting file loop: file %d\n", file_i);


    // load the next fresfile but *don't* make the offsets
    DPRINTF ("main(): calling loadFresnelFile().\n");
    loadFresnelFile(pfresPatt[0], fresfilelist[file_i]);

    DPRINTF ("main(): allocating memory and calling makeKernel().\n");
    initPfresPatt(pfresPatt[0], ptsa, i_vel);

    // print a headre
    if (file_i == 0) {
      fprintf (fp, "# Created by hideKBO\n");
      fprintf (fp, "# %s \telongation: %.1f\n", LELONG, 180.0*elong/PI);
      fprintf (fp, "# %s \tlambda_Low: \t%g\n", LLAMBLO, pfresPatt[0]->lambLo);
      fprintf (fp, "# %s \tlambda_High: \t%g\n", LLAMBHI, pfresPatt[0]->lambHi);
      fprintf (fp, "# %s \tNum_lambdas: \t%d\n", LNLAMBDA, pfresPatt[0]->nLambda);
      fprintf (fp, "# %s \tTeff: \t%g\n", LTEFF, pfresPatt[0]->Teff);
      fprintf (fp, "# %s \tProjected_Star_Radius: \t%g\n", LRSTAR, pfresPatt[0]->RStar);
      fprintf (fp, "# %s \tNum_point_sources_for_star: \t%d\n", LNSTARS, pfresPatt[0]->nStars);
      fprintf (fp, "# Number_added: %d\n", Nfresfiles);
      fprintf (fp, "# \n");
      fprintf (fp, "# No.  Index  time   rad   AU   offset  AR   angle\n");
    }



    // add the KBO    
    DPRINTF ("main(): calling addKBO()\n");
    addKBO(addlist, ptsa, pfresPatt[0], pars->NperCycle, i_vel,NO_CENTRE,RANDOM);

    // print info about the event to a header
    DPRINTF ("main(): outputting event records to file\n");
    for (i=0; i<pars->NperCycle; i++)
      fprintf (fp, "# %3d %6d %.4f   %5.0f %5.0f %5.0f  %6.2f %5.0f\n", 
	       file_i, addlist[i], addlist[i]*ptsa->dt, 
	       pfresPatt[0]->aa, pfresPatt[0]->AU, pfresPatt[0]->offset, 
	       pfresPatt[0]->aspRat, pfresPatt[0]->angle);

  }  // end for file_i < Nfresfiles

  // print the added timeseries
  DPRINTF ("main(): printing the new timeseries to file\n");
  for (i=0;i<pts->n;i++) 
    fprintf (fp, "%.4f %.4f\n", ptsa->tI[i].t, ptsa->tI[i].I);
  
  DPRINTF ("main(): closing file handle\n");
  closefile(fp);
  

  // free the memory
  DPRINTF ("main(): freeing memory.\n");
  free( pfresPatt[0]->I[i_vel] );

  DPRINTF ("main(): exiting\n");
  return(0);
}
  
