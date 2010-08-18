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

  int i, j, k, l, q;
  int Npfres, Nfresfiles;
  int icentre, nCorrHit, nChiHit, nHit;
  //int file_i;
  HIT corrhits[MAX_HITS], chihits[MAX_HITS], recovered[MAX_HITS];
  char *paramfile, *tsfile, *thisfile;
  char hitsfile[MAX_FILENAME], statsfile[MAX_FILENAME];
  time_t tseed = time(NULL);
  STAT stats;

  PARAMLIST params;
  PARAMLIST *pars;
  pars = &params;
  TS ts, tsa, *pts, *ptsa;
  pts = &ts;        // pure time series   (ptsa has added kbo occ events)
  ptsa = &tsa;
  FLOAT elong, offset, incl;


  srand(tseed);   // initialize random number generator with current time.


  // read in command line arguments
  DPRINTF("main(): Reading in Command line args.\n");
  thisfile = argv[0];
  if ( argc != 6 ) {
    fprintf(stderr,"Usage: %s elongation incl paramFile timeseries offset\n", thisfile);
    exit(EXIT_FAILURE);
  } else {
    elong = RAD * atof(argv[1]);
    incl = RAD * atof(argv[2]);
    paramfile = argv[3];
    tsfile = argv[4];
    offset = atof(argv[5]);
  }


  ///////////////////////////////////////////////////////////
  // load parameters and timeseries from files
  DPRINTF ("main(): calling loadParameters().\n");
  loadParameters(paramfile, pars);  
  DPRINTF ("main(): calling loadTimeseries().\n");
  loadTimeseries(tsfile, pts, ptsa, elong, incl);



  //////////////////////////////////////////////////////////
  // prepare to load the fresnel diffraction patterns
  DPRINTF ("main(): calling getFresnelFileList().\n");
  char *fresfilelist[MAX_FRESFILES];
  Nfresfiles = getFresnelFileList(fresfilelist, pars);

  // allocate memory for the pfresPatt structure
  DPRINTF ("main(): allocating pfresPatt memory.\n");
  //SHADOW *pfresPatt[ (pars->Noffset+1) ];
  //allocateFresPattMem(pfresPatt, pars->Noffset+1);
  SHADOW *pfresPatt[ ((Nfresfiles)*(pars->Noffset+1)) ];
  allocateFresPattMem(pfresPatt, ((Nfresfiles)*(pars->Noffset+1)) );
  
  
  int addlist[pars->NperCycle];
  SERIES chi;
  SERIES *pchi = &chi;
  pchi->value = det_malloc(sizeof(FLOAT)*pts->n);
  
  FLOAT cormags[pars->NperCycle * pars->cycles];
  FLOAT chimags[pars->NperCycle * pars->cycles];
  FILE *fp_hits, *fp_stats; //, *fp_test;
  
  
  ///////////////////////////////////////////////////////////////////
  // main loop for adding KBOs and finding them
  ///////////////////////////////////////////////////////////////////
  
  DPRINTF ("main(): opening output files.\n");
  sprintf(hitsfile, "%s.Bhits", tsfile);
  fp_hits = openfile(hitsfile, "w");
  printHitHeader(fp_hits);
  
  sprintf(statsfile, "%s.Bstats", tsfile);
  fp_stats = openfile(statsfile, "w");
  printBStatsHeader(fp_stats);
  
  //for (file_i=0; file_i<Nfresfiles; file_i++) {
    
  //DPRINTF ("main(): Starting file loop\n");

  // load the next fresfile and make the offsets
  Npfres = 1;
  DPRINTF ("main(): calling loadFresnelFile().\n");
  loadFresnelFiles(pfresPatt, pars);
  DPRINTF ("main(): callin offsetFresPatts().\n");
  Npfres = offsetFresPatts(pfresPatt, pars, Nfresfiles, offset);
  
  
  // make the convolution/xchi detection kernels
  DPRINTF ("main(): allocating memory and calling makeKernel().\n");
  for (i=0; i<N_VELOCITIES; i++) {
    for (j=0; j<Npfres; j++) {  // ignore return n value for makeKernel
      
      initPfresPatt(pfresPatt[j], ptsa, i);
      
      if (! pfresPatt[j]->is_all_zero)  // don't bother for 'nosignal'
	makeKernel(NULL, &icentre, pfresPatt[j], pts, 
		   i, NO_RANDOM, OVERWRITE); 
      
    }
  }
  
  for (i=0; i<N_VELOCITIES; i++) {      
    for (j=0; j<Npfres; j++) {
      
      for (q=0; q<Npfres; q++) {
	
	pfresPatt[j]->nRecover[i] = 0;
	for (k=0; k<(pars->cycles); k++) {
	  
	  // skip this one if the kernel had no signal in it (nothin to find)
	  if (pfresPatt[j]->is_all_zero)
	    continue;
	  
	  // reset the timeseries intensity values
	  for (l=0; l<pts->n; l++)
	    ptsa->tI[l].I = pts->tI[l].I;
	  
	  DPRINTF ("main(): calling addKBO()  i=%d  j=%d  k=%d.\n", i, j, k);
	  addKBO(addlist, ptsa, pfresPatt[q], pars->NperCycle, i, NO_CENTRE, RANDOM);
	  ptsa->FFTisdone = 0; // reset for newly added objects
	  
	  DPRINTF ("main(): calling xcorrelate()\n");
	  nCorrHit = xcorrelate(corrhits,ptsa,pfresPatt[j],i,pars->corrThresh, NO_DUMP, pars->norm);
	  
	  //nChiHit  = xchi(chihits, ptsa, pfresPatt[j], i, pars->chiThresh);
	  //nHit     = merge(corrhits, chihits, nCorrHit, nChiHit);
	  DPRINTF ("main(): calling xchi2().\n");
	  nChiHit  = xchi2(chihits, ptsa, pfresPatt[j], i, pchi);
	  DPRINTF ("main(): calling merge2().\n");
	  nHit     = merge2(corrhits, chihits, nCorrHit, nChiHit, pchi, pars);
	  DPRINTF ("main(): calling addVerify().\n"); 
	  pfresPatt[j]->nRecover[i] = addVerify(addlist,pars->NperCycle, corrhits, chihits, &nHit, recovered, pfresPatt[j], ptsa);
	  
	  
	  // record anything found that we didn't add (detection or false pos.)
	  for (l=0; l<nHit; l++)
	    printHit(fp_hits, pfresPatt[j], corrhits[l], chihits[l]);
	  
	}   // end  for k < cycles
	
	DPRINTF ("main(): loop over k-cycles done\n");
	
	// get stats on the detection strengths of the recovered objects
	///  need to put values in an array to call stats functions
	for (k=0; k<(pfresPatt[j]->nRecover[i]); k++) {
	  cormags[k] = recovered[k].cormag;
	  chimags[k] = recovered[k].chimag;
	}
	
	
	DPRINTF ("main(): Starting calculating mean and rms  (n=%d)\n", pfresPatt[j]->nRecover[i]);
	stats = meanRms_d(cormags, pfresPatt[j]->nRecover[i], pfresPatt[j]->nRecover[i]);
	pfresPatt[j]->meanCor[i] = stats.mean;
	pfresPatt[j]->rmsCor[i] = stats.rms;

	stats = meanRms_d(chimags, pfresPatt[j]->nRecover[i], pfresPatt[j]->nRecover[i]);
	pfresPatt[j]->meanChi[i] = stats.mean;
	pfresPatt[j]->rmsChi[i] = stats.rms;

	// write a summary to the stats file
	DPRINTF ("main(): calling printStats().\n");
	printBStats(fp_stats, pfresPatt[j], pfresPatt[q], corrhits[0], chihits[0], pars, i);
	
      }  // end for q < Npfres

      fprintf (stderr, "done %d/%d\n", j+1, Npfres);

    }    // end for j < Npfres
  }      // end for i < N_VELOCITIES
  
  
  //}  // end for file_i < Nfresfiles
  
  DPRINTF ("main(): closing output files.\n");
  closefile(fp_hits);
  closefile(fp_stats);
  
  fprintf (stderr, "\n");
  
  // have a look at output for debug purposes
  // debugOutput(pfresPatt, pars, pts, ptsa, Npfres);
  
  // free the memory
  DPRINTF ("main(): freeing memory.\n");
  for (j=0; j<Npfres; j++) {
    for (i=0; i<N_VELOCITIES; i++) 
      free( pfresPatt[j]->I[i] );
    free( pfresPatt[j] );
  }

  DPRINTF ("main(): exiting\n");
  return(0);
}
  
