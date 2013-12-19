/*            */
/* Made with makeScript, Tue Jun 28, 2005  10:57:29 DST */
/* Host: kuiper */
/* Working Directory: /1/home/bickersj/sandbox/cdiffracSim  */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <getopt.h>
#include <assert.h>

#include "detection.h"
#include "statistics.h"

/*   globals    */

static void usage() {
  fprintf(stdout, "Usage: komplete [options] elongation incl paramFile timeseries\n");
  fprintf(stdout, "  -h            = this message\n");
  fprintf(stdout, "  -r readnoise  = noise added after event planting\n");
  fprintf(stdout, "  -x maxOffset  = override 2.5Fsu offset of kernels\n");
}

/*  Main body of code */

int main(int argc, char *argv[])
{

  int i, j, k, l;
  int Npfres, Nfresfiles;
  int icentre, nCorrHit, nChiHit, nHit;
  int file_i;
  HIT corrhits[MAX_HITS], chihits[MAX_HITS], recovered[MAX_HITS];
  char *paramfile, *tsfile;
  char hitsfile[MAX_FILENAME], statsfile[MAX_FILENAME];
  char *thisfile;
  time_t tseed = time(NULL);
  STAT stats;
  
  char c;
  FLOAT rdnoise = 0.0;
  FLOAT max_offset = 0.0;

  PARAMLIST params;
  PARAMLIST *pars;
  pars = &params;
  TS ts, tsa, *pts, *ptsa;
  pts = &ts;        // pure time series   (ptsa has added kbo occ events)
  ptsa = &tsa;
  FLOAT elong, incl;


  srand(tseed);   // initialize random number generator with current time.


  // read in command line arguments

  while ( (c=getopt(argc, argv, "hr:x:")) != -1 ) {
    switch(c) {
    case 'h':
      usage();
      exit(EXIT_FAILURE);
      break;
    case 'r':
      rdnoise = atof(optarg);
      break;
    case 'x':
      max_offset = atof(optarg);
      break;
    default:
      break;
    }
  }

  DPRINTF("main(): Reading in Command line args.\n");
  thisfile = argv[0];
  if ( argc-optind != 4 ) {
    usage();
    exit(EXIT_FAILURE);
  } else {
    elong = RAD * atof(argv[optind]);
    incl = RAD * atof(argv[optind+1]);
    paramfile = argv[optind+2];
    tsfile = argv[optind+3];
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

  int addlist[pars->NperCycle];
  SERIES chi;
  SERIES *pchi = &chi;
  pchi->value = det_malloc(sizeof(FLOAT)*pts->n);
  for ( i=0; i<pts->n; i++ ) { pchi->value[i] = 0.0; }

  // '2' req'd as sometimes more values are found than were added
  //  --> fringes show up as hits, or hits coincidently overlap
  FLOAT cormags[2 * pars->NperCycle * pars->cycles];
  FLOAT chimags[2 * pars->NperCycle * pars->cycles];
  FILE *fp_hits, *fp_stats; //, *fp_test;


  ///////////////////////////////////////////////////////////////////
  // main loop for adding KBOs and finding them
  ///////////////////////////////////////////////////////////////////

  DPRINTF ("main(): opening output files.\n");
  snprintf(hitsfile, MAX_FILENAME, "%s.hits", tsfile);
  fp_hits = openfile(hitsfile, "w");
  printHitHeader(fp_hits);

  snprintf(statsfile, MAX_FILENAME,"%s.stats", tsfile);
  fp_stats = openfile(statsfile, "w");
  printStatsHeader(fp_stats);

  Npfres = 1;
  for (file_i=0; file_i<Nfresfiles; file_i++) {

    DPRINTF ("main(): Starting file loop\n");

    // allocate memory for the pfresPatt structure
    DPRINTF ("main(): allocating pfresPatt memory.\n");
    SHADOW *pfresPatt[ (pars->Noffset+1) ];
    allocateFresPattMem(pfresPatt, pars->Noffset+1);
    
    // load the next fresfile and make the offsets
    Npfres = 1;
    DPRINTF ("main(): calling loadFresnelFile().\n");
    loadFresnelFile(pfresPatt[0], fresfilelist[file_i]);
    pfresPatt[0]->is_all_zero = 0; // want to try it, even if zero
    if ( max_offset ) { pfresPatt[0]->maxOffset = max_offset; }
    DPRINTF ("main(): callin offsetFresPatts().\n");
    Npfres = offsetFresPatts(pfresPatt, pars, Npfres, AUTO_OFFSET);


    // make the convolution/xchi detection kernels
    DPRINTF ("main(): allocating memory and calling makeKernel().\n");
    for (j=0; j<Npfres; j++) {  // ignore return n value for makeKernel
      for (i=0; i<N_VELOCITIES; i++) {

	initPfresPatt(pfresPatt[j], ptsa, i);

	if (! pfresPatt[j]->is_all_zero) { // don't bother for 'nosignal'
	  makeKernel(NULL, &icentre, pfresPatt[j], pts, 
		     i, NO_RANDOM, OVERWRITE); 
	}

      }
    }

    for (j=0; j<Npfres; j++) {
      for (i=0; i<N_VELOCITIES; i++) {      
	
	pfresPatt[j]->nRecover[i] = 0;
	for (k=0; k<(pars->cycles); k++) {
	  
	  // skip this one if the kernel had no signal in it (nothin to find)
	  //printf("zero: %d\n", pfresPatt[j]->is_all_zero);
	  if (pfresPatt[j]->is_all_zero) {
	    continue;
	  }
	  
	  // reset the timeseries intensity values
	  for (l=0; l<pts->n; l++) {
	    ptsa->tI[l].I = pts->tI[l].I;
	  }
	  DPRINTF ("main(): calling addKBO()  i=%d  j=%d  k=%d.\n", i, j, k);
	  addKBO(addlist, ptsa, pfresPatt[j], pars->NperCycle, i, NO_CENTRE, RANDOM);


	  if (rdnoise) {
	    addReadnoise(rdnoise, ptsa);
	  }
	  
	  DPRINTF ("main(): calling xcorrelate()\n");
	  nCorrHit = xcorrelate(corrhits,ptsa,pfresPatt[j],i,
				pars->corrThresh,NO_DUMP, pars->norm);

	  // need to reset this for newly added events
	  if (ptsa->FFTisdone) {
	    free(ptsa->fI);
	    ptsa->FFTisdone = 0;
	  }

	  //nChiHit  = xchi(chihits, ptsa, pfresPatt[j], i, pars->chiThresh);
	  //nHit     = merge(corrhits, chihits, nCorrHit, nChiHit);
	  DPRINTF ("main(): calling xchi2().\n");
	  //nChiHit  = xchi3(chihits, ptsa, pfresPatt[j], i, pchi, corrhits, nCorrHit);
	  nChiHit  = xchi2(chihits, ptsa, pfresPatt[j], i, pchi);
	  DPRINTF ("main(): calling merge2().\n");
	  nHit     = merge2(corrhits, chihits, nCorrHit, nChiHit, pchi, pars);
	  DPRINTF ("main(): calling addVerify().\n"); 

	  pfresPatt[j]->nRecover[i] = addVerify(addlist,pars->NperCycle, corrhits, chihits, &nHit, recovered, pfresPatt[j], ptsa);
	  
	  // record anything found that we didn't add (detection or false pos.)
	  for (l=0; l<nHit; l++) {
	    printHit(fp_hits, pfresPatt[j], corrhits[l], chihits[l]);
	  }
	  
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
	printStats(fp_stats, pfresPatt[j], corrhits[0], chihits[0], pars, i);


	free(pfresPatt[j]->I[i]);
      }      // end for i < N_VELOCITIES
      free(pfresPatt[j]->xI);
      free(pfresPatt[j]);
    }    // end for j < Npfres

    fprintf (stderr, "%s %d/%d\n", 
	     fresfilelist[file_i], file_i+1, Nfresfiles);

  }  // end for file_i < Nfresfiles

  DPRINTF ("main(): closing output files.\n");
  closefile(fp_hits);
  closefile(fp_stats);
  
  fprintf (stderr, "\n");
  
  // have a look at output for debug purposes
  // debugOutput(pfresPatt, pars, pts, ptsa, Npfres);
  
  // free the memory
  free(pchi->value);
  //DPRINTF ("main(): freeing memory.\n");
  //for (j=0; j<Npfres; j++) {
  //  free( pfresPatt[j] );
  //}

  freeFresnelFileList(fresfilelist, Nfresfiles);

  free(pts->tI);
  free(ptsa->tI);
  if (pts->FFTisdone) {
    free( pts->fI );
  }
  if (ptsa->FFTisdone) {
    free( ptsa->fI );
  }


  DPRINTF ("main(): exiting\n");
  return(0);
}
  
