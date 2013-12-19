/*            */
/* Made with makeScript, Tue Jun 28, 2005  10:57:29 DST */
/* Host: kuiper */
/* Working Directory: /1/home/bickersj/sandbox/cdiffracSim  */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <assert.h>
#include <getopt.h>

#include "detection.h"
#include "statistics.h"

/*   globals    */



/*  Main body of code */

int main(int argc, char *argv[])
{

  int i, j, l;
  int Npfres, Nfresfiles;
  int icentre, nCorrHit, nChiHit, nHit, total_hits=0;
  int file_i;
  HIT corrhits[MAX_HITS], chihits[MAX_HITS];
  char *paramfile, *tsfile;
  char hitsfile[MAX_FILENAME], lengthfile[MAX_FILENAME];
  char *thisfile;
  time_t tseed = time(NULL);

  HIT *allhits;
  allhits = det_malloc(MAX_HITS*sizeof(HIT));

  int nallhits = 0;

  PARAMLIST params;
  PARAMLIST *pars;
  pars = &params;
  TS ts, tsa, *pts, *ptsa;
  pts = &ts; // pure time series   (ptsa has added kbo occ events)
  ptsa = &tsa;

  FLOAT elong, incl;

  srand(tseed);   // init rand num generator with current time.


  char c;
  int useKO = 0;
  while ( (c=getopt(argc, argv, "k")) != -1 ) {
    switch (c) {
    case 'k':
      useKO = 1;
      break;
    default:
      break;
    }
  }

  // read in command line arguments
  DPRINTF("main(): Reading in Command line args.\n");
  thisfile = argv[0];
  if ( argc - optind != 4 ) {
    fprintf(stderr,"Usage: %s elongation incl paramFile timeseries\n", thisfile);
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

  // allocate memory for the pfresPatt structure
  DPRINTF ("main(): allocating pfresPatt memory.\n");
  SHADOW *pfresPatt[ (pars->Noffset+1) ];
  allocateFresPattMem(pfresPatt, pars->Noffset+1);

  SERIES chi;
  SERIES *pchi = &chi;
  pchi->value = det_malloc(sizeof(FLOAT)*pts->n);
  

  ///////////////////////////////////////////////////////////////////
  // main loop for adding KBOs and finding them
  ///////////////////////////////////////////////////////////////////

  DPRINTF ("main(): opening output files.\n");
  FILE *fp_hits; //, *fp_test;
  snprintf(hitsfile, MAX_FILENAME, "%s.hits", tsfile);
  fp_hits = openfile(hitsfile, "w");
  printHitHeader(fp_hits);

  FILE *fp_len; //, *fp_test;
  snprintf (lengthfile, MAX_FILENAME, "%s.len", tsfile);
  fp_len = openfile(lengthfile, "w");
  printLenHeader(fp_len);

  Npfres = 1;
  for (file_i=0; file_i<Nfresfiles; file_i++) {

    DPRINTF ("main(): Starting file loop\n");

    // load the next fresfile and make the offsets
    Npfres = 1;
    DPRINTF ("main(): calling loadFresnelFile().\n");
    loadFresnelFile(pfresPatt[0], fresfilelist[file_i]);
    DPRINTF ("main(): calling offsetFresPatts().\n");
    Npfres = offsetFresPatts(pfresPatt, pars, Npfres, AUTO_OFFSET);


    // make the convolution/xchi detection kernels
    DPRINTF ("main(): allocating memory and calling makeKernel().\n");
    for (i=0; i<N_VELOCITIES; i++) {
      for (j=0; j<Npfres; j++) {  // ignore return n value for makeKernel
	
	initPfresPatt(pfresPatt[j], ptsa, i);
	
	if (! pfresPatt[j]->is_all_zero) { // don't bother for 'nosignal'
	  makeKernel(NULL, &icentre, pfresPatt[j], pts, i, 
		     NO_RANDOM, OVERWRITE); 
	}

      }
    }


    for (i=0; i<N_VELOCITIES; i++) {      
      for (j=0; j<Npfres; j++) {
		
	// skip this one if the kernel had no signal in it (nothin to find)
	if (pfresPatt[j]->is_all_zero)
	  continue;
	
	DPRINTF ("main(): calling xcorrelate()\n");
	if ( (useKO)  &&  (ptsa->nKO > 0 ) ) {
	  nCorrHit = xcorrelateKO(corrhits,ptsa,pfresPatt[j],i,
				pars->corrThresh,NO_DUMP,pars->norm);
	} else { 
	  nCorrHit = xcorrelate(corrhits,ptsa,pfresPatt[j],i,
				pars->corrThresh,NO_DUMP,pars->norm);
	}
	
	//nChiHit  = xchi(chihits, ptsa, pfresPatt[j], i, pars->chiThresh);
	//nHit     = merge(corrhits, chihits, nCorrHit, nChiHit);
	DPRINTF ("main(): calling xchi2().\n");
	//nChiHit  = xchi2(chihits, ptsa, pfresPatt[j], i, pchi);
	nChiHit  = xchi3(chihits, ptsa, pfresPatt[j], i, pchi, corrhits, nCorrHit);
	DPRINTF ("main(): calling merge2().\n");
	nHit     = merge2(corrhits, chihits, nCorrHit, nChiHit, pchi, pars);
	total_hits += nHit;

	/* realloc allhits */
	allhits = det_realloc(allhits,(1+nallhits+nHit)*sizeof(HIT));
	/* put stuff in allhits */
	for (l=0;l<nHit;l++) {
	  allhits[nallhits+l].index = corrhits[l].index;
	}
	nallhits += nHit;

	// record anything found that we didn't add (detection or false pos.)
	for (l=0; l<nHit; l++) {
	  printHit(fp_hits, pfresPatt[j], corrhits[l], chihits[l]);
	}
	
	if (i==0  &&  pfresPatt[j]->offset == 0.0) {
	  printLen(fp_len, ptsa, pfresPatt[j]);
	}

      }    // end for j < Npfres
    }      // end for i < N_VELOCITIES
    
    fprintf (stderr, "%s %d/%d\n", 
	     fresfilelist[file_i], file_i+1, Nfresfiles);

    for (j=0; j<Npfres; j++) {
      for (i=0; i<N_VELOCITIES; i++) {
	free( pfresPatt[j]->I[i] );
      }
    }
    

  }  // end for file_i < Nfresfiles

  DPRINTF ("main(): closing output files.\n");
  closefile(fp_hits);
  
  fprintf (stderr, "\n");  

  /* ---  write a log file  ---- */
  FILE *fp_log; //, *fp_test;
  char logfile[MAX_FILENAME];
  snprintf (logfile, MAX_FILENAME, "%s.log", tsfile);
  fp_log = openfile(logfile, "w");
  fprintf(fp_log, "# tsfile  N   Length(s)   Hits   RMS  Kurtosis\n");
  fprintf(fp_log, "%s   %d  %8.1f   %d   %.4f  %.4f\n", 
	  tsfile, pts->n, pts->n*pts->dt, 
	  countHits(allhits,nallhits), pts->rms, pts->kurt);
  closefile(fp_log);

  // free the memory
  DPRINTF ("main(): freeing memory.\n");
  free(pchi->value);
  for (j=0; j<Npfres; j++) {
    free( pfresPatt[j]->xI );
    free( pfresPatt[j] );
  }
  free(allhits);
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
  
