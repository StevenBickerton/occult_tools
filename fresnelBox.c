/*            */
/* Made with makeScript, Tue Jun 21, 2005  15:03:06 DST */
/* Host: kuiper */
/* Working Directory: /1/home/bickersj/sandbox/cdiffracSim  */
/*                                                          */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

#include "detection.h"


#define MAXLINE 128*sizeof(char)
#define MAX_NLAMBDA 1024
#define MAX_FRES_INT 40000

int main(int argc, char *argv[])
{


  // initializations 
  //  char paramfile[MAX_FILENAME], thisfile[MAX_FILENAME], paramLine[80];
  char *thisfile, paramLine[80];
  char paramfile[MAX_FILENAME];
  FLOAT value_tmp;

  FLOAT offset = 0;
  FLOAT RStarProj = 100;
  FLOAT AU  = 100.0;
  int order = 5;
  int Nptsrc = 1;
  FLOAT maxX00 = 20000;
  FLOAT x00Step = 50;
  FLOAT lambLo = 4.0e-7;
  FLOAT lambHi = 7.0e-7;
  int nLambda = 1;
  FLOAT aa;
  FLOAT aspectRatio = 1.0;
  int dump = 0;

  int i,j;
  FILE *fp;

  snprintf (paramfile, MAX_FILENAME, "%s/%s", getenv("HOME"), ".fresnelparams");

  // read in command line arguments
  thisfile  =  argv[0];
  if ( argc < 2 || argc > 4 ) {
    printf("Usage: %s aa paramFile [dump]\n", thisfile);
    exit(EXIT_FAILURE);
  } else {
    aa = atof(argv[1]);
    if ( argc >= 3 )
      snprintf(paramfile, MAX_FILENAME, "%s", argv[2]);
    if ( argc == 4 )
      dump = atoi(argv[3]);
  }

  // open file pointer
  if ( (fp = fopen(paramfile,"r")) == NULL ) {
    perror("Opening file paramFile");
    exit(EXIT_FAILURE);
  }



  // read in parameters
  value_tmp = -1;
  while ( fgets(paramLine, MAXLINE, fp) != NULL ) {


    // get rid of any comment lines or blank lines
    if ( (paramLine[0]=='#') || ( strlen(paramLine)<2 ) )
      continue;

#if defined(USE_FLOAT)
#define FMT "%f %f %d %f %f %f %d %f %f %f %d"
#else
#define FMT "%le %le %d %le %le %le %d %le %le %le %d"
#endif

    if ( sscanf (paramLine, FMT, &lambLo, &lambHi, &nLambda, &maxX00, &x00Step, &RStarProj, &Nptsrc, &AU, &offset, &aspectRatio, &order) != 11 ) {
      perror("reading parameter file (parameters missing?)");
      exit(EXIT_FAILURE);
    }
  }

  //printf ("Lo%g Hi%g n%d max%g Step%g R%g N%d AU%g off%g asp%.2g ord%d\n", lambLo, lambHi, nLambda, maxX00, x00Step, RStarProj, Nptsrc, AU, offset, aspectRatio, order);

  // close file pointer
  if ( fclose(fp) < 0 ) {
    perror("close file");
    exit(EXIT_FAILURE);
  }

  /* add the stellar disk radius as a buffer for smoothing */
  maxX00 += RStarProj;

  /* ************************************************************
   *Divide up the ellipse into small boxes  
   *
   * *********************************************************** */

  FLOAT smaj = aa;
  FLOAT smin = aa / aspectRatio;
  FLOAT k_val = sqrt(smaj*smaj - smin*smin) / smaj;

  FLOAT segLength = ( 1.0/pow(2,order) ) * (1.0/smaj) * PIHALF * sqrt( (smaj*smaj + smin*smin)/2.0 );
  int Nsegment = 1;
  FLOAT arclength = 0;
  int pow2order = pow(2,order);
  FLOAT dtheta = PIHALF/200.0;

  struct box {
    FLOAT xur;
    FLOAT yur;
    FLOAT xll;
    FLOAT yll;
    int order;
  }  coords[pow2order];


  FLOAT r;

  int ord;
  int ordtmp;
  FLOAT theta;
  int mod;

  // integrate along the arc and staor the x,y values as they're encountered
  int nbox = 0;
  for (theta=0; theta<PIHALF; theta+=dtheta) {

    // how far we've integrated so far
    arclength += sqrt( 1.0 - k_val*k_val*sin(theta)*sin(theta) ) * dtheta;
    
    // if we've just gone one 2^n'th of the distance 0..PIHALF,
    //    then we're at the coords of one of the segments
    if (  ( ( arclength - Nsegment*segLength ) > 0.0 ) && Nsegment < pow2order ) {
      r = sqrt ( (smaj*smaj * smin*smin) / ( (smaj*smaj*sin(theta)*sin(theta)) + (smin*smin*cos(theta)*cos(theta)) ) );


      //  step through the orders and see which order this segment is
      //    ya, .. this works nicely.
      ord = 1;
      for (ordtmp=1; ordtmp<=order; ordtmp++) {
	ord = ordtmp;
	mod = (int) ( pow2order / pow(2,ord) );
	if ( (Nsegment % mod ) == 0)
	  break;
      }

      // put the upper-right and lower-left coords in the structure
      coords[nbox].order = ord;
      coords[nbox].xur = r*cos(theta);
      coords[nbox].yur = r*sin(theta);
      coords[nbox].xll = -r*cos(theta);
      coords[nbox].yll = -r*sin(theta);

      Nsegment++;
      nbox++;

    }

  }


  /*  Done divide of ellipse   ******************************** */




  /* ************************************************************
   *  Set the lower-left corner for each point so a small box is defined
   *
   * *********************************************************** */

  FLOAT xAct,yAct,xll,yll, xComp, yComp;
  int currentOrder;

  for (i=0; i<nbox; i++) {

    currentOrder = coords[i].order;
    xAct = coords[i].xur;
    yAct = coords[i].yur;
    xll = coords[i].xll;
    yll = coords[i].yll;

    for (j=0; j<nbox; j++) {
      
      if (coords[j].order >= currentOrder )
	continue;

      xComp = coords[j].xur;
      yComp = coords[j].yur;
      
      if ( xComp < xAct  && xComp > coords[i].xll )
	coords[i].xll = xComp;
      if ( yComp < yAct  && yComp > coords[i].yll )
	coords[i].yll = yComp;
      
    }

  }

  /*  Done lower-left x,y   ******************************** */



  /*  ************************************************************
   *   Mirror upper right quadrant to the other three  
   *      ... careful not to 'over-copy' boxes which straddle an axis
   *
   *  ************************************************************ */

  struct cbox {
    FLOAT x;
    FLOAT y;
    FLOAT dx;
    FLOAT dy;
  } cens[4*pow2order];

  FLOAT x,y,dx,dy;

  j = 0;
  for (i=0; i<nbox; i++) {
    
    // get centres and widths for this box

    x = (coords[i].xur+coords[i].xll)/2.0;
    y = (coords[i].yur+coords[i].yll)/2.0;
    dx = (coords[i].xur-coords[i].xll)/2.0;
    dy = (coords[i].yur-coords[i].yll)/2.0;
    cens[j].x = x;   cens[j].y = y;   cens[j].dx = dx;    cens[j++].dy = dy;

    // if the box in extirely in the 1st quadrant, copy to the other 3
    if (coords[i].xll>0 && coords[i].yll>0) {
      cens[j].x = -x;   cens[j].y = -y;   cens[j].dx = dx;  cens[j++].dy = dy;
      cens[j].x = x;    cens[j].y = -y;   cens[j].dx = dx;  cens[j++].dy = dy;
      cens[j].x = -x;   cens[j].y = y;    cens[j].dx = dx;  cens[j++].dy = dy;
      
    // if the bos straddles the xaxis in quads 1->4, copy to 2->3
    } else if (coords[i].xll>0 && coords[i].yll<0 ) {
      cens[j].x = -x;   cens[j].y = y;    cens[j].dx = dx;  cens[j++].dy = dy;
      
      // if the box straddles the yaxis in quads 1->2, copy to 3->4
    } else if (coords[i].xll<0 && coords[i].yll>0 ) {
      cens[j].x = x;    cens[j].y = -y;   cens[j].dx = dx;  cens[j++].dy = dy;
    }
    
  }
  nbox = j;


  /*  Done mirroring to other quadrants  ********************** */

    



  /* ************************************************************
   *  Build the background star
   *
   * *********************************************************** */

  // Ro is the radius of the background star in arcsec
  FLOAT Ro = RStarProj;

  int boxesPerSide = (int) sqrt(Nptsrc*4.0/PI);
  FLOAT boxWidth = 2.0 * Ro / ((FLOAT) boxesPerSide);

  // start just inside the circle
  // .. i did this becuase the point source is supposed to be at the center
  //    of the sub-box.  Otherwise, the center of a box would lie on the edge
  //    of the circle and be divided in half.
  FLOAT Rx = Ro - boxWidth/2.0;

  FLOAT xPt, yPt;

  // this is about 4/pi times bigger than necessary, but i was having a
  //   problem with overflows
  struct starCoords {
    FLOAT xPt;
    FLOAT yPt;
  } sources[boxesPerSide*boxesPerSide];

  int nStars=0;

  for (xPt=-Rx; xPt<=Ro; xPt+=boxWidth) {
    for (yPt=-Rx; yPt<=Ro; yPt+=boxWidth) {
      if ( (xPt*xPt + yPt*yPt) <= Ro*Ro ) {
	sources[nStars].xPt = xPt;
	sources[nStars].yPt = yPt;
	nStars++;
	//printf ("xPt: %f   yPt: %f\n",xPt,yPt);
      }
    }
  }

  /*  Done building background star   ******************************** */


  //for (i=0; i<nbox; i++) {
  //  printf ("%g %g %g %g\n", cens[i].x, cens[i].y, cens[i].dx, cens[i].dy );
  //}


  /* ************************************************************
   * Select lambda values
   *
   * *********************************************************** */

  FLOAT lambdaStep;

  if (nLambda>1) 
    lambdaStep = (lambHi - lambLo) / (nLambda -1);
  else 
    lambdaStep = (lambHi - lambLo);
  
  if (nLambda==1)
    lambLo = lambLo + 0.5*(lambHi - lambLo);

  FLOAT lambdas[MAX_NLAMBDA];

  for (i=0; i<nLambda; i++) {
    lambdas[i] = lambLo + i*lambdaStep;
    //printf ("%g %g %g\n",lambLo,lambHi,lambdas[i]);
  }

  /*  Done selecting lambda values   ******************************** */







  /* ************************************************************
   * compute the fresnel integrals 
   *
   * *********************************************************** */

  FLOAT C[MAX_FRES_INT+1];
  FLOAT S[MAX_FRES_INT+1];

  FLOAT intC, intS, Cint_prev, Sint_prev, Cprev, Sprev;
  FLOAT alphasqr;
  FLOAT dalpha = 0.001;


  //  FILE *fp_cornu;
  //if ( (fp_cornu = fopen("cornu_spiral", "w")) < 0 ) {
  //  perror ("opening output file: cornu_spiral");
  //  exit (EXIT_FAILURE);
  //} 
  
  Cprev = 0;
  Sprev = 0;
  for (j=0; j<MAX_FRES_INT; j++) {

    alphasqr = (FLOAT) j*j/1e6;
    intC = cos(PIHALF*alphasqr);
    intS = sin(PIHALF*alphasqr);
    
    Cint_prev = (j>0) ? C[j-1] : 0;
    Sint_prev = (j>0) ? S[j-1] : 0;

    // triangle rule integration
    C[j] = Cint_prev + 0.5*(Cprev+intC)*dalpha;
    S[j] = Sint_prev + 0.5*(Sprev+intS)*dalpha;
    //C[j] = Cint_prev + intC*dalpha;
    //S[j] = Sint_prev + intS*dalpha;

    Cprev = intC;
    Sprev = intS;

    //fprintf (fp_cornu, "%.6f %.6f %.6f\n", j/1000.0, C[j], S[j]);
  }
  C[MAX_FRES_INT] = 0.5;  // these are actually for j->infty
  S[MAX_FRES_INT] = 0.5;

  // if ( fclose(fp_cornu) < 0 ) {
  //  perror ("closing output file: cornu_spiral");
  // exit (EXIT_FAILURE);
  //}
  
  
  /*  Done computing fresnel integrals   ******************************** */









  /* ************************************************************
   * Build the diffraction patter for each box
   *
   * *********************************************************** */


  FLOAT z = AU * 1.5e11;
  FLOAT A = 0.5;
  FLOAT k;
  int Ncompon = nStars * nLambda;
  int ptSrcCount;
  FLOAT xStar, yStar;
  FLOAT lamb;

  FLOAT y00, y0, x00, x0, xCen, yCen;
  FLOAT xi0, eta0, eta1sign, eta2sign;
  FLOAT xi1sign = 1.0;
  FLOAT xi2sign = 1.0;
  FLOAT eta1f, eta2f, xi1f, xi2f;
  int   xi1, eta1, xi2, eta2;
  FLOAT Ceta, Seta, Cxi, Sxi;

  int nMaxX00 = maxX00 / x00Step + 2;
  int m,q;
  FLOAT Ur[nMaxX00], Ui[nMaxX00], Ihole[nMaxX00], Idisk[nMaxX00];
  FLOAT I0;

  int boxCount;
  int lambCount = 0;
  for (i=0; i<nLambda; i++) {
    
    lamb = lambdas[i];
    k = TWOPI / lamb;
    xi0 = sqrt (k / (PI*z));
    eta0 = xi0;
    y00 = offset;      
    xStar = yStar = 0.0;


    
    // 00 coords are inthe global system
    //  0 coords are wrt to the aperture.
    
    // see Goodman 'introduction to Fourier optics', pp. 70-74
    //    for explanation of diffraction calcs
    
    boxCount = 0;
    for (m=0; m<nbox; m++) {
      
      xCen = cens[m].x;       yCen = cens[m].y; 
      dx   = cens[m].dx; 	dy   = cens[m].dy;
      y0 = y00 + yStar - yCen;
      
      
      eta1f = -eta0 * (dy + y0);
      eta2f =  eta0 * (dy - y0); 
      eta1sign = (eta1f > 0.0) ? (1.0) : (-1.0); 
      eta2sign = (eta2f > 0.0) ? (1.0) : (-1.0); 
      eta1 = (int) fabs(1000*eta1f);
      eta2 = (int) fabs(1000*eta2f);
      
      // need to make sure we won't get a seg fault for large etaN
      //   Ceta and Seta are very small, so approximating to 0.5 (infty lim)
      //   prevents seg faults for rare values that go beyond the array.
      eta1 = (eta1<MAX_FRES_INT )? (eta1): (MAX_FRES_INT);
      eta2 = (eta2<MAX_FRES_INT )? (eta2): (MAX_FRES_INT);
      
      Ceta = eta2sign*C[eta2] - eta1sign*C[eta1];
      Seta = eta2sign*S[eta2] - eta1sign*S[eta1];
      
      
      /*  massive debugging
	  printf ("x: %.3f  y: %.3f  dx: %.3f  dy: %.3f  y0: %.3f  yStar: %.3f\n", xCen, yCen, dx, dy, y0, yStar);
	  printf ("\tz:  %.3e  AU:  %.3e   eta0: %.3e  xi0: %.3e\n", z, AU, eta0, xi0);
	  printf ("\tk:  %.3e  lamb:  %.3e  TWOPI: %.5e\n", k, lamb, TWOPI);
	  printf ("\tC2: %.3e  C1:  %.3e   S2: %.3e  S1:  %.3e\n", C[eta2], C[eta1], S[eta2], S[eta1]);
	  printf ("\teta1: %g  eta2: %g  ceta: %g  seta: %g\n", eta1f, eta2f, Ceta, Seta);	    
      */
      
      for (q=0; q<nMaxX00; q++) {
	
	//DB(1.0);
	x00 = q*x00Step;
	
	if ( boxCount == 0) {
	  Ur[q] = 0.0;
	  Ui[q] = 0.0;
	}
	
	x0 = x00 + xStar - xCen;
	
	xi1f = -xi0 * (dx + x0) * 1000;
	xi2f =  xi0 * (dx - x0) * 1000;
	xi1sign = (xi1f > 0) ? (1.0) : (-1.0);
	xi2sign = (xi2f > 0) ? (1.0) : (-1.0);
	xi1 = (int) fabs(xi1f);
	xi2 = (int) fabs(xi2f);
	
	// need to make sure we won't get a seg fault for large xiN
	xi1 = (xi1<MAX_FRES_INT) ? (xi1) : (MAX_FRES_INT);
	xi2 = (xi2<MAX_FRES_INT) ? (xi2) : (MAX_FRES_INT);
	Cxi = xi2sign*C[xi2] - xi1sign*C[xi1];
	Sxi = xi2sign*S[xi2] - xi1sign*S[xi1];
	
	// if these appear to be switched, have a look for a factor
	//  of -i in the constant coefficient ... you'll see it
	Ur[q] +=   A * (Cxi*Seta + Sxi*Ceta);
	Ui[q] +=   A * (Sxi*Seta - Cxi*Ceta);
	
	/*  Debugging
	    if ( x00 > 660 && x00 < 710) {
	    printf ("\t xi1: %g  xi2: %g  cxi: %g  sxi: %g\n", xi1f, xi2f, Cxi, Sxi);
	    }
	*/
	
	//DB(2.0);
      }
      boxCount++;
    }
    
    for (q=0; q<nMaxX00; q++) {
      
      if (! Ihole[q] )
	Ihole[q] = 0;
      if (! Idisk[q] )
	Idisk[q] = 0;
      
      I0 = ( Ur[q]*Ur[q] + Ui[q]*Ui[q] );
      Ihole[q] += I0;
      Idisk[q] += ( 1.0 - 2.0*Ur[q] + I0 ); 
      
    }    
    lambCount++;
  }
  
  /*  Done building diffraction pattern   ******************************** */
  
  /* smooth over the stellar disk */
  FLOAT IholeSum[nMaxX00], IdiskSum[nMaxX00];
  for (i=0;i<nMaxX00;i++) {
    IholeSum[i]=0.0;
    IdiskSum[i]=0.0;
  }
  maxX00 -= RStarProj;
  nMaxX00 =  maxX00 / x00Step + 1;
  ptSrcCount = 0;
  for (j=0; j<nStars; j++) {
    
    xStar = sources[j].xPt;
    yStar = sources[j].yPt;

    for (i=0; i<nMaxX00; i++) {

      /* get the distance and its index */
      FLOAT dx = xStar - (x00Step*i);
      FLOAT dy = yStar;
      FLOAT dist = sqrt( dx*dx + dy*dy );
      int i_dist = (int) (dist / x00Step);
      FLOAT frac = (dist - i_dist*x00Step)/x00Step;
      IholeSum[i] += Ihole[i_dist] + frac*(Ihole[i_dist] - Ihole[i_dist+1]);
      IdiskSum[i] += Idisk[i_dist] + frac*(Idisk[i_dist] - Idisk[i_dist+1]);
	

    }
    

  }
  ptSrcCount++;      


  
  /* ***********************************************************************
   *    OUTPUT 
   * 
   * ******************************************************************* */
  
  char outfile[MAX_FILENAME];

  int aaint = (int) aa;
  int AUint = (int) AU;
  snprintf (outfile, MAX_FILENAME, "fres-%05d_%05d", aaint, AUint);
  if ( (fp = fopen(outfile, "w")) < 0 ) {
    perror ("opening output file");
    exit (EXIT_FAILURE);
  }
  
  fprintf (fp, "# %s \tKBO_radius: \t%g\n", LAA, aa);
  fprintf (fp, "# %s \tlambda_Low: \t%g\n", LLAMBLO, lambLo);
  fprintf (fp, "# %s \tlambda_High: \t%g\n", LLAMBHI, lambHi);
  fprintf (fp, "# %s \tNum_lambdas: \t%d\n", LNLAMBDA, nLambda);
  fprintf (fp, "# %s \tmax_dist_from_shadow_centre(m): \t%g\n", LMAXX00, maxX00);
  fprintf (fp, "# %s \tstep_size: \t%g\n", LX00STEP, x00Step);
  fprintf (fp, "# %s \tProjected_Star_Radius: \t%g\n", LRSTAR, RStarProj);
  fprintf (fp, "# %s \tNum_point_sources_for_star: \t%d\n", LNSTARS, nStars);
  fprintf (fp, "# %s \tdist_to_KBO_(AU): \t%g\n", LAU, AU);
  fprintf (fp, "# %s \toffset_from_line_of_sight: \t%g\n", LOFFSET, offset);
  fprintf (fp, "# %s \taspectRatio_of_KBO_occulter(x/y): \t%.2g\n", LASPRAT, aspectRatio);
  fprintf (fp, "# %s order_of_approximation: \t%d\n", LORDER, order);

  for (q=0; q<nMaxX00; q++) {
    x00 = x00Step*q;
    fprintf(fp,"%.0f %.8f %.8f\n",x00,IholeSum[q]/Ncompon,IdiskSum[q]/Ncompon);
  }

  if ( fclose(fp) < 0 ) {
    perror ("closing output file");
    exit (EXIT_FAILURE);
  }


  // dump star coords and box coord info
  if ( dump ) {
    
    // stars
      snprintf (outfile, MAX_FILENAME, "fres-%05d_%05d.stars", aaint, AUint);
    if ( (fp = fopen(outfile, "w")) < 0 ) {
      perror ("opening star-dump output file");
      exit (EXIT_FAILURE);
    }
    fprintf (fp, "# %s \tProjected_Star_Radius: \t%g\n", LRSTAR, RStarProj);
    for (j=0; j<nStars; j++)
      fprintf (fp, "%.3f %.3f\n", sources[j].xPt, sources[j].yPt);
    
    if ( fclose(fp) < 0 ) {
      perror ("closing star-dump output file");
      exit (EXIT_FAILURE);
    }
    

    // boxes
    snprintf (outfile, MAX_FILENAME, "fres-%05d_%05d.boxes", aaint, AUint);
    if ( (fp = fopen(outfile, "w")) < 0 ) {
      perror ("opening box-dump output file");
      exit (EXIT_FAILURE);
    }
    for (m=0; m<nbox; m++)
      fprintf (fp, "%.3f %.3f %.3f %.3f\n", cens[m].x, cens[m].y, cens[m].dx, cens[m].dy);
    
    if ( fclose(fp) < 0 ) {
      perror ("closing box-dump output file");
      exit (EXIT_FAILURE);
    }
            

  }
  
  /* ***  Done Output  ********************************************** */
  
    
  return 0;
  
}
