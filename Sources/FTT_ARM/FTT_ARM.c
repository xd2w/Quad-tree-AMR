#include <stdio.h>
#include <stdlib.h>
#include "ftt.h"
#include "pfplib.h"

// FTT: Fully Threaded Tree DATA Sturcture


int main(int argc, char *argv[])
{

/* start pfplib for parameter readings */
  getPars(argc, argv);


  int it, itNb, tmax, iCell, ndata, tplot;
  FILE* fp;
  printf("--------------------------------\n");
  printf("Fully Threaded Tree Algorithms for Adaptive Refinement Methods Modified\n");
  printf("--------------------------------\n");

  dfetch("ellipse_A", &init_VOF_coefs[0]);
  dfetch("ellipse_B", &init_VOF_coefs[1]);
  dfetch("ellipse_C", &init_VOF_coefs[2]);


  // init_VOF_coefs[0] = 1 ; // a
  // init_VOF_coefs[1] = 0.1; // b
  // init_VOF_coefs[2] = 2; // c

  // a*X^2 + b*X*Y + c*Y^2 = radius^2

  initialize();
  fttStatistics();
  // exit(1);
  tmax = 10;
  tplot = 1;
  ifetch("tmax", &tmax);
  ifetch("tplot", &tplot);
//  plotFTT(0);
//  plotSFC(0);
  initVOF(0);
  // plotFTTInterf(0); 
  //seekCell(0,.343, .473); printf("vof %g\n", vof[55]); exit(1);
  // exit(1);
  for(itNb=0; itNb<=tmax;itNb++)
  {
    printf("\nFully Threaded Tree Data Structure: iteration = %d\n", itNb);

    if(itNb % tplot == 0)
    {
      ndata = itNb/tplot;
      plotFTT(ndata);
      plotSFC(ndata);
      plotHilbertSFC(ndata); 
      plotFTTInterf(ndata);
    }
    plic();
    reMesh(itNb);
    fttStatistics();
  }  
  return 1;
}
