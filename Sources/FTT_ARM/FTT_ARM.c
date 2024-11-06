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
  FILE *fp;
  printf("--------------------------------\n");
  printf("Fully Threaded Tree Algorithms for Adaptive Refinement Methods\n");
  printf("--------------------------------\n");

  initialize();
  fttStatistics();
  tmax = 10;
  tplot = 1;
  ifetch("tmax", &tmax);
  ifetch("tplot", &tplot);
  plotFTT(0);
  plotSFC(0);
  // initVOF(0);
  plotFTTInterf(0);
  // seekCell(0,.343, .473); printf("vof %g\n", vof[55]); exit(1);

  exit(1);

  for (itNb = 0; itNb <= tmax; itNb++)
  {
    printf("\nFully Threaded Tree Data Structure: iteration = %d\n", itNb);

    if (itNb % tplot == 0)
    {
      // all below for plotting
      ndata = itNb / tplot;
      plotFTT(ndata);
      plotSFC(ndata);
      plotHilbertSFC(ndata);
      plotFTTInterf(ndata);
    }
    plic();          // calculating VOF + Flagging cells
    reMesh(itNb);    // re organises the mesh
    fttStatistics(); // counts the number of cells and octs
  }
  return 1;
}
