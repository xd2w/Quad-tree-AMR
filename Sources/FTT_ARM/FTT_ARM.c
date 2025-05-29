#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "ftt.h"
#include "pfplib.h"
#include "nrutil.h"

// FTT: Fully Threaded Tree DATA Sturcture

int main(int argc, char *argv[])
{

  /* start pfplib for parameter readings */
  getPars(argc, argv);

  int it, itNb, tmax, iCell, ndata, tplot, flipped;
  FILE *fp;
  printf("--------------------------------\n");
  printf("Fully Threaded Tree Algorithms for Adaptive Refinement Methods Modified\n");
  printf("--------------------------------\n");

  dfetch("ellipse_A", &init_VOF_coefs[0]);
  dfetch("ellipse_B", &init_VOF_coefs[1]);
  dfetch("ellipse_C", &init_VOF_coefs[2]);

  fp = fopen("DATA/runtime.data", "w");
  if (fp == NULL)
  {
    printf("Error opening file DATA/runtime.data\n");
    exit(1);
  }

  initialize();
  // fttStatistics();
  // exit(1);
  tmax = 10;
  tplot = 1;
  ifetch("tmax", &tmax);
  ifetch("tplot", &tplot);
  //  plotFTT(0);
  //  plotSFC(0);
  initVOF(0);
  // plotFTTInterf(0);

  fprintf(fp, "Max Level = %d, tmax = %d, tplot = %d\n\n", maxLevel, tmax, tplot);
  fprintf(fp, "itNb, t_total, dt, #Cell\n");
  fprintf(fp, "--------------------------------\n");

  setPLICPramForAll();

  initTheoreticalInterf();

  t_total = 0;
  flipped = 0;

  for (itNb = 0; itNb <= tmax; itNb++)
  {
    printf("\nFully Threaded Tree Data Structure: iteration = %d\n", itNb);
    setTimeStep();

    if (t_total + global_dt >= 2 && t_total <= 2)
    {
      global_dt = 2 - t_total;
      // computeVelocityAtLeaves();
      // setPLICPramForAll();
      plic();
      setPLICPramForAll();
      advTheoreticalInterf();

      reMesh(itNb);
      setPLICPramForAll();
      t_total += global_dt;
      printf("t_total = %g\n", t_total);
      fprintf(fp, "%d, %f, %f, %d\n", itNb, t_total, global_dt, numberOfCells);
      itNb++;

      ndata = 200;
      plotFTT(ndata);
      plotSFC(ndata);
      plotHilbertSFC(ndata);
      plotFTTInterf(ndata);
      plotVOF(ndata);
      // printf("finished plotting vof\n");
      // plotCellGradAtIntf(ndata);
      plotCurvatureAtLeafCells(ndata);
      plotTheoreticalInterf(ndata);
    }
    if (t_total + global_dt >= 4 && t_total <= 4)
    {
      global_dt = 4 - t_total;
      // computeVelocityAtLeaves();
      // setPLICPramForAll();
      plic();
      setPLICPramForAll();
      advTheoreticalInterf();

      reMesh(itNb);
      setPLICPramForAll();
      t_total += global_dt;
      printf("t_total = %g\n", t_total);

      ndata = 400;
      plotFTT(ndata);
      plotSFC(ndata);
      plotHilbertSFC(ndata);
      plotFTTInterf(ndata);
      plotVOF(ndata);
      // printf("finished plotting vof\n");
      // plotCellGradAtIntf(ndata);
      plotCurvatureAtLeafCells(ndata);
      plotTheoreticalInterf(ndata);
      return 0;
    }

    printf("setTimeStep: completed successfully dt = %g\n", global_dt);

    // constructHilbertIndex();
    // printf("constructHilbertIndex: completed successfully \n");
    // fttStatistics();

    if (itNb % tplot == 0)
    {
      ndata = itNb / tplot;
      plotFTT(ndata);
      // plotSFC(ndata);
      // plotHilbertSFC(ndata);
      plotFTTInterf(ndata);
      // plotVOF(ndata);
      // printf("finished plotting vof\n");
      // plotCellGradAtIntf(ndata);
      plotCurvatureAtLeafCells(ndata);
      plotTheoreticalInterf(ndata);
      // fttStatistics();
    }
    // computeVelocityAtLeaves();
    // setPLICPramForAll();
    plic();
    setPLICPramForAll();
    advTheoreticalInterf();

    reMesh(itNb);
    setPLICPramForAll();
    t_total += global_dt;
    printf("t_total = %g\n", t_total);
    // fprintf(fp, "itNb, t_total, dt, #Cell\n");
    fprintf(fp, "%d, %f, %f, %d\n", itNb, t_total, global_dt, numberOfCells);
  }
  fclose(fp);
  return 0;
}
