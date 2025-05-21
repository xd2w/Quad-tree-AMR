#include <stdio.h>
#include <stdlib.h>
#include "ftt.h"
#include "nrutil.h"

void plotFTTCell(int iCell, FILE *fp)
{
  int iOct, cLv;
  Real xLeft, yLeft, xRight, yRight;

  // do not draw if cell is not a leaf
  //  if(cellFlag[iCell]) return;

  iOct = iCell / cellNumberInOct;
  cLv = octLv[iOct];
  // printf("iCell %d iOct %d cLv %d\n", iCell, iOct, cLv);
  xLeft = xCell[iCell];
  yLeft = yCell[iCell];
  xRight = xLeft + dxCell[cLv];
  yRight = yLeft + dyCell[cLv];
  // printf("x %g y %g dx %g dy %g\n", xLeft, yLeft, dxCell[cLv], dyCell[cLv]);
  fprintf(fp, "%g %g\n", xLeft, yLeft);
  fprintf(fp, "%g %g\n", xRight, yLeft);
  fprintf(fp, "%g %g\n", xRight, yRight);
  fprintf(fp, "%g %g\n", xLeft, yRight);
  fprintf(fp, "%g %g\n\n", xLeft, yLeft);

  return;
}

void plotFTT(int ndata)
{
  int iCell, i, i1, i2, i3;
  char fname[] = "DATA/mesh.000";
  char fsfcv[] = "DATA/sfcv.000";
  char fpGradName[] = "DATA/fgrd.000";
  FILE *fp, *fsfci, *fpGrad;
  i = ndata;
  i1 = i % 10;
  i /= 10;
  i2 = i % 10;
  i /= 10;
  i3 = i % 10;
  fname[10] = '0' + i3;
  fname[11] = '0' + i2;
  fname[12] = '0' + i1;

  fpGradName[10] = '0' + i3;
  fpGradName[11] = '0' + i2;
  fpGradName[12] = '0' + i1;

  fp = fopen(fname, "w");
  fpGrad = fopen(fpGradName, "w");
  for (iCell = 0; iCell < numberOfCells; iCell++)
  {
    if (cellChOct[iCell] == 0)
    {
      plotFTTCell(iCell, fp);
      plotCellGrad_4x4(iCell, fpGrad);
      // plotCellGradSmoothed_4x4(iCell, fpGrad);
    }
    //    if(cellChOct[iCell]==0 && vof[iCell]>0.99999) plotFlagCell(iCell, fp);
  }
  fclose(fp);
  fclose(fpGrad);
  printf("plot ftt mesh: numberOfCells %d\n", numberOfCells);
  // exit(1);
  return;
}

void plotFTTAtLevel(int ndata, int level)
{
  int iCell, i, i1, i2, i3;
  char fname[] = "DATA/mesh_lev0.000";
  FILE *fp;
  i = ndata;
  i1 = i % 10;
  i /= 10;
  i2 = i % 10;
  i /= 10;
  i3 = i % 10;
  fname[13] = '0' + level;
  fname[15] = '0' + i3;
  fname[16] = '0' + i2;
  fname[17] = '0' + i1;

  fp = fopen(fname, "w");
  for (iCell = 0; iCell < numberOfCells; iCell++)
  {
    if (octLv[iCell / 4] == level)
    {
      plotFTTCell(iCell, fp);
    }
    //    if(cellChOct[iCell]==0 && vof[iCell]>0.99999) plotFlagCell(iCell, fp);
  }
  fclose(fp);
  //  printf("plot ftt mesh: numberOfCells %d\n",numberOfCells);
  return;
}

void plotFlagFTT(int ndata)
{
  int iCell, i, i1, i2, i3;
  char fname[] = "DATA/fCell.000";
  FILE *fp;
  i = ndata;
  i1 = i % 10;
  i /= 10;
  i2 = i % 10;
  i /= 10;
  i3 = i % 10;
  fname[11] = '0' + i3;
  fname[12] = '0' + i2;
  fname[13] = '0' + i1;
  fp = fopen(fname, "w");

  for (iCell = 0; iCell < numberOfCells; iCell++)
  {
    if (cellFlag[iCell] && cellChOct[iCell] == 0)
    {
      plotFlagCell(iCell, fp);
    }
  }

  fclose(fp);
}
void plotFlagCellsAtLevel(int ndata, int level)
{
  int iCell, iOct, i, i1, i2, i3;
  char fname[] = "DATA/fCell.000";
  FILE *fp;
  i = ndata;
  i1 = i % 10;
  i /= 10;
  i2 = i % 10;
  i /= 10;
  i3 = i % 10;
  fname[11] = '0' + i3;
  fname[12] = '0' + i2;
  fname[13] = '0' + i1;
  fp = fopen(fname, "w");

  for (iCell = 0; iCell < numberOfCells; iCell++)
  {
    iOct = iCell / cellNumberInOct;
    if (cellFlag[iCell] && octLv[iOct] == level)
    {
      plotFlagCell(iCell, fp);
    }
  }

  fclose(fp);
}

void plotFlagCell(int iCell, FILE *fp)
{
  int iOct, cLv;
  Real xLeft, yLeft, xRight, yRight;

  // do not draw if cell is not a leaf
  //  if(cellFlag[iCell]) return;

  iOct = iCell / cellNumberInOct;
  cLv = octLv[iOct];
  /*   if(cLv != maxLevel)
     {
       printf("plotFlagCell: iCell %d iOct %d cLv %d\n", iCell, iOct, cLv);
     }
  */
  xLeft = xCell[iCell];
  yLeft = yCell[iCell];
  xRight = xLeft + dxCell[cLv];
  yRight = yLeft + dyCell[cLv];
  // printf("x %g y %g dx %g dy %g\n", xLeft, yLeft, dxCell[cLv], dyCell[cLv]);
  fprintf(fp, "%g %g\n", xLeft, yLeft);
  fprintf(fp, "%g %g\n\n", xRight, yRight);
  fprintf(fp, "%g %g\n", xRight, yLeft);
  fprintf(fp, "%g %g\n\n", xLeft, yRight);
}
