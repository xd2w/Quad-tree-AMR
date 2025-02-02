#include <stdio.h>
#include <stdlib.h>
#include "ftt.h"
#include "nrutil.h"

void plotFTTOct(int iOct, FILE *fp)
{
  int iCell, cLv;
  Real xLeft, yLeft, xRight, yRight;
  Real dx, dy;

  iCell = cellNumberInOct*iOct;
  cLv = octLv[iOct];

  xLeft = xCell[iCell];
  yLeft = yCell[iCell];
  dx = dxCell[cLv];
  dy = dyCell[cLv];

  fprintf(fp,"%g %g\n", xLeft, yLeft+dy);
  fprintf(fp,"%g %g\n\n", xLeft+2.0*dx, yLeft+dy);
  fprintf(fp,"%g %g\n", xLeft+dx, yLeft);
  fprintf(fp,"%g %g\n\n", xLeft+dx, yLeft+2.0*dy);
}

void drawFTTOct(int iOct)
{
  int iCell, cLv;
  Real xLeft, yLeft, xRight, yRight;
  Real dx, dy;

  iCell = cellNumberInOct*iOct;
  cLv = octLv[iOct];

  xLeft = xCell[iCell];
  yLeft = yCell[iCell];
  dx = dxCell[cLv];
  dy = dyCell[cLv];

  printf("%g %g\n", xLeft, yLeft+dy);
  printf("%g %g\n\n", xLeft+2.0*dx, yLeft+dy);
  printf("%g %g\n", xLeft+dx, yLeft);
  printf("%g %g\n\n", xLeft+dx, yLeft+2.0*dy);
}

void plotFlagOct(int iOct, FILE *fp)
{
  int iCell, cLv;
  Real xLeft, yLeft, xRight, yRight;
  Real dx, dy;

  iCell = cellNumberInOct*iOct;
  cLv = octLv[iOct];

  xLeft = xCell[iCell];
  yLeft = yCell[iCell];
  dx = dxCell[cLv];
  dy = dyCell[cLv];

  fprintf(fp,"%g %g\n", xLeft, yLeft);
  fprintf(fp,"%g %g\n\n", xLeft+2.0*dx, yLeft+2.0*dy);
  fprintf(fp,"%g %g\n", xLeft+2.0*dx, yLeft);
  fprintf(fp,"%g %g\n\n", xLeft, yLeft+2.0*dy);
}

void plotOctMesh(int ndata)
{
  int iOct, i, i1, i2, i3;
  char fname[] = "DATA/Mesh.000";
  FILE *fp;
  i = ndata;
  i1 = i % 10; i /= 10;
  i2 = i % 10; i /= 10;
  i3 = i % 10;
  fname[10] ='0'+i3;
  fname[11] ='0'+i2;
  fname[12] ='0'+i1;
  fp =  fopen(fname, "w");
  fprintf(fp,"%g %g\n", -Lx, -Ly);
  fprintf(fp,"%g %g\n", 3*Lx, -Ly);
  fprintf(fp,"%g %g\n", 3*Lx, 3*Ly);
  fprintf(fp,"%g %g\n", -Lx, 3*Ly);
  fprintf(fp,"%g %g\n", -Lx, -Ly);

  for(iOct=0; iOct<numberOfOcts; iOct++)
  {
    plotFTTOct(iOct, fp);
  }
  fclose(fp);
}
void plotFlagOcts(int ndata)
{
  int iOct, i, i1, i2, i3;
  char fname[] = "DATA/fOct.000";
  FILE *fp;
  i = ndata;
  i1 = i % 10; i /= 10;
  i2 = i % 10; i /= 10;
  i3 = i % 10;
  fname[10] ='0'+i3;
  fname[11] ='0'+i2;
  fname[12] ='0'+i1;
  fp =  fopen(fname, "w");
  for(iOct=0; iOct<numberOfOcts; iOct++)
  {
    if(octFlag[iOct]) plotFlagOct(iOct, fp);
  }
}
void plotFlagOctsAtLevel(int ndata, int level)
{
  int iOct, i, i1, i2, i3;
  char fname[] = "DATA/fOct.000";
  FILE *fp;
  i = ndata;
  i1 = i % 10; i /= 10;
  i2 = i % 10; i /= 10;
  i3 = i % 10;
  fname[10] ='0'+i3;
  fname[11] ='0'+i2;
  fname[12] ='0'+i1;
  fp =  fopen(fname, "w");
  for(iOct=0; iOct<numberOfOcts; iOct++)
  {
    if(octFlag[iOct] && octLv[iOct]==level) plotFlagOct(iOct, fp);
  }
}

void plotNgFlagOcts(int ndata)
{
  int iOct, i, i1, i2, i3;
  char fname[] = "DATA/nOct.000";
  FILE *fp;
  i = ndata;
  i1 = i % 10; i /= 10;
  i2 = i % 10; i /= 10;
  i3 = i % 10;
  fname[10] ='0'+i3;
  fname[11] ='0'+i2;
  fname[12] ='0'+i1;
  fp =  fopen(fname, "w");
  for(iOct=0; iOct<numberOfOcts; iOct++)
  {
    if(octFlag[iOct]==0 && octLv[iOct] == maxLevel) plotFlagOct(iOct, fp);
  }
}
void plotNgFlagOctsAtLevel(int ndata, int level)
{
  int iOct, i, i1, i2, i3;
  char fname[] = "DATA/nOct.000";
  FILE *fp;
  i = ndata;
  i1 = i % 10; i /= 10;
  i2 = i % 10; i /= 10;
  i3 = i % 10;
  fname[10] ='0'+i3;
  fname[11] ='0'+i2;
  fname[12] ='0'+i1;
  fp =  fopen(fname, "w");
  for(iOct=0; iOct<numberOfOcts; iOct++)
  {
    if(octFlag[iOct]==0 && octLv[iOct] == level) plotFlagOct(iOct, fp);
  }
}

