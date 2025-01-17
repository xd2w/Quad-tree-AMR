#include <stdio.h>
#include <stdlib.h>
#include "ftt.h"
#include "nrutil.h"

void drawNgbCells(int iCell)
{
  int ngbCell, southCell, northCell;
  FILE* fp;
  fp = fopen("ngbCells","w");

  printf("draw cell %d neighbour cells\n", iCell);

  ngbCell = cellNb[0][iCell];
  if(ngbCell >= 0) plotFTTCell(ngbCell, fp);

  southCell = cellNb[2][ngbCell];
  if(southCell >= 0) plotFTTCell(southCell, fp);

  northCell = cellNb[3][ngbCell];
  if(northCell >= 0) plotFTTCell(northCell, fp);

  ngbCell = cellNb[1][iCell];
  if(ngbCell >= 0) plotFTTCell(ngbCell, fp);

  southCell = cellNb[2][ngbCell];
  if(southCell >= 0) plotFTTCell(southCell, fp);

  northCell = cellNb[3][ngbCell];
  if(northCell >= 0) plotFTTCell(northCell, fp);


  ngbCell = cellNb[2][iCell];
  if(ngbCell >= 0) plotFTTCell(ngbCell, fp);

  ngbCell = cellNb[3][iCell];
  if(ngbCell >= 0) plotFTTCell(ngbCell, fp);

  fclose(fp);


}

