#include <stdio.h>
#include <stdlib.h>
#include "ftt.h"
#include "nrutil.h"

void drawPrCells(int iCell)
{
  int iOct, prCell;
  FILE* fp;
  fp = fopen("prCells","w");
  printFTTCell(iCell);
  plotFTTCell(iCell, fp);

  iOct = iCell/4;
  while(iOct > 0)
  { 
    prCell = octPrCell[iOct];
    printf("parent cell %d\n", prCell);
    plotFTTCell(prCell, fp);
    iOct = prCell/4;
  }
  fclose(fp);
}
