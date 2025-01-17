#include <stdio.h>
#include <stdlib.h>
#include "ftt.h"
#include "nrutil.h"

void drawChCells(int iCell, FILE* fp)
{
  int iOct, chCell;

  plotFTTCell(iCell, fp);
  iOct = cellChOct[iCell];
  if(iOct>0)
  {
    chCell = iOct*4;
//    printf("iCell %d (%d) chCell %d (%d) %d (%d) %d (%d) %d (%d)\n", 
          iCell, cellFlag[iCell], 
          chCell,cellFlag[chCell], chCell+1,cellFlag[chCell+1], 
//          chCell+2,cellFlag[chCell+2], chCell+3, cellFlag[chCell+3]);
    drawChCells(chCell, fp);
    chCell += 1;
    drawChCells(chCell, fp);
    chCell += 1;
    drawChCells(chCell, fp);
    chCell += 1;
    drawChCells(chCell, fp);
  } 

}


