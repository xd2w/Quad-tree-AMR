#include <stdio.h>
#include <stdlib.h>
#include "ftt.h"
#include "nrutil.h"
#include "box.h"

void initVOF(int itNb)
{
  int iCell, iLv, iOct, iFlag, cellNumber;
 
// flag the leave cells
  cellNumber = numberOfCells;
  for(iCell=0; iCell<cellNumber; iCell++)
  {
    cellFlag[iCell] = 0;
  }
  for(iCell=0; iCell<cellNumber; iCell++)
  {
    if(cellChOct[iCell]==0)
    {
      vof[iCell] = computeVOF(iCell, itNb);
    }
  }
  restrField(vof);

  return;
}
