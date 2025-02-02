#include <stdio.h>
#include <stdlib.h>
#include "ftt.h"
#include "nrutil.h"

void flagOct(void)
{
  int iCell, iLv, iOct, iFlag, octCell;

  for(iCell=0; iCell<numberOfCells; iCell++)
  {
    iOct = iCell/cellNumberInOct;
    octCell = 4*iOct;
    if(cellFlag[iCell])
    {
      cellFlag[octCell  ] = 1;
      cellFlag[octCell+1] = 1;
      cellFlag[octCell+2] = 1;
      cellFlag[octCell+3] = 1;
    }
  }
}
