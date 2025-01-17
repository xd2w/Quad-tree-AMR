#include <stdio.h>
#include <stdlib.h>
#include "ftt.h"
#include "nrutil.h"

void markVOF(void)
{
  int iCell, iLv, iOct, iFlag, chCell;
  Real fraction;
  const Real MinVof = 1.0e-12, MaxVof = 1.0 - MinVof;

  for(iCell=0; iCell<numberOfCells; iCell++)
  {
    cellMark[iCell] = cellFlag[iCell];
    fraction = vof[iCell];
    cellFlag[iCell] =  0;
    if(fraction > MinVof && fraction < MaxVof) cellFlag[iCell] = 1;
  }
  for(iCell=0; iCell<numberOfCells; iCell++)
  {
    if(cellMark[iCell] != cellFlag[iCell])
    {
      printf("error in markVOF\n");
      printf("cell %d vof %g\n", iCell, vof[iCell]);
      printf("cellMark %d cellFlag %d\n", cellMark[iCell], cellFlag[iCell]);
      printFTTCell(iCell); 
      exit(1);
    }
  }
}
