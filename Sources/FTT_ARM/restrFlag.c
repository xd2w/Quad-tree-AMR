#include <stdio.h>
#include <stdlib.h>
#include "ftt.h"
#include "nrutil.h"

void restrFlag(void)
{
  int it;

 // return;
  for(it=maxLevel-1; it>=0; it--)
  {
    restrFlagAtLevel(it);
  }  
}

void restrFlagAtLevel(int lev)
{
  int iCell, iLv, iOct, iFlag, chCell;

  for(iCell=0; iCell<numberOfCells; iCell++)
  {
    iOct = iCell/cellNumberInOct;
    iLv = octLv[iOct];
    chCell = 4*cellChOct[iCell];
    if(iLv == lev && chCell != 0)
    { 
      if(cellFlag[chCell]) cellFlag[iCell] = 1;
      chCell++;
      if(cellFlag[chCell]) cellFlag[iCell] = 1;
      chCell++;
      if(cellFlag[chCell]) cellFlag[iCell] = 1;
      chCell++;
      if(cellFlag[chCell]) cellFlag[iCell] = 1;
    }
  }
  printf("level %d flag %d\n", lev, cellFlag[0]);
}
