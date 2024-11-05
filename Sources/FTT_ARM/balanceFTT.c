#include <stdio.h>
#include <stdlib.h>
#include "ftt.h"
#include "pfplib.h"

void balanceFTT(int balanceLevel, int levelDrop)
{
  int iCell, iLv, iOct, iFlag, cellNumber;
  int iBalance;
 
  iBalance=0;
  ifetch("iBalance", &iBalance);
  if(iBalance == 0) return;
  
  cellNumber = numberOfCells;
  //printf("\n\n*** balanceFTT at level %d *** \n\n", balanceLevel);
  for(iCell=0; iCell<cellNumber; iCell++)
  {
    cellFlag[iCell] = 0;
    if(cellChOct[iCell] && cellType[iCell] == 0 ) cellFlag[iCell] = 1;
  }
  for(iCell=0; iCell<cellNumber; iCell++)
  {
    iOct = iCell/cellNumberInOct;
    iLv = octLv[iOct];
/* cellChOct[iCell] > 0 => iCell is not a leaf */
    if(balanceLevel == iLv && cellFlag[iCell])
    {
//        printf("balanceLevel %d iLv %d iCell %d\n", 
//                balanceLevel, iLv, iCell);
//        printFTTCell(iCell);
      balanceFTTxP(iCell, levelDrop);
      balanceFTTyN(iCell, levelDrop);
      balanceFTTxN(iCell, levelDrop);
      balanceFTTyP(iCell, levelDrop);
    } 
  }
 // printf("end of balance tree at level %d\n", balanceLevel);
}
