#include <stdio.h>
#include <stdlib.h>
#include "ftt.h"
#include "nrutil.h"
#include "box.h"
//extern int oldCell;

/* check the level of x negative-direction neigbour. If the neigbour is
at level less than levelDrop, split the neigbour. */
void balanceFTTxN(int tmpCell, int levelDrop)
{
  int it;
  int tmpOct, tmpPrCell, iLv, tmpLv, tmpIndex;
  int localIndex[50];
  tmpOct = tmpCell/4;
  tmpPrCell = octPrCell[tmpOct];
  tmpLv = octLv[tmpOct];
  tmpIndex = tmpCell%4;
  if(tmpIndex % 2 == 1)
  {
    tmpCell--;
    if(levelDrop == 0 && cellChOct[tmpCell] == 0) splitCell(tmpCell);
    return;
  }
  iLv = tmpLv;
  localIndex[iLv] = tmpIndex;
//  iLv --;
  while(iLv >0 && (tmpIndex % 2 == 0))
  {
    tmpCell = octPrCell[tmpOct];
    tmpOct = tmpCell/4;
    iLv = octLv[tmpOct];
    localIndex[iLv] = tmpIndex = tmpCell%4;
  }
  if(tmpIndex % 2 == 0) return; // left boundary cell, do nothing
  if( localIndex[iLv]== 1)
  {
    tmpCell = tmpOct*4;
  }
  else
  {
    tmpCell = tmpOct*4+2;
  }
  if(cellChOct[tmpCell] == 0) splitCell(tmpCell);
  //printFTTCell(tmpCell);  

while(iLv < tmpLv-levelDrop)
{  
  tmpOct = cellChOct[tmpCell];
  iLv ++;
  if(localIndex[iLv])
  {
    tmpCell = tmpOct*4+3;
  }
  else
  {
    tmpCell = tmpOct*4+1;
  }
  if(cellChOct[tmpCell] == 0) splitCell(tmpCell);
}
}
