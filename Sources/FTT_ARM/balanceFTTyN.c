#include <stdio.h>
#include <stdlib.h>
#include "ftt.h"
#include "nrutil.h"
#include "box.h"

extern int localIndex[50];
void balanceFTTyN(int tmpCell, int levelDrop)
{
  int it;
  int tmpOct, tmpPrCell, iLv, tmpLv, tmpIndex;
 
  tmpOct = tmpCell/4;
  tmpPrCell = octPrCell[tmpOct];
  tmpLv = octLv[tmpOct];
  tmpIndex = tmpCell%4;
  if(tmpIndex > 1)
  {
    tmpCell -= 2;
    if(levelDrop == 0 && cellChOct[tmpCell] == 0) splitCell(tmpCell);
    return;
  }
  iLv = tmpLv;
  localIndex[iLv] = tmpIndex;
//  iLv --;
  while(iLv >0 && (tmpIndex < 2))
  {
    tmpCell = octPrCell[tmpOct];
    tmpOct = tmpCell/4;
    iLv = octLv[tmpOct];
    localIndex[iLv] = tmpIndex = tmpCell%4;
  }
  if(tmpIndex < 2) return; // bottom boundary cell, do nothing
  if( localIndex[iLv]== 2)
  {
    tmpCell = tmpOct*4;
  }
  else
  {
    tmpCell = tmpOct*4+1;
  }
  if(cellChOct[tmpCell] == 0) splitCell(tmpCell);
  //printFTTCell(tmpCell);  

while(iLv < tmpLv-levelDrop)
{  
  tmpOct = cellChOct[tmpCell];
  iLv ++;
  if(localIndex[iLv]==0)
  {
    tmpCell = tmpOct*4+2;
  }
  else
  {
    tmpCell = tmpOct*4+3;
  }
  if(cellChOct[tmpCell] == 0) splitCell(tmpCell);
}
}
