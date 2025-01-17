#include <stdio.h>
#include <stdlib.h>
#include "ftt.h"
#include "nrutil.h"
#include "box.h"

extern int localIndex[50];
void balanceFTTyP(int tmpCell, int levelDrop)
{
  int it;
  int tmpOct, tmpPrCell, iLv, tmpLv, tmpIndex;

  tmpOct = tmpCell/4;
  tmpPrCell = octPrCell[tmpOct];
  tmpLv = octLv[tmpOct];
  tmpIndex = tmpCell%4;
//  printFTTCell(tmpPrCell);
  //printf("tmpOct %d tmpLv %d\n", tmpOct, tmpLv);
  //printf("tmpCell %d local index %d\n\n", tmpCell, tmpIndex);
  if(tmpIndex < 2)
  {
    tmpCell += 2;
    if(levelDrop == 0 && cellChOct[tmpCell] == 0) splitCell(tmpCell);
    return;
  }

  iLv = tmpLv;
  localIndex[iLv] = tmpIndex;
 // iLv --;
  while(iLv >0 && (tmpIndex > 1))
  {
    tmpCell = octPrCell[tmpOct];
    tmpOct = tmpCell/4;
    iLv = octLv[tmpOct];
    localIndex[iLv] = tmpIndex = tmpCell%4;
    //printf("tmpOct %d iLv %d\n", tmpOct, iLv);
    //printf("tmpCell %d local index %d\n\n", tmpCell, tmpIndex);
  }
  //printf("localIndex = %d at level %d\n", localIndex[iLv], iLv); 
  //printFTTCell(tmpCell);  
  if(tmpIndex > 1) return; // top boundary cell, do nothing
  if( localIndex[iLv]== 0)
  {
    tmpCell = tmpOct*4+2;
  }
  else
  {
    tmpCell = tmpOct*4+3;
  }
  if(cellChOct[tmpCell] == 0) splitCell(tmpCell);
  //printFTTCell(tmpCell);  

while(iLv < tmpLv-levelDrop)
{  
  //printf("iLv %d tmpLv = %d\n", iLv,  tmpLv);
  tmpOct = cellChOct[tmpCell];
  iLv ++;
  //printf("localIndex %d at level %d\n", localIndex[iLv], iLv);
  if(localIndex[iLv]==2)
  {
    tmpCell = tmpOct*4;
  }
  else
  {
    tmpCell = tmpOct*4+1;
  }
  if(cellChOct[tmpCell] == 0) splitCell(tmpCell);
}
}
