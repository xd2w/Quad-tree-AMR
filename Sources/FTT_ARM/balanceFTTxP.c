#include <stdio.h>
#include <stdlib.h>
#include "ftt.h"
#include "nrutil.h"
#include "box.h"

int localIndex[50];
/* check the level of x positive-direction neigbour. If the neigbour is
at level  less than levelDrop, split the neigbour. */
void balanceFTTxP(int tmpCell, int levelDrop)
{
  int it;
  int tmpOct, tmpPrCell, iLv, tmpLv, tmpIndex;
  int iCell = tmpCell;
  tmpOct = tmpCell/4;
  tmpPrCell = octPrCell[tmpOct];
  tmpLv = octLv[tmpOct];
  tmpIndex = tmpCell%4; 
/* on south-west or north-west cell */
  if(tmpIndex % 2 == 0)
  {
    tmpCell++;
    /* tmpCell: east nieghbour.  */
    if(levelDrop == 0 && cellChOct[tmpCell] == 0) splitCell(tmpCell);
    return;
  }
/* on south-east or north-east cell */
  iLv = tmpLv;
  localIndex[iLv] = tmpIndex;
//  iLv --;
/* clim up until reaching west-south or west-north  cell*/
  while(iLv >0 && (tmpIndex % 2 == 1))
  {
    tmpCell = octPrCell[tmpOct];
    tmpOct = tmpCell/4;
    iLv = octLv[tmpOct];
    localIndex[iLv] = tmpIndex = tmpCell%4;
  }
  if(tmpIndex % 2 == 1) return; // right boundary cell, do nothing
  
  if( localIndex[iLv]== 0) /* south-west cell */
  {
    tmpCell = tmpOct*4+1; /* south-est cell */
  }
  else
  /* localIndex[iLv]== 2, north-west cell */
  {
    tmpCell = tmpOct*4+3; /* north-est cell */
  }
  if(cellChOct[tmpCell] == 0) splitCell(tmpCell);
  tmpOct = cellChOct[tmpCell];
  //printFTTCell(tmpCell); 
  while(iLv < tmpLv-levelDrop)
  {  
    tmpOct = cellChOct[tmpCell];
    iLv ++;
    if(localIndex[iLv]==1) /* south-est cell */
    {
      tmpCell = tmpOct*4; /* south-west cell */
    }
    else /*localIndex[iLv]=3, north-est cell */
    {
      tmpCell = tmpOct*4+2; /*north-west cell */
    }
    if(cellChOct[tmpCell] == 0) splitCell(tmpCell);
  }
}
