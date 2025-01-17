#include <stdio.h>
#include <stdlib.h>
#include "ftt.h"
#include "nrutil.h"

/* 
Elstablish the neighbourhood among cells.
0 -> west neighbour, 1 -> east neighbour
2 -> south neighbour, 3 -> north neighbour
Negative number neighbour means no neighbour.
*/
/*
A cell is allowed to have neighbour at lower level.
If a cell's neighbour is at a level lower, then
the neighbour's neighbour is not the cell itself,
but cell's ancestor at the same leve as the neighbour.
*/
/*
In a balanced tree, a cell's neighbour can be at maximun
one level lower.
*/
void establishNb(void)
{
  int iLv;

  establishNbAtBase();
  for(iLv=1; iLv <= maxLevel; iLv++)
  {
    establishNbAtLevel(iLv);
  }
}
 
void establishNbAtBase(void)
{
  int iCell, iLv;

  for(iCell=0; iCell<maxNumberOfCells; iCell++)
  {
    cellNb[0][iCell] = maxNumberOfCells;
    cellNb[1][iCell] = maxNumberOfCells;
    cellNb[2][iCell] = maxNumberOfCells;
    cellNb[3][iCell] = maxNumberOfCells;
  }
/* initial the neighbourhood data at the first oct, at level 0 */
// cell 0
  cellNb[0][0] = -1;
  cellNb[1][0] =  1;
  cellNb[2][0] = -1;
  cellNb[3][0] =  2;
// cell 1
  cellNb[0][1] =  0;
  cellNb[1][1] = -1;
  cellNb[2][1] = -1;
  cellNb[3][1] =  3;
// cell 2
  cellNb[0][2] = -1;
  cellNb[1][2] =  3;
  cellNb[2][2] =  0;
  cellNb[3][2] = -1;
// cell 3
  cellNb[0][3] =  2; 
  cellNb[1][3] = -1; 
  cellNb[2][3] =  1; 
  cellNb[3][3] = -1; 
  return;
}

/* Neighbourhood has been established at nghLv-1 */
void establishNbAtLevel(int nghLv)
{
  int iLv, iOct, octCell, prCell, nghOct, nghPrCell, nghCell;

  for(iOct=1; iOct<numberOfOcts; iOct++)
  {
    iLv = octLv[iOct];
    if(iLv == nghLv)
    {
      octCell = 4*iOct; 
      prCell = octPrCell[iOct];
// sibbling neighbourhood
      cellNb[0][octCell+1] = octCell  ;
      cellNb[0][octCell+3] = octCell+2;
      cellNb[1][octCell  ] = octCell+1;
      cellNb[1][octCell+2] = octCell+3;
      cellNb[2][octCell+2] = octCell  ;
      cellNb[2][octCell+3] = octCell+1;
      cellNb[3][octCell  ] = octCell+2;
      cellNb[3][octCell+1] = octCell+3;
// west neighbour
      nghPrCell = cellNb[0][prCell];
      if(nghPrCell>=0) // parent  has neighbour
      {
        nghOct = cellChOct[nghPrCell];
        if(nghOct) // parent neighbour has children
        {
          nghCell= 4*nghOct;
          cellNb[0][octCell  ] = nghCell+1;
          cellNb[0][octCell+2] = nghCell+3;
        }
        else // parent neighbour has no child
        {
          cellNb[0][octCell  ] = nghPrCell;
          cellNb[0][octCell+2] = nghPrCell;
        } 
      }
      else // parent  has no neighbour
      {
        cellNb[0][octCell  ] = -1;
        cellNb[0][octCell+2] = -1;
      } 
// east neighbour 
      nghPrCell = cellNb[1][prCell];
      if(nghPrCell>=0) // parent  has a neighbour
      {
        nghOct = cellChOct[nghPrCell];
        if(nghOct) //  parent neighbour has children
        {
          nghCell= 4*nghOct;
          cellNb[1][octCell+1] = nghCell;
          cellNb[1][octCell+3] = nghCell+2;
        }
        else // parent neighbour has no child
        {
          cellNb[1][octCell+1] = nghPrCell;
          cellNb[1][octCell+3] = nghPrCell;
        }
      }
      else // parent  has no neighbour
      {
        cellNb[1][octCell+1] = -1;
        cellNb[1][octCell+3] = -1;
      }
// south neighbour
      nghPrCell = cellNb[2][prCell];
      if(nghPrCell>=0)
      {
        nghOct = cellChOct[nghPrCell];
        if(nghOct)
        {
          nghCell= 4*nghOct;
          cellNb[2][octCell  ] = nghCell+2;
          cellNb[2][octCell+1] = nghCell+3;
        }
        else
        {
          cellNb[2][octCell  ] = nghPrCell;
          cellNb[2][octCell+1] = nghPrCell;
        }
      }
      else
      {
        cellNb[2][octCell  ] = -1;
        cellNb[2][octCell+1] = -1;
      }
// north neighbour
      nghPrCell = cellNb[3][prCell];
      if(nghPrCell>=0)
      {
        nghOct = cellChOct[nghPrCell];
        if(nghOct)
        {
          nghCell= 4*nghOct;
          cellNb[3][octCell+2] = nghCell  ;
          cellNb[3][octCell+3] = nghCell+1;
        }
        else
        {
          cellNb[3][octCell+2] = nghPrCell;
          cellNb[3][octCell+3] = nghPrCell;
        }
      }
      else
      {
        cellNb[3][octCell+2] = -1;
        cellNb[3][octCell+3] = -1;
      }
    }
  }
} 
