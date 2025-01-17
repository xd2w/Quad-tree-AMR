#include <stdio.h>
#include <stdlib.h>
#include "ftt.h"
#include "nrutil.h"

/* 
Check the neighbourhood among cells.
0 -> west neighbour, 1 -> est neighbour
2 -> south neighbour, 3 -> neighbour
Negative number neighbour means no neighbour.
*/

void bcOctTree(void)
{
// west boundary
  bcOctTreeAtCell(6, 1);
  bcOctTreeAtCell(12, 1);
// est boundary
  bcOctTreeAtCell(10, 0);
  bcOctTreeAtCell(16, 0);
// south boundary
  bcOctTreeAtCell(5, 3);
// north boundary
  bcOctTreeAtCell(15, 2);
// south-west corner
  bcOctTreeAtCell(4, 1);
// north-west corner
  bcOctTreeAtCell(14, 1);
// south-east corner
  bcOctTreeAtCell(8, 0);
// north-east corner
  bcOctTreeAtCell(18, 0);

}

void bcOctTreeAtCell(int iCell, int direction)
{ 
  int iOct, chOct, chCell, octCell, ngbOct, ngbCell, ngbLv, iLv;
  static int dirShift1[4] = {0, 1, 0, 2};
  static int dirShift2[4] = {2, 3, 1, 3};
  chOct = cellChOct[iCell];
  if(chOct)
  { 
    octCell = 4*chOct;
    chCell = octCell+dirShift1[direction];
    bcOctTreeAtCell(chCell, direction);
    chCell = octCell+dirShift2[direction];;
    bcOctTreeAtCell(chCell,direction);
  }
  else
  { 
    iOct = iCell/4;
    ngbCell = cellNb[direction][iCell];
    ngbOct = ngbCell/4;
    iLv = octLv[iOct];
    ngbLv = octLv[ngbOct];
    if(iLv > ngbLv)
    { 
      printf("mistake in neighbourhood at direction %d\n", direction); exit(1);
    }
    // neighbour cell has child 
    if(cellChOct[ngbCell])
    { 
      splitCell(iCell);
    }
  }
}

