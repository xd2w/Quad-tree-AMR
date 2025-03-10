#include <stdio.h>
#include <stdlib.h>
#include "ftt.h"
#include "nrutil.h"
/* seek for the leave cell in which point (x, y) is */
/* start from a cell in which the point is and the leave
  cell is the smallest child */
int seekCell(int iCell, Real x, Real y)
{
  int iOct, chCell, within, cLv;
  Real dx, dy;
  iOct = cellChOct[iCell];

  printf("%g %g is in cell %d type %d\n", x, y, iCell, cellType[iCell]);
  if (iOct == 0)
  {
    printf("*****************************\n");
    printf("%g %g is in cell %d\n", x, y, iCell);
    printf("*****************************\n");
    printFTTCell(iCell);
    return iCell;
  }
  cLv = octLv[iOct];
  dx = dxCell[cLv];
  dy = dyCell[cLv];

  chCell = iOct * 4;
  if (xCell[chCell] <= x && x < xCell[chCell] + dx &&
      yCell[chCell] <= y && y < yCell[chCell] + dy)
  {
    seekCell(chCell, x, y);
  }
  chCell++;
  if (xCell[chCell] <= x && x < xCell[chCell] + dx &&
      yCell[chCell] <= y && y < yCell[chCell] + dy)
  {
    seekCell(chCell, x, y);
  }
  chCell++;
  if (xCell[chCell] <= x && x < xCell[chCell] + dx &&
      yCell[chCell] <= y && y < yCell[chCell] + dy)
  {
    seekCell(chCell, x, y);
  }
  chCell++;
  if (xCell[chCell] <= x && x < xCell[chCell] + dx &&
      yCell[chCell] <= y && y < yCell[chCell] + dy)
  {
    seekCell(chCell, x, y);
  }
  // return -1;
}
