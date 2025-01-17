#include <stdio.h>
#include <stdlib.h>
#include "ftt.h"
#include "nrutil.h"

/* NOTE BOOK (VI), page 12-13 */
/* delect the cell with cellFlag = 0 and replaced by the last cell with cellFlag !=0 */
/* update the parent-child relationship */

void garbageCollection(void)
{
  int id, iCell, iLv, iOct, iFlag, octCell, chOct, child, prCell;
  FILE *fp;

  for (iOct = 1; iOct < numberOfOcts; iOct++)
  {
    if (cellChOct[0] == 0)
    {
      printf("cell 0, cellChOct %d\n", cellChOct[0]);
      printf("iOct = %d\n", iOct);
      exit(1);
    }
    octCell = 4 * iOct;
    if (cellFlag[octCell] == 0)
    {
      numberOfOcts--;
      iCell = numberOfOcts * 4;
      while (cellFlag[iCell] == 0)
      {
        numberOfOcts--;
        iCell = numberOfOcts * 4;
      }
      // have found the last oct with octFlag !=0
      if (iOct < numberOfOcts)
      {
        prCell = octPrCell[numberOfOcts];
        cellChOct[prCell] = iOct;
        prCell = octPrCell[iOct];
        cellChOct[prCell] = 0;
        octLv[iOct] = octLv[numberOfOcts];
        octPrCell[iOct] = octPrCell[numberOfOcts];

        child = cellChOct[octCell];
        chOct = cellChOct[octCell] = cellChOct[iCell];
        octPrCell[chOct] = octCell;
        octPrCell[child] = iCell;
        child = cellChOct[octCell + 1];
        chOct = cellChOct[octCell + 1] = cellChOct[iCell + 1];
        octPrCell[chOct] = octCell + 1;
        octPrCell[child] = iCell + 1;
        child = cellChOct[octCell + 2];
        chOct = cellChOct[octCell + 2] = cellChOct[iCell + 2];
        octPrCell[chOct] = octCell + 2;
        octPrCell[child] = iCell + 2;
        child = cellChOct[octCell + 3];
        chOct = cellChOct[octCell + 3] = cellChOct[iCell + 3];
        octPrCell[chOct] = octCell + 3;
        octPrCell[child] = iCell + 3;
        // copy coordinates
        xCell[octCell] = xCell[iCell];
        yCell[octCell] = yCell[iCell];
        xCell[octCell + 1] = xCell[iCell + 1];
        yCell[octCell + 1] = yCell[iCell + 1];
        xCell[octCell + 2] = xCell[iCell + 2];
        yCell[octCell + 2] = yCell[iCell + 2];
        xCell[octCell + 3] = xCell[iCell + 3];
        yCell[octCell + 3] = yCell[iCell + 3];
        // copy vof
        // vof[octCell  ] = vof[iCell  ];
        // vof[octCell+1] = vof[iCell+1];
        // vof[octCell+2] = vof[iCell+2];
        // vof[octCell+3] = vof[iCell+3];

        vx[octCell + 0] = vx[iCell + 0];
        vx[octCell + 1] = vx[iCell + 1];
        vx[octCell + 2] = vx[iCell + 2];
        vx[octCell + 3] = vx[iCell + 3];

        vy[octCell + 0] = vy[iCell + 0];
        vy[octCell + 1] = vy[iCell + 1];
        vy[octCell + 2] = vy[iCell + 2];
        vy[octCell + 3] = vy[iCell + 3];

        cellFlag[octCell] = 1;
        cellFlag[octCell + 1] = 1;
        cellFlag[octCell + 2] = 1;
        cellFlag[octCell + 3] = 1;
        cellFlag[iCell] = 0;
        cellFlag[iCell + 1] = 0;
        cellFlag[iCell + 2] = 0;
        cellFlag[iCell + 3] = 0;
      }
    }
  }
  if (iOct > numberOfOcts)
    numberOfOcts++;
  numberOfCells = numberOfOcts * cellNumberInOct;

  /* If a child is flaged, so is its parent. It is not true
  other way around. Delete children connection if children is not flaged */
  for (iCell = 0; iCell < numberOfCells; iCell++)
  {
    iOct = cellChOct[iCell];
    octCell = 4 * iOct;
    if (cellFlag[octCell] == 0)
      cellChOct[iCell] = 0;
  }

  for (iCell = 0; iCell < numberOfCells; iCell++)
  {
    iOct = iCell / 4;
    iLv = octLv[iOct];
  }
}
