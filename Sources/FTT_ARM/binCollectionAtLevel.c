#include <stdio.h>
#include <stdlib.h>
#include "ftt.h"
#include "nrutil.h"

/* NOTE BOOK (VI), page 12-13 */
/* delect the oct with octlFlag = 0 at level "level"
and replaced by the last cell with octFlag !=0 at level "level"
or not at level "level"*/
/* update the parent-child relationship */

// does garbage collection of octs at a given level
// gets rid of holes in the array of octs (and consiquently cells)
void binCollectionAtLevel(int level)
{
  int id, iCell, iLv, iOct, iFlag, octCell, chOct, child, prCell, octNb;
  FILE *fp;

  octNb = numberOfOcts;
  // printf("before binCollection at level %d: numberOfOcts = %d numberOfCells %d\n", level, numberOfOcts,numberOfCells);
  for (iOct = 1; iOct < numberOfOcts; iOct++)
  {
    octCell = 4 * iOct;
    if (octFlag[iOct] == 0 && octLv[iOct] == level)
    {
      //      printf("iOct = %d level = %d\n", iOct, octLv[iOct]);
      numberOfOcts--;
      while (octFlag[numberOfOcts] == 0 && octLv[numberOfOcts] == level)
      {
        prCell = octPrCell[numberOfOcts];
        cellChOct[prCell] = 0;
        // printf("numberOfOcts %d parent Cell %d\n", numberOfOcts, prCell);
        numberOfOcts--;
      }
      iCell = numberOfOcts * 4;
      // have found the last oct with octFlag !=0
      //      printf("iOct %d octCell %d iCell %d\n", iOct, octCell, iCell);
      if (iOct < numberOfOcts)
      {
        //        printf("numberOfOcts = %d level = %d\n", numberOfOcts, octLv[numberOfOcts]);
        prCell = octPrCell[numberOfOcts];
        // printf("numberOfOcts %d parent Cell %d\n", numberOfOcts, prCell);
        cellChOct[prCell] = iOct;
        prCell = octPrCell[iOct];
        cellChOct[prCell] = 0;
        octLv[iOct] = octLv[numberOfOcts];
        octPrCell[iOct] = octPrCell[numberOfOcts];

        // filling hole created by unrefining with cells at the back
        // octCell -> Cell of IOct and iCell -> Cell of numberOfOcts
        // swapping all the data of those 2 cell

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
        vof[octCell] = vof[iCell];
        vof[octCell + 1] = vof[iCell + 1];
        vof[octCell + 2] = vof[iCell + 2];
        vof[octCell + 3] = vof[iCell + 3];

        cellFlag[octCell] = cellFlag[iCell];
        cellFlag[octCell + 1] = cellFlag[iCell + 1];
        cellFlag[octCell + 2] = cellFlag[iCell + 2];
        cellFlag[octCell + 3] = cellFlag[iCell + 3];
        octFlag[iOct] = octFlag[numberOfOcts];

        cellType[octCell] = cellType[iCell];
        cellType[octCell + 1] = cellType[iCell + 1];
        cellType[octCell + 2] = cellType[iCell + 2];
        cellType[octCell + 3] = cellType[iCell + 3];

        // turn off the child oct pointer of the deleted Oct
        //        prCell = octPrCell[numberOfOcts];
        //        cellChOct[prCell] = 0;
      }
    }
  }
  if (iOct > numberOfOcts)
    numberOfOcts++;
  numberOfCells = numberOfOcts * cellNumberInOct;

  /* If a child is flaged, so is its parent. It is not true
  other way around. Delete children connection if children is not flaged */
  /*
    for(iCell=0; iCell<numberOfCells; iCell++)
    {
      iOct=cellChOct[iCell];
      octCell = 4*iOct;
      if(cellFlag[octCell] == 0) cellChOct[iCell] = 0;
    }

    for(iCell=0; iCell<numberOfCells; iCell++)
    {
      iOct = iCell/4;
      iLv = octLv[iOct];
    }
  */

  // printf("after binCollection at level %d: numberOfOcts = %d numberOfCells %d\n",level, numberOfOcts, numberOfCells);
  // printf("%d oct deleted\n", octNb-numberOfOcts);
}
