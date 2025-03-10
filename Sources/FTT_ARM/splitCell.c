#include <stdio.h>
#include <stdlib.h>
#include "ftt.h"
/* creat a new Oct: 4 children cells and their SW corner coordinates */
/* establish relation of parent cell (iCell) and child Oct (numberOfOcts) */

void splitCell(int iCell)
{
  int cLv, iOct, chCell, cellTp;
  Real xLeft, yLeft;
  Real dx2, dy2, cellVOF;

  /* number of Oct: cell number / 4 */
  iOct = iCell / cellNumberInOct;
  cLv = octLv[iOct]; /* level of 4 cells in an oct is the same */
                     //  printf("split cell #%d at level %d\n", iCell, cLv);
                     //  printFTTCell(iCell);
  /* establish parent Cell <--> child Oct */
  cellChOct[iCell] = numberOfOcts; // iCell's child Oct is numberOfOcts
  octPrCell[numberOfOcts] = iCell; // Oct numberOfOcts's parent cell is iCell
  octLv[numberOfOcts] = cLv + 1;   /* increase the oct level by 1 */

  /* define children cells */
  // chCell = 0: SW; 1: SE; 2: NW; 3:NE
  dx2 = 0.5 * dxCell[cLv];
  dy2 = 0.5 * dyCell[cLv];
  xLeft = xCell[iCell];
  yLeft = yCell[iCell];
  chCell = numberOfCells;
  cellTp = cellType[iCell];
  /* children cells have no child Oct */
  cellChOct[chCell] = 0;
  cellChOct[chCell + 1] = 0;
  cellChOct[chCell + 2] = 0;
  cellChOct[chCell + 3] = 0;
  cellType[chCell] = cellTp;
  cellType[chCell + 1] = cellTp;
  cellType[chCell + 2] = cellTp;
  cellType[chCell + 3] = cellTp;
  xCell[chCell] = xLeft;
  yCell[chCell] = yLeft;
  xCell[chCell + 1] = xLeft + dx2;
  yCell[chCell + 1] = yLeft;
  xCell[chCell + 2] = xLeft;
  yCell[chCell + 2] = yLeft + dy2;
  xCell[chCell + 3] = xLeft + dx2;
  yCell[chCell + 3] = yLeft + dy2;
  cellVOF = vof[iCell];
  // printf("cellVOF %g\n", cellVOF); exit(1); // vof not initiated
  vof[chCell] = cellVOF;
  vof[chCell + 1] = cellVOF;
  vof[chCell + 2] = cellVOF;
  vof[chCell + 3] = cellVOF;

  // Real list[4];
  // getChildVOF(0, list, iCell);
  // vof[chCell] = list[0];
  // vof[chCell + 1] = list[1];
  // vof[chCell + 2] = list[2];
  // vof[chCell + 3] = list[3];

  numberOfOcts++;
  if (numberOfOcts > maxNumberOfOcts)
  {
    printf("numberOfOcts > maxNumberOfOcts in splitCell\n");
    printf("allocated memory exceeded\n");
    exit(1);
  }
  numberOfCells = numberOfOcts * cellNumberInOct;

  return;
}
