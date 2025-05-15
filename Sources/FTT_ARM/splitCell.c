#include <stdio.h>
#include <stdlib.h>
#include <math.h>
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

  Real list[4], checkSum;
  getChildVOF(0, list, iCell);
  checkSum = 0;

  vof[chCell] = list[0];
  checkSum += vof[chCell];
  vof[chCell + 1] = list[1];
  checkSum += vof[chCell + 1];
  vof[chCell + 2] = list[2];
  checkSum += vof[chCell + 2];
  vof[chCell + 3] = list[3];
  checkSum += vof[chCell + 3];

  if (fabs(checkSum - 4 * cellVOF) > 1e-5)
  {
    printf("*******************************************\n");
    printf("splitCell.c/splitCell_smart: vof miss-match\n\n");
    exit(1);
  }

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

void splitCell_smart(int iCell)
{
  int cLv, iOct, chCell, cellTp;
  Real xLeft, yLeft;
  Real dx2, dy2, cellVOF;
  Real list[4], checkSum, cc[3][3];

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
  // vof[chCell] = cellVOF;
  // vof[chCell + 1] = cellVOF;
  // vof[chCell + 2] = cellVOF;
  // vof[chCell + 3] = cellVOF;
  cellFlag[chCell] = cellFlag[iCell];
  cellFlag[chCell + 1] = cellFlag[iCell];
  cellFlag[chCell + 2] = cellFlag[iCell];
  cellFlag[chCell + 3] = cellFlag[iCell];

  getChildVOF(0, list, iCell);
  checkSum = 0;

  vof[chCell] = list[0];
  checkSum += vof[chCell];
  vof[chCell + 1] = list[1];
  checkSum += vof[chCell + 1];
  vof[chCell + 2] = list[2];
  checkSum += vof[chCell + 2];
  vof[chCell + 3] = list[3];
  checkSum += vof[chCell + 3];

  getCellNgbVOF(chCell, cc);
  setPLICPramForOne(chCell, cc);
  getCellNgbVOF(chCell + 1, cc);
  setPLICPramForOne(chCell + 1, cc);
  getCellNgbVOF(chCell + 2, cc);
  setPLICPramForOne(chCell + 2, cc);
  getCellNgbVOF(chCell + 3, cc);
  setPLICPramForOne(chCell + 3, cc);

  if (fabs(checkSum - 4 * cellVOF) > 1e-5)
  {
    printf("*******************************************\n");
    printf("splitCell.c/splitCell_smart: vof miss-match\n\n");
    printf("original vof = %f\n", vof[iCell]);
    printf("split into :\n");
    printf("\t%f\t%f\n", list[2], list[3]);
    printf("\t%f\t%f\n", list[0], list[1]);
    exit(1);
  }

  numberOfOcts++;
  if (numberOfOcts > maxNumberOfOcts)
  {
    printf("numberOfOcts > maxNumberOfOcts in splitCell\n");
    printf("allocated memory exceeded\n");
    exit(1);
  }
  numberOfCells = numberOfOcts * cellNumberInOct;

  // printf("cell[%d] split\n", iCell);

  // establish nb + balance recursively
  int i, dir, dest, nbChOct, prNbCell, prprNbCell, leftCell, rightCell;
  for (i = 0; i < 4; i++)
  {
    for (dir = 0; dir < 4; dir++)
    {
      dest = morton_lookup[dir][i];
      if (dest < 4)
      { // siblings
        cellNb[dir][chCell + i] = chCell + dest;
      }
      else
      { // outer neighbours
        prNbCell = cellNb[dir][iCell];

        if (octLv[prNbCell / 4] == octLv[iCell / 4])
        {
          if (cellChOct[prNbCell])
          {
            nbChOct = cellChOct[prNbCell / 4];
            cellNb[dir][chCell + i] = 4 * nbChOct + (dest - 4);
          }
          else
          {
            cellNb[dir][chCell + i] = prNbCell;
          }
        }
        else
        {
          // printf("error unbalanced \n");
          // exit(1);
          if (cellChOct[prNbCell] == 0)
          {
            // printf("recur split\n");
            // exit(0);
            splitCell_smart(prNbCell);
            nbChOct = cellChOct[prNbCell];
            cellNb[dir][chCell + i] = 4 * nbChOct + (dest - 4);
          }
          else
          {
            nbChOct = cellChOct[prNbCell / 4];
            cellNb[dir][chCell + i] = 4 * nbChOct + (dest - 4);
            // printf("neighbour relation error\n");
            // exit(1);
          }
        }
      }
    }
    // leftCell = cellNb[0][iCell];
    // rightCell = cellNb[1][iCell];

    // if (octLv[leftCell / 4] != cLv)
    // {
    //   printf("error");
    //   exit(1);
    // }

    // if (octLv[rightCell / 4] != cLv)
    // {
    //   printf("error\n");
    //   exit(1);
    // }

    // prNbCell = cellNb[3][leftCell];
    // if (octLv[prNbCell / 4] > cLv)
    //   splitCell_smart(prNbCell);

    // prNbCell = cellNb[2][leftCell];
    // if (octLv[prNbCell / 4] > cLv)
    //   splitCell_smart(prNbCell);

    // prNbCell = cellNb[3][rightCell];
    // if (octLv[prNbCell / 4] > cLv)
    //   splitCell_smart(prNbCell);

    // prNbCell = cellNb[2][rightCell];
    // if (octLv[prNbCell / 4] > cLv)
    //   splitCell_smart(prNbCell);
  }

  return;
}
