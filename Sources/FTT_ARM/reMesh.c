#include <stdio.h>
#include <stdlib.h>
#include "ftt.h"
#include <math.h>
#include "pfplib.h"

extern void splitFlagCellsAtLevel(int level);
extern void oct_PrCellFlagAtLvel(int level);
extern void refineToKappa(void);
extern void refineToKappaAtLevel(int level);
extern void cleanToKappaAtLevel(int level);

void reMesh(int itNb)
{
  int level;

  setCellInt1DZero(cellFlag);
  setOctInt1DZero(octFlag);
  flagInterfCells();
  // refineToKappa();
  // flagInterfLeaves();
  // printf("flagInterfLeavs done successfully \n");
  // // exit(0);

  // for (level = maxLevel; level > maxLevel - 1; level--)
  for (level = maxLevel; level > minLevel; level--)
  {
    // cleanToKappaAtLevel(level);
    setOctInt1DZeroAtLevel(octFlag, level);
    cell_OctFlagAtLevel(level);
    setCellInt1DZeroAtLevel(cellFlag, level);
    propagateOctFlagAtLevel(level);
    propagateOctFlagAtLevel(level);
    // refineToKappaAtLevel(level);
    // refineToKappaAtLevel(level - 1);
    splitFlagCellsAtLevel(level - 1);
    binCollectionAtLevel(level);
    oct_PrCellFlagAtLvel(level);
    establishNb();
    checkNb();
  }
  restrField(vof);

  // plotFTTInterf(200);
  // exit(0);

  /*
    drawNgbCells(510);
    printCellNgbVOF(510);
    level=maxLevel-1;
    setOctInt1DZeroAtLevel(octFlag, level);
    cell_OctFlagAtLevel(level);
    plotFlagOctsAtLevel(0, level);
      propagateOctFlagAtLevel(level);
    plotFlagOctsAtLevel(19, level);
      plotFlagCellsAtLevel(19, level-1);
      propagateOctFlagAtLevel(level);
      splitFlagCellsAtLevel(level-1);
    plotFlagOctsAtLevel(20, level);
      plotFlagCellsAtLevel(22, level-1);
      plotNgFlagOctsAtLevel(20, level);
      binCollectionAtLevel(level);
      plotFlagCellsAtLevel(23, level-1);
  */
}

void refineToKappa(void)
{
  int iCell, iOct, iLv, cOct, i, j, dir, dest, prNbCell;
  Real fraction, kappa, cc[6][6], list[4], ccp[3][3], checkSum, refine_th;
  int currentNumberOfCells = numberOfCells;
  dfetch("refine_threshold", &refine_th);
  for (iCell = 0; iCell < currentNumberOfCells; iCell++)
  {
    cellFlag[iCell] = 0;
    if (cellChOct[iCell] == 0)
    {
      fraction = vof[iCell];
      if (fraction > 0.0 && fraction < 1.0)
      {
        iOct = iCell / cellNumberInOct;
        iLv = octLv[iOct];
        getCellNgbVOF_6x6(iOct, cc);
        kappa = kappaBarickALELike(iCell, cc);
        if (log(kappa + 1) > refine_th * (iLv))
        {
          balanceCellsAround(iCell);
          splitCell_smart(iCell);

          // getCellNgbVOF(iCell, ccp);
          // for (i = 0; i < 3; i++)
          // {
          //   for (j = 0; j < 3; j++)
          //   {
          //     printf("\t%f,", ccp[j][i]);
          //   }
          //   printf("\n");
          // }
          // printf("\n");

          // getChildVOF(0, list, iCell);
          // checkSum = 0;
          // for (i = 0; i < cellNumberInOct; i++)
          // {
          //   vof[4 * cOct + i] = list[i];
          //   checkSum += list[i];
          //   for (dir = 0; dir < 4; dir++)
          //   {
          //     dest = morton_lookup[dir][4 * cOct + i];
          //     if (dest < 4)
          //     {
          //       cellNb[dir][4 * cOct + i] = 4 * cOct + dest;
          //     }
          //     else
          //     {
          //       prNbCell = cellNb[dir][iCell];
          //       if (octLv[prNbCell / 4] > octLv[iCell / 4])
          //       {
          //         // to make later
          //         cellNb[dir][4 * cOct + i] = prNbCell;
          //       }
          //       cellNb[dir][4 * cOct + i] = prNbCell;
          //     }
          //   }
          // }

          // // printf("sum = %f\n", 0.25 * checkSum);
          // if (fabs(0.25 * checkSum - vof[iCell]) > 1e-4)
          // {
          //   printf("*************\n");
          //   printf("vof missmatch\n\n");
          //   exit(1);
          // }
        }
      }
    }
  }
}

void balanceCell(int iCell, int iLv)
{
  if (octLv[iCell / 4] < iLv && cellChOct[iCell] == 0)
  {
    splitCell_smart(iCell);
    // splitCell(iCell);
  }
}

void refineToKappaAtLevel(int level)
{
  int iCell, iOct, iLv, cOct, i, j, dir, dest, prNbCell;
  Real fraction, kappa, cc[6][6], list[4], ccp[3][3], checkSum;
  int currentNumberOfCells = numberOfCells;
  Real tol = 1e-2, refine_th;
  int chCell, ngbOctCell;

  dfetch("refine_threshold", &refine_th);
  for (iCell = 0; iCell < currentNumberOfCells; iCell++)
  {
    // cellFlag[iCell] = 0;
    if (cellChOct[iCell] == 0 && octLv[iCell / 4] == level)
    {
      fraction = vof[iCell];
      if (fraction > tol && fraction < 1.0 - tol)
      {
        iOct = iCell / cellNumberInOct;
        iLv = octLv[iOct];
        cellFlag[iCell] = 1;

        getCellNgbVOF_6x6(iOct, cc);
        kappa = kappaBarickALELike(iCell, cc);
        // kappa = kappaHF(iCell, cc);
        // if (-log(kappa) < log(8 * dxCell[iLv]))
        if (log(kappa + 1) > refine_th * (iLv))
        {
          balanceCellsAround(iCell);
          splitCell_smart(iCell);

          // chCell = 4 * cellChOct[iCell];

          // // west cell
          // ngbOctCell = cellNb[0][iCell];
          // balanceCell(ngbOctCell, iLv + 1);
          // // south-west cell
          // balanceCell(cellNb[2][ngbOctCell], iLv + 1);
          // // north-west cell
          // balanceCell(cellNb[3][ngbOctCell], iLv + 1);

          // // east cell
          // ngbOctCell = cellNb[1][iCell];
          // balanceCell(ngbOctCell, iLv + 1);
          // // south-east cell
          // balanceCell(cellNb[2][ngbOctCell], iLv + 1);
          // // north-east cell
          // balanceCell(cellNb[3][ngbOctCell], iLv + 1);

          // // south cell
          // balanceCell(cellNb[2][iCell], iLv + 1);

          // // north cell
          // balanceCell(cellNb[3][iCell], iLv + 1);
        }
      }
    }
  }
}

void cleanToKappaAtLevel(int level)
{
  int iCell, iOct, iLv, cOct, i, j, dir, dest, prNbCell;
  Real fraction, kappa, cc[6][6], list[4], ccp[3][3], checkSum;
  int currentNumberOfCells = numberOfCells;
  Real tol = 1e-2, refine_th;
  int chCell, ngbOctCell;

  dfetch("refine_threshold", &refine_th);
  for (iCell = 0; iCell < currentNumberOfCells; iCell++)
  {
    if (cellChOct[iCell] == 0 && octLv[iCell / 4] == level && cellFlag[iCell] == 1)
    {
      fraction = vof[iCell];
      if (fraction > tol && fraction < 1.0 - tol)
      {
        iOct = iCell / cellNumberInOct;
        iLv = octLv[iOct];
        cellFlag[iCell] = 1;

        getCellNgbVOF_6x6(iOct, cc);
        kappa = kappaBarickALELike(iCell, cc);
        // kappa = kappaHF(iCell, cc);
        // if (-log(kappa) > -log(8 * dxCell[iLv - 1]))
        if (log(kappa) < refine_th * (iLv - 2))
        {
          cellFlag[iCell] = 0;
        }
      }
    }
  }
}

void balanceCellsAround(int iCell)
{
  int ngbOctCell, prCell;
  Real list[4];
  int iLv;

  prCell = iCell;
  iLv = octLv[iCell / 4];

  // west cell
  ngbOctCell = cellNb[0][prCell];
  balanceCell(ngbOctCell, iLv);
  // south-west cell
  balanceCell(cellNb[2][ngbOctCell], iLv);
  // north-west cell
  balanceCell(cellNb[3][ngbOctCell], iLv);

  // east cell
  ngbOctCell = cellNb[1][prCell];
  balanceCell(ngbOctCell, iLv);
  // south-east cell
  balanceCell(cellNb[2][ngbOctCell], iLv);
  // north-east cell
  balanceCell(cellNb[3][ngbOctCell], iLv);

  // south cell
  balanceCell(cellNb[2][prCell], iLv);

  // north cell
  balanceCell(cellNb[3][prCell], iLv);
}

// void refine_next_streaming_cells(int iCell)
// {
//   int iOct, prCell, ngbOctCell;
//   Real list[4];
//   int iLv;

//   prCell = iCell;
//   iLv = octLv[iCell / 4];

//   // west cell
//   ngbOctCell = cellNb[0][prCell];
//   balanceCell(ngbOctCell, iLv);

//   // east cell
//   ngbOctCell = cellNb[1][prCell];
//   balanceCell(ngbOctCell, iLv);

//   // south cell
//   ngbOctCell = cellNb[2][prCell];
//   balanceCell(ngbOctCell, iLv);

//   // north cell
//   ngbOctCell = cellNb[3][prCell];
//   balanceCell(ngbOctCell, iLv);
// }

// void _balanceCellsAround(int iOct)
// { // generates 6x6 grid of data from neighbours of octs
//   int ngbOct, ngbOctCell, southOct, southOctCell, northOct, northOctCell, prCell;
//   Real list[4];
//   int iLv;

//   prCell = octPrCell[iOct];
//   iLv = octLv[iOct];

//   // west cell
//   ngbOctCell = cellNb[0][prCell];
//   balanceCell(ngbOctCell, iLv);

//   if (octLv[ngbOctCell / 4] == iLv - 1)
//   {
//     if (morton_lookup[2][ngbOctCell % 4] > 3)
//     {
//       // south-west cell
//       southOctCell = cellNb[2][ngbOctCell];
//       balanceCell(southOctCell, iLv);
//     }
//     if (morton_lookup[3][ngbOctCell % 4] > 3)
//     {
//       // south-west cell
//       northOctCell = cellNb[3][ngbOctCell];
//       balanceCell(southOctCell, iLv);
//     }
//   }
//   else
//   {
//     // south-west cell
//     southOctCell = cellNb[2][ngbOctCell];
//     balanceCell(southOctCell, iLv);

//     // north-west cell
//     northOctCell = cellNb[3][ngbOctCell];
//     balanceCell(northOctCell, iLv);
//   }

//   // east cell
//   ngbOctCell = cellNb[1][prCell];
//   balanceCell(ngbOctCell, iLv);

//   if (octLv[ngbOctCell / 4] == iLv - 1)
//   {
//     if (morton_lookup[2][ngbOctCell % 4] > 3)
//     {
//       // south-east cell
//       southOctCell = cellNb[2][ngbOctCell];
//       balanceCell(southOctCell, iLv);
//     }
//     if (morton_lookup[3][ngbOctCell % 4] > 3)
//     {
//       // south-east cell
//       northOctCell = cellNb[3][ngbOctCell];
//       balanceCell(southOctCell, iLv);
//     }
//   }
//   else
//   {
//     // south-east cell
//     southOctCell = cellNb[2][ngbOctCell];
//     balanceCell(southOctCell, iLv);

//     // north-east cell
//     northOctCell = cellNb[3][ngbOctCell];
//     balanceCell(northOctCell, iLv);
//   }

//   // south cell
//   ngbOctCell = cellNb[2][prCell];
//   balanceCell(ngbOctCell, iLv);

//   // north cell
//   ngbOctCell = cellNb[3][prCell];
//   balanceCell(ngbOctCell, iLv);
// }
