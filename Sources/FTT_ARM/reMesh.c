#include <stdio.h>
#include <stdlib.h>
#include "ftt.h"
#include <math.h>

extern void splitFlagCellsAtLevel(int level);
extern void oct_PrCellFlagAtLvel(int level);
extern void refineToKappa(void);
extern void refineToKappaAtLevel(int level);

void reMesh(int itNb)
{
  int level;

  setCellInt1DZero(cellFlag);
  setOctInt1DZero(octFlag);
  // flagInterfCells();
  // refineToKappa();
  // flagInterfLeaves();
  plotFTT(200);
  // plotFlagFTT(200);
  printf("flagInterfLeavs done successfully \n");

  //  for(level=maxLevel; level>maxLevel-1; level--)
  for (level = maxLevel; level > minLevel; level--)
  {
    // refineToKappaAtLevel(level);
    setOctInt1DZeroAtLevel(octFlag, level);
    cell_OctFlagAtLevel(level);
    setCellInt1DZeroAtLevel(cellFlag, level);
    propagateOctFlagAtLevel(level);
    // propagateOctFlagAtLevel(level);
    splitFlagCellsAtLevel(level - 1);
    binCollectionAtLevel(level);
    oct_PrCellFlagAtLvel(level);
    printf("before establish NB \n");
    establishNb();
    printf("establish NB done \n");
    checkNb();
    printf("checkNB done \n");
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
  Real fraction, kappa, cc[6][6], list[4], ccp[3][3], checkSum;
  int currentNumberOfCells = numberOfCells;
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
        if (kappa > 0.7 * (iLv - minLevel - 3))
        {
          splitCell(iCell);
          cOct = cellChOct[iCell];

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

          getChildVOF(0, list, iCell);
          checkSum = 0;
          for (i = 0; i < cellNumberInOct; i++)
          {
            vof[4 * cOct + i] = list[i];
            checkSum += list[i];
            for (dir = 0; dir < 4; dir++)
            {
              dest = morton_lookup[dir][4 * cOct + i];
              if (dest < 4)
              {
                cellNb[dir][4 * cOct + i] = 4 * cOct + dest;
              }
              else
              {
                prNbCell = cellNb[dir][iCell];
                if (octLv[prNbCell / 4] > octLv[iCell / 4])
                {
                  // to make later
                  cellNb[dir][4 * cOct + i] = prNbCell;
                }
                cellNb[dir][4 * cOct + i] = prNbCell;
              }
            }
          }

          // printf("sum = %f\n", 0.25 * checkSum);
          if (fabs(0.25 * checkSum - vof[iCell]) > 1e-4)
          {
            printf("*************\n");
            printf("vof missmatch\n\n");
            exit(1);
          }
        }
      }
    }
  }
}

void refineToKappaAtLevel(int level)
{
  int iCell, iOct, iLv, cOct, i, j, dir, dest, prNbCell;
  Real fraction, kappa, cc[6][6], list[4], ccp[3][3], checkSum;
  int currentNumberOfCells = numberOfCells;
  for (iCell = 0; iCell < currentNumberOfCells; iCell++)
  {
    cellFlag[iCell] = 0;
    if (cellChOct[iCell] == 0 && octLv[iCell / 4] == level)
    {
      fraction = vof[iCell];
      if (fraction > 0.0 && fraction < 1.0)
      {
        iOct = iCell / cellNumberInOct;
        iLv = octLv[iOct];
        getCellNgbVOF_6x6(iOct, cc);
        kappa = kappaBarickALELike(iCell, cc);
        if (kappa > 0.7 * (iLv - minLevel - 3))
        {
          splitCell(iCell);
          cOct = cellChOct[iCell];
          octFlag[cOct] = 1;

          getChildVOF(0, list, iCell);
          checkSum = 0;
          for (i = 0; i < cellNumberInOct; i++)
          {
            vof[4 * cOct + i] = list[i];
            checkSum += list[i];
          }

          // printf("sum = %f\n", 0.25 * checkSum);
          if (fabs(0.25 * checkSum - vof[iCell]) > 1e-4)
          {
            printf("*************\n");
            printf("vof missmatch\n\n");
            exit(1);
          }
        }
      }
    }
  }
}