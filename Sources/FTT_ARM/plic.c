#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "ftt.h"
#include "nrutil.h"

// calculated VOF of flagged cells in both x and y dirs
void plic(void)
{

  // x-sweep
  flagInterfCells();
  propagateFlag(0);
  propagateFlag(0);
  // plotFlagFTT(); exit(1);
  computeXVOF();

  // y-sweep
  flagInterfCells();
  propagateFlag(1);
  propagateFlag(1);
  computeYVOF();
}
void computeXVOF(void)
{
  int iCell, leftNgb, rightNgb;
  Real fraction, MinVof = 1.0e-16, MaxVof = 1.0 - MinVof;

  /* compute vof1, vof2, vof3 on flaged cells */
  for (iCell = 0; iCell < numberOfCells; iCell++)
  {
    if (cellFlag[iCell])
    {
      octTreeXSwp(iCell);
    }
  }
  /* add vof1, vof2, vof3 on marked cells which are one layer less than
     flaged cells */
  for (iCell = 0; iCell < numberOfCells; iCell++)
  {
    if (cellMark[iCell])
    {
      leftNgb = cellNb[0][iCell];
      rightNgb = cellNb[1][iCell];
      fraction = work3[leftNgb] + work2[iCell] + work1[rightNgb];
      if (fraction > MaxVof)
        fraction = 1.0;
      if (fraction < MinVof)
        fraction = 0.0;
      vof[iCell] = fraction;
    }
  }
}
void computeYVOF(void)
{
  int iCell, leftNgb, rightNgb;
  Real fraction, MinVof = 1.0e-16, MaxVof = 1.0 - MinVof;

  /* compute vof1, vof2, vof3 on flaged cells */
  for (iCell = 0; iCell < numberOfCells; iCell++)
  {
    if (cellFlag[iCell])
    {
      octTreeYSwp(iCell);
    }
  }
  /* add vof1, vof2, vof3 on marked cells which are one layer less than
     flaged cells */
  for (iCell = 0; iCell < numberOfCells; iCell++)
  {
    if (cellMark[iCell])
    {
      leftNgb = cellNb[2][iCell];
      rightNgb = cellNb[3][iCell];
      fraction = work3[leftNgb] + work2[iCell] + work1[rightNgb];
      if (fraction > MaxVof)
        fraction = 1.0;
      if (fraction < MinVof)
        fraction = 0.0;
      vof[iCell] = fraction;
    }
  }
}

// copy from one array to another
void copyCellInt1D(Int1D from, Int1D to)
{
  int iCell;
  for (iCell = 0; iCell < numberOfCells; iCell++)
  {
    to[iCell] = from[iCell];
  }
}

// flag cells at the interface of 2 fluids
void flagInterfCells(void)
{
  int iCell, iOct, iLv;
  Real fraction;
  for (iCell = 0; iCell < numberOfCells; iCell++)
  {
    cellFlag[iCell] = 0;
    iOct = iCell / cellNumberInOct;
    iLv = octLv[iOct];
    if (iLv == maxLevel)
    {
      fraction = vof[iCell];
      if (fraction > 0.0 && fraction < 1.0)
      {
        cellFlag[iCell] = 1;
      }
    }
  }
}

/* propagate flag in direction dir
(flags the neighboring cells of the flagged cells in direction x or y)
*/
void propagateFlag(int dir)
{
  int iCell, ngbCell;

  copyCellInt1D(cellFlag, cellMark);
  for (iCell = 0; iCell < numberOfCells; iCell++)
  {
    if (cellMark[iCell]) // cellMark is Copy of cellFlag
    {
      if (dir == 0)
      {
        // left 0  right 1
        ngbCell = cellNb[0][iCell];
        cellFlag[ngbCell] = 1;
        ngbCell = cellNb[1][iCell];
        cellFlag[ngbCell] = 1;
      }
      else
      {
        // down 2  up 3
        ngbCell = cellNb[2][iCell];
        cellFlag[ngbCell] = 1;
        ngbCell = cellNb[3][iCell];
        cellFlag[ngbCell] = 1;
      }
    }
  }
}
