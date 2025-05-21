#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "ftt.h"
#include "nrutil.h"

void plic(void)
{
  // calcPlicPramForAll();
  // x-sweep
  flagInterfLeaves();
  propagateFlag(0);
  // plotFlagFTT(100);
  propagateFlag(0);
  // plotFlagFTT(101);
  // exit(1);
  computeXVOF();
  restrField(vof);

  // y-sweep
  flagInterfLeaves();
  propagateFlag(1);
  // plotFlagFTT(200);
  propagateFlag(1);
  // plotFlagFTT(201);
  // exit(0);
  computeYVOF();
  restrField(vof);
}
void computeXVOF(void)
{
  int h, iCell, leftNgb, rightNgb, leftLv, rightLv, myLv;
  Real fraction, MinVof = 1.0e-16, MaxVof = 1.0 - MinVof;
  Real flux1, flux2, flux3;

  /* compute vof1, vof2, vof3 on flaged cells */
  for (iCell = cellHilb[h]; h < numberOfCells; h++)
  {
    if (cellFlag[iCell])
    {
      octTreeXSwp(iCell);
    }
  }
  /* add vof1, vof2, vof3 on marked cells which are one layer less than
     flaged cells */
  for (iCell = cellHilb[h]; h < numberOfCells; h++)
  {
    if (cellMark[iCell])
    {
      fraction = temp_vof[iCell];
      if (fraction > MaxVof)
        fraction = 1.0;
      if (fraction < MinVof)
        fraction = 0.0;
      vof[iCell] = fraction;
    }
    temp_vof[iCell] = 0;
  }
}
void computeYVOF(void)
{
  int h, iCell, leftNgb, rightNgb;
  Real fraction, MinVof = 1.0e-16, MaxVof = 1.0 - MinVof;

  /* compute vof1, vof2, vof3 on flaged cells */
  for (iCell = cellHilb[h]; h < numberOfCells; h++)
  {
    if (cellFlag[iCell])
    {
      octTreeYSwp(iCell);
    }
  }
  /* add vof1, vof2, vof3 on marked cells which are one layer less than
     flaged cells */
  for (iCell = cellHilb[h]; h < numberOfCells; h++)
  {
    if (cellMark[iCell])
    {
      fraction = temp_vof[iCell];
      if (fraction > MaxVof)
        fraction = 1.0;
      if (fraction < MinVof)
        fraction = 0.0;
      vof[iCell] = fraction;
    }
    temp_vof[iCell] = 0;
  }
}
void copyCellInt1D(Int1D from, Int1D to)
{
  int iCell;
  for (iCell = 0; iCell < numberOfCells; iCell++)
  {
    to[iCell] = from[iCell];
  }
}

void flagInterfCells(void)
{
  int iCell, iOct, iLv;
  Real fraction;
  for (iCell = 0; iCell < numberOfCells; iCell++)
  {
    cellFlag[iCell] = 0;
    iOct = iCell / cellNumberInOct;
    iLv = octLv[iOct];
    // if(iLv == maxLevel)
    {
      fraction = vof[iCell];
      if (fraction > 0.0 && fraction < 1.0)
      {
        cellFlag[iCell] = 1;
      }
    }
  }
}

void flagInterfLeaves(void)
{
  int iCell, iOct, iLv;
  Real fraction;
  for (iCell = 0; iCell < numberOfCells; iCell++)
  {
    cellFlag[iCell] = 0;
    if (cellChOct[iCell] == 0)
    {
      fraction = vof[iCell];
      if (fraction > 0.0 && fraction < 1.0)
      {
        cellFlag[iCell] = 1;
      }
    }
  }
}

/* propagate flag in direction dir */
void propagateFlag(int dir)
{
  int iCell, ngbCell;

  copyCellInt1D(cellFlag, cellMark);
  for (iCell = 0; iCell < numberOfCells; iCell++)
  {
    if (cellMark[iCell])
    {
      if (dir == 0)
      {
        ngbCell = cellNb[0][iCell];
        if (cellChOct[ngbCell] == 0)
        {
          cellFlag[ngbCell] = 1;
        }
        else
        {
          cellFlag[4 * cellChOct[ngbCell] + 1] = 1;
          cellFlag[4 * cellChOct[ngbCell] + 3] = 1;
        }

        ngbCell = cellNb[1][iCell];
        if (cellChOct[ngbCell] == 0)
        {
          cellFlag[ngbCell] = 1;
        }
        else
        {
          cellFlag[4 * cellChOct[ngbCell] + 0] = 1;
          cellFlag[4 * cellChOct[ngbCell] + 2] = 1;
        }
      }
      else
      {
        ngbCell = cellNb[2][iCell];
        if (cellChOct[ngbCell] == 0)
        {
          cellFlag[ngbCell] = 1;
        }
        else
        {
          cellFlag[4 * cellChOct[ngbCell] + 2] = 1;
          cellFlag[4 * cellChOct[ngbCell] + 3] = 1;
        }

        ngbCell = cellNb[3][iCell];
        if (cellChOct[ngbCell] == 0)
        {
          cellFlag[ngbCell] = 1;
        }
        else
        {
          cellFlag[4 * cellChOct[ngbCell] + 0] = 1;
          cellFlag[4 * cellChOct[ngbCell] + 1] = 1;
        }
      }
    }
  }
}

/* propagate flag in direction dir */
void propagateFlagAtLevel(int dir, int level)
{
  int h, iCell, ngbCell;

  copyCellInt1D(cellFlag, cellMark);
  for (iCell = 0; iCell < numberOfCells; iCell++)
  {
    if (cellMark[iCell] && octLv[iCell / 4] == level)
    {
      if (dir == 0)
      {
        ngbCell = cellNb[0][iCell];
        if (cellChOct[ngbCell] == 0)
        {
          cellFlag[ngbCell] = 1;
        }
        else
        {
          cellFlag[4 * cellChOct[ngbCell] + 1] = 1;
          cellFlag[4 * cellChOct[ngbCell] + 3] = 1;
        }

        ngbCell = cellNb[1][iCell];
        if (cellChOct[ngbCell] == 0)
        {
          cellFlag[ngbCell] = 1;
        }
        else
        {
          cellFlag[4 * cellChOct[ngbCell] + 0] = 1;
          cellFlag[4 * cellChOct[ngbCell] + 2] = 1;
        }
      }
      else
      {
        ngbCell = cellNb[2][iCell];
        if (cellChOct[ngbCell] == 0)
        {
          cellFlag[ngbCell] = 1;
        }
        else
        {
          cellFlag[4 * cellChOct[ngbCell] + 2] = 1;
          cellFlag[4 * cellChOct[ngbCell] + 3] = 1;
        }

        ngbCell = cellNb[3][iCell];
        if (cellChOct[ngbCell] == 0)
        {
          cellFlag[ngbCell] = 1;
        }
        else
        {
          cellFlag[4 * cellChOct[ngbCell] + 0] = 1;
          cellFlag[4 * cellChOct[ngbCell] + 1] = 1;
        }
      }
    }
  }
}

void initialise_temp_vof_to(Real val)
{
  int iCell;
  for (iCell = 0; iCell < numberOfCells; iCell++)
  {
    temp_vof[iCell] = val;
  }
}