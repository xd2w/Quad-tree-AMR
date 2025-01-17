#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "ftt.h"
#include "nrutil.h"

void plic(void)
{
  flagInterfCells();
  propagateFlag(0);
  copyOctInt(cellFlag, cellMark);
  propagateFlag(0);
  plotFlagFTT();
  computeXVOF();
}
void computeXVOF()
{
  int iCell, leftNgb, rightNgb;
  Real fraction, MinVof = 1.0e-16, MaxVof = 1.0 - MinVof;
 
/* compute vof1, vof2, vof3 on flaged cells */ 
  for(iCell=0; iCell<numberOfCells; iCell++)
  {
    if(cellFlag[iCell])
    {
      octTreeXSwp(iCell);
    }
  }
/* add vof1, vof2, vof3 on marked cells which are one layer less than 
   flaged cells */
  for(iCell=0; iCell<numberOfCells; iCell++)
  {
    if(cellMark[iCell])
    {
      leftNgb = cellNb[0][iCell];
      rightNgb = cellNb[1][iCell];
      fraction = work3[leftNgb]+work2[iCell]+work1[rightNgb];
      if(fraction > MaxVof) fraction = 1.0; 
      if(fraction < MinVof) fraction = 0.0; 
      vof[iCell] = fraction;
    }
  }

}
void  copyOctInt(Int1D from, Int1D to)
{
  int iCell;
  for(iCell=0; iCell<numberOfCells; iCell++)
  {
    to[iCell] = from[iCell];
  }
}

void flagInterfCells(void)
{
  int iCell, iOct, iLv;
  Real fraction;

  for(iCell=0; iCell<numberOfCells; iCell++)
  {
    cellFlag[iCell] = 0;
    iOct = iCell/cellNumberInOct;
    iLv = octLv[iOct];
    if(iLv == maxLevel)
    {
      fraction = vof[iCell];
      if(fraction > 0.0 && fraction <1.0)
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

  copyOctInt(cellFlag, cellMark);
  for(iCell=0; iCell<numberOfCells; iCell++)
  {
    if(cellMark[iCell])
    {
      if(dir == 0)
      {
        ngbCell = cellNb[0][iCell];
        cellFlag[ngbCell] = 1;
        ngbCell = cellNb[1][iCell];
        cellFlag[ngbCell] = 1;
      }
      else
        ngbCell = cellNb[2][iCell];
        cellFlag[ngbCell] = 1;
        ngbCell = cellNb[3][iCell];
        cellFlag[ngbCell] = 1;
      {
      } 
    } 
  }
}
