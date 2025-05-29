#include <stdio.h>
#include <stdlib.h>
#include "ftt.h"
#include "nrutil.h"

void setCellInt1DZero(Int1D val)
{
  int iCell;
  for (iCell = 0; iCell < numberOfCells; iCell++)
  {
    val[iCell] = 0;
  }
}
void setCellInt1DZeroAtLevel(Int1D val, int level)
{
  int iCell, iOct;
  for (iCell = 0; iCell < numberOfCells; iCell++)
  {
    iOct = iCell / cellNumberInOct;
    if (octLv[iOct] == level)
      val[iCell] = 0;
  }
}

void setOctInt1DZero(Int1D val)
{
  int iOct;
  for (iOct = 0; iOct < numberOfOcts; iOct++)
  {
    val[iOct] = 0;
  }
}
void setOctInt1DZeroAtLevel(Int1D val, int level)
{
  int iOct;
  for (iOct = 0; iOct < numberOfOcts; iOct++)
  {
    if (octLv[iOct] == level)
      val[iOct] = 0;
  }
}

void splitFlagCells(void)
{
  int h, iCell, chOct;
  for (iCell = cellHilb[h = 0]; h < numberOfCells; iCell = cellHilb[++h])
  {
    if (cellChOct[iCell] == 0 && cellFlag[iCell])
    {
      splitCell(iCell);
      chOct = cellChOct[iCell];
      octFlag[chOct] = 1;
    }
  }
}
void splitFlagCellsAtLevel(int level)
{
  int h, iCell, chOct, iOct;
  for (iCell = cellHilb[h = 0]; h < numberOfCells; iCell = cellHilb[++h])
  {
    iOct = iCell / cellNumberInOct;
    if (cellChOct[iCell] == 0 && cellFlag[iCell] && octLv[iOct] == level)
    {
      splitCell(iCell);
      chOct = cellChOct[iCell];
      octFlag[chOct] = 1;
    }
  }
}

void cell_OctFlag(void)
{
  int iCell, iOct;

  for (iOct = 0; iOct < numberOfOcts; iOct++)
  {
    //    octFlag[iOct] = 0;
    iCell = cellNumberInOct * iOct;
    if (cellFlag[iCell] == 1)
      octFlag[iOct] = 1;
    iCell++;
    if (cellFlag[iCell] == 1)
      octFlag[iOct] = 1;
    iCell++;
    if (cellFlag[iCell] == 1)
      octFlag[iOct] = 1;
    iCell++;
    if (cellFlag[iCell] == 1)
      octFlag[iOct] = 1;
  }
}
void cell_OctFlagAtLevel(int level)
{
  int iCell, iOct;

  for (iOct = 0; iOct < numberOfOcts; iOct++)
  {
    //    octFlag[iOct] = 0;
    if (octLv[iOct] == level)
    {
      iCell = cellNumberInOct * iOct;
      if (cellFlag[iCell] == 1)
        octFlag[iOct] = 1;
      iCell++;
      if (cellFlag[iCell] == 1)
        octFlag[iOct] = 1;
      iCell++;
      if (cellFlag[iCell] == 1)
        octFlag[iOct] = 1;
      iCell++;
      if (cellFlag[iCell] == 1)
        octFlag[iOct] = 1;
    }
  }
}

void oct_PrCellFlag(void)
{
  int iOct, prCell;
  for (iOct = 0; iOct < numberOfOcts; iOct++)
  {
    prCell = octPrCell[iOct];
    cellFlag[prCell] = 1;
  }
}
void oct_PrCellFlagAtLvel(int level)
{
  int iOct, prCell;

  for (iOct = 0; iOct < numberOfOcts; iOct++)
  {
    if (octLv[iOct] == level)
    {
      prCell = octPrCell[iOct];
      cellFlag[prCell] = 1;
    }
  }
}

void copyOctInt1D(Int1D from, Int1D to)
{
  int iOct;
  for (iOct = 0; iOct < numberOfOcts; iOct++)
  {
    to[iOct] = from[iOct];
  }
}
void propagateOctFlag(void)
{
  int iCell, iOct, prCell, chOct, ngbCell;

  copyOctInt1D(octFlag, octMark);

  for (iOct = 0; iOct < numberOfOcts; iOct++)
  {
    if (octMark[iOct])
    {
      prCell = octPrCell[iOct];
      {
        ngbCell = cellNb[0][prCell];
        chOct = cellChOct[ngbCell];
        if (chOct > 0)
        {
          octFlag[chOct] = 1;
        }
        else
        {
          cellFlag[ngbCell] = 1;
        }
        ngbCell = cellNb[1][prCell];
        chOct = cellChOct[ngbCell];
        if (chOct > 0)
        {
          octFlag[chOct] = 1;
        }
        else
        {
          cellFlag[ngbCell] = 1;
        }
        ngbCell = cellNb[2][prCell];
        chOct = cellChOct[ngbCell];
        if (chOct > 0)
        {
          octFlag[chOct] = 1;
        }
        else
        {
          cellFlag[ngbCell] = 1;
        }
        ngbCell = cellNb[3][prCell];
        chOct = cellChOct[ngbCell];
        if (chOct > 0)
        {
          octFlag[chOct] = 1;
        }
        else
        {
          cellFlag[ngbCell] = 1;
        }
      }
    }
  }
}

void propagateOctFlagAtLevel(int level)
{
  int iCell, iOct, prCell, chOct, ngbCell;

  copyOctInt1D(octFlag, octMark);

  for (iOct = 0; iOct < numberOfOcts; iOct++)
  {
    if (octMark[iOct] && octLv[iOct] == level)
    {
      prCell = octPrCell[iOct];

      ngbCell = cellNb[0][prCell];
      if (ngbCell > 0)
      {
        chOct = cellChOct[ngbCell];
        if (chOct > 0)
        {
          octFlag[chOct] = 1;
        }
        else
        {
          cellFlag[ngbCell] = 1;
        }
      }
      ngbCell = cellNb[1][prCell];
      chOct = cellChOct[ngbCell];
      if (chOct > 0)
      {
        octFlag[chOct] = 1;
      }
      else
      {
        cellFlag[ngbCell] = 1;
      }
      ngbCell = cellNb[2][prCell];
      if (ngbCell > 0)
      {
        chOct = cellChOct[ngbCell];
        if (chOct > 0)
        {
          octFlag[chOct] = 1;
        }
        else
        {
          cellFlag[ngbCell] = 1;
        }
      }
      ngbCell = cellNb[3][prCell];
      if (ngbCell > 0)
      {
        chOct = cellChOct[ngbCell];
        if (chOct > 0)
        {
          octFlag[chOct] = 1;
        }
        else
        {
          cellFlag[ngbCell] = 1;
        }
      }
    }
  }
}
