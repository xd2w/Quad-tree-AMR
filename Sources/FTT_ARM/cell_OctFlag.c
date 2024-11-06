#include <stdio.h>
#include <stdlib.h>
#include "ftt.h"
#include "nrutil.h"

// set all val of cell to 0
void setCellInt1DZero(Int1D val)
{
  int iCell;
  for (iCell = 0; iCell < numberOfCells; iCell++)
  {
    val[iCell] = 0;
  }
}

// set all val of cell to 0 at a given level
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

// sets all val of oct to 0
void setOctInt1DZero(Int1D val)
{
  int iOct;
  for (iOct = 0; iOct < numberOfOcts; iOct++)
  {
    val[iOct] = 0;
  }
}

// set all val of Oct to 0 at a given level
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
  int iCell, chOct;
  for (iCell = 0; iCell < numberOfCells; iCell++)
  {
    // no kid + flagged
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
  int iCell, chOct, iOct;
  for (iCell = 0; iCell < numberOfCells; iCell++)
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

// if any cell is flagged, then oct is flagged
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

// if any cell is flagged, then oct is flagged at a given level
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

// flags all prCells of octs
void oct_PrCellFlag(void)
{
  int iOct, prCell;
  for (iOct = 0; iOct < numberOfOcts; iOct++)
  {
    prCell = octPrCell[iOct];
    cellFlag[prCell] = 1;
  }
}

// flags all prCells of octs at a given level
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

// copy from one Oct values to another array
void copyOctInt1D(Int1D from, Int1D to)
{
  int iOct;
  for (iOct = 0; iOct < numberOfOcts; iOct++)
  {
    to[iOct] = from[iOct];
  }
}

// if oct is flagged the the neighbouring oct/cells are flagged
void propagateOctFlag(void)
{
  int iCell, iOct, prCell, chOct, ngbCell;

  copyOctInt1D(octFlag, octMark);

  for (iOct = 0; iOct < numberOfOcts; iOct++)
  {
    if (octMark[iOct]) // if oct is flagged
    {
      prCell = octPrCell[iOct];
      {
        ngbCell = cellNb[0][prCell]; // left (west) neighbor of parent cell
        chOct = cellChOct[ngbCell];  // child oct of left neighbor
        if (chOct > 0)               // if child exist?
        {
          octFlag[chOct] = 1; // flag the child oct as well
        }
        else // if not
        {
          cellFlag[ngbCell] = 1; // flag the neighbor of the paretn neighbor
        }
        // ... for rest of neighbors
        ngbCell = cellNb[1][prCell]; // right (east)
        chOct = cellChOct[ngbCell];
        if (chOct > 0)
        {
          octFlag[chOct] = 1;
        }
        else
        {
          cellFlag[ngbCell] = 1;
        }
        ngbCell = cellNb[2][prCell]; // south
        chOct = cellChOct[ngbCell];
        if (chOct > 0)
        {
          octFlag[chOct] = 1;
        }
        else
        {
          cellFlag[ngbCell] = 1;
        }
        ngbCell = cellNb[3][prCell]; // north
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

/* if oct is flagged then the neighbouring oct/cells are flagged
 * only at a givel level
 */
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
