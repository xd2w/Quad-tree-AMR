#include <stdio.h>
#include <stdlib.h>
#include "ftt.h"
#include "nrutil.h"

void plotSFC(int ndata)
{
  int iCell, i, i1, i2, i3;
  char fsfcv[] = "DATA/sfcv.000";
  FILE *fsfc;
  i = ndata;
  i1 = i % 10;
  i /= 10;
  i2 = i % 10;
  i /= 10;
  i3 = i % 10;
  fsfcv[10] = '0' + i3;
  fsfcv[11] = '0' + i2;
  fsfcv[12] = '0' + i1;

  fsfc = fopen(fsfcv, "w");
  for (iCell = 0; iCell < cellNumberInOct; iCell++)
  {
    plotFTTCellSFC(iCell, fsfc);
  }
  fclose(fsfc);

  return;
}

void plotFTTCellSFC(int iCell, FILE *fsfc)
{
  int iOct, cLv;
  Real xLeft, yLeft, xRight, yRight;

  // do not draw if cell is not a leaf
  //  if(cellFlag[iCell]) return;

  iOct = iCell / cellNumberInOct;
  cLv = octLv[iOct];
  // printf("iCell %d iOct %d cLv %d\n", iCell, iOct, cLv);
  xLeft = xCell[iCell];
  yLeft = yCell[iCell];

  Real xx = xLeft + .5 * dxCell[cLv];
  Real yy = yLeft + .5 * dyCell[cLv];

  /* only plot on leaves of the tree */
  iOct = cellChOct[iCell];
  if (iOct == 0)
  {
    fprintf(fsfc, "%g %g\n", xx, yy);
  }
  else
  {
    iCell = iOct * cellNumberInOct;
    plotFTTCellSFC(iCell, fsfc);
    plotFTTCellSFC(iCell + 1, fsfc);
    plotFTTCellSFC(iCell + 2, fsfc);
    plotFTTCellSFC(iCell + 3, fsfc);
  }
  return;
}
void plotHilbertSFC(int ndata)
{
  int iCell, iOct, i, i1, i2, i3;
  char hilbv[] = "DATA/hilb.000";
  FILE *fhil;
  i = ndata;
  i1 = i % 10;
  i /= 10;
  i2 = i % 10;
  i /= 10;
  i3 = i % 10;
  hilbv[10] = '0' + i3;
  hilbv[11] = '0' + i2;
  hilbv[12] = '0' + i1;

  fhil = fopen(hilbv, "w");

  iOct = 0;
  iCell = iOct * cellNumberInOct;
  int hType = 0;
  for (int cell = 0; cell < cellNumberInOct; cell++)
  {
    int ic = hilbert_map[hType][cell];
    cellHilb[iCell + cell] = hilbert_production[hType][cell];
    // printf("prod %d\n", hilbert_production[0][ic]);
    plotFTTCellHilbert(iCell + ic, fhil);
  }

  fclose(fhil);

  return;
}

void plotFTTCellHilbert(int iCell, FILE *fsfc)
{
  int iOct, cLv;
  Real xLeft, yLeft, xRight, yRight;

  // do not draw if cell is not a leaf
  //  if(cellFlag[iCell]) return;

  int hType = cellHilb[iCell];
  iOct = iCell / cellNumberInOct;
  cLv = octLv[iOct];
  if (iOct == 4)
    printf("iCell %d iOct %d cLv %d hType %d\n", iCell, iOct, cLv, hType);
  xLeft = xCell[iCell];
  yLeft = yCell[iCell];

  Real xx = xLeft + .5 * dxCell[cLv];
  Real yy = yLeft + .5 * dyCell[cLv];

  /* only plot on leaves of the tree */
  iOct = cellChOct[iCell];
  if (iOct == 0)
  {
    fprintf(fsfc, "%g %g\n", xx, yy);
  }
  else
  {
    iCell = iOct * cellNumberInOct;
    for (int cell = 0; cell < cellNumberInOct; cell++)
    {
      int ic = hilbert_map[hType][cell];
      cellHilb[iCell + ic] = hilbert_production[hType][ic];
      plotFTTCellHilbert(iCell + ic, fsfc);
    }
  }
  return;
}

int _cHilbup(int iCell, int index, int hType)
{ // recursive method of hilbert curve production
  int i, ic, iOct, n;
  if (cellChOct[iCell] == 0)
  {
    return 1;
  }
  else
  {
    iOct = cellChOct[iCell];
    n = index;
    for (i = 0; i < 4; i++)
    {
      ic = hilbert_map[hType][i];
      cellMark[4 * iOct + ic] += n;
      n += _cHilbup(4 * iOct + ic, n, hilbert_production[hType][i]);
    }
    // cellMark[iCell] = n;
    // return n;
  }
}

int _cHilbdown(int iCell, int hType)
{ // recursive method of hilbert curve production
  int i, ic, iOct, n;
  if (cellChOct[iCell] == 0)
  {
    return 1;
  }
  else
  {
    iOct = cellChOct[iCell];
    n = cellMark[iCell];
    for (i = 0; i < 4; i++)
    {
      ic = hilbert_map[hType][i];
      cellMark[4 * iOct + ic] += n;
    }
  }
}

// int _cHilb2(int iCell, int **index, int hType)
// { // recursive method of hilbert curve production
//   int i, ic, iOct, n;
//   if (cellChOct[iCell] == 0)
//   {
//     return cellMark[iCell] = index++;
//   }
//   else
//   {
//     iOct = octPrCell[iCell];
//     n = index;
//     for (i = 0; i < 4; i++)
//     {
//       ic = hilbert_map[hType][i];
//       cellMark[4 * iOct + ic] += n;
//       n += _cHilb(4 * iOct + ic, n, hilbert_production[hType][i]);
//     }
//   }
// }

void _cHilbParallel(void)
{ // dynamic programming method of hilbert curve production
  // paralaisable easily but inefficent
  int i, ic, iOct, n, temp;

  for (i = 0; i < numberOfCells; i++)
  {
    cellMark[i] = 0;
    if (cellChOct[i] == 0)
      cellMark[i] = 1;
  }

  for (int level = maxLevel; level > 0; level--)
  {
    for (iOct = 0; iOct < numberOfOcts; iOct++)
    {
      if (octLv[iOct] == level)
      {
        n = 0;
        for (i = 0; i < cellNumberInOct; i++)
        {
          n += cellMark[4 * iOct + i];
        }
        cellMark[octPrCell[iOct]] = n;
      }
    }
  }

  cellHilb[0] = 0;
  cellMark[0] = 0;

  for (int level = 1; level < maxLevel; level++)
  {
    for (iOct = 1; iOct < numberOfOcts; iOct++)
    {
      if (octLv[iOct] == level)
      {
        n = cellMark[octPrCell[iOct]];
        for (i = 0; i < cellNumberInOct; i++)
        {
          ic = hilbert_map[cellHilb[octPrCell[iOct]]][i];
          cellHilb[4 * iOct + ic] = hilbert_production[cellHilb[octPrCell[iOct]]][i];
          temp = cellMark[4 * iOct + ic];
          cellMark[4 * iOct + ic] += n;
          n += temp;
        }
      }
    }
  }
}

void constructHilbertIndex(void)
{
  int iCell, iOct, i, ic, level;
  // _cHilbParallel();
  setCellInt1DZero(cellMark);
  _cHilbup(0, 0, 0);
  // _cHilbdown(0, 0);
  printf("%d\n\n", cellMark[1]);
  int numLeaves = 0;
  for (iCell = 0; iCell < numberOfCells; iCell++)
  {
    if (cellChOct[iCell] == 0)
    {
      printf("%d\n", cellMark[iCell]);
      if (cellMark[iCell] > numberOfCells)
      {
        printf("somthing wrong\n");
        exit(1);
      }
      cellHilb[cellMark[iCell]] = iCell;
      numLeaves++;
    }
  }
  printf("num leaves : %d\n", numberOfCells);
  exit(1);
}
