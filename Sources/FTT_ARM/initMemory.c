#include <stdio.h>
#include <stdlib.h>
#include "ftt.h"
#include "nrutil.h"

/*
 Four adjacent cells form an Oct.
 Initialize the first Oct which is composed of the first four cells.
 This Oct is at level 1.  Each cell is given coordinates at the
 south-west corner but no child Oct: cellChOct[0,1,2,3] = 0.
 The first Oct will never have parent cell: octPrCell[0] is not defined.
*/

void initMemory(void)
{
  int iCell, iLv;
  Real dx, dy;

  /* global parameters defined FTT */
  maxNumberOfCells = maxNumberOfOcts * cellNumberInOct;
  dxCell = dvector(0, maxLevel);
  dyCell = dvector(0, maxLevel);
#if (ocTree) /* 3D */
  dzOct = dvector(0, maxLevel);
#endif

  /* ftt oct variables */
  octFlag = ivector(0, maxNumberOfOcts);
  octMark = ivector(0, maxNumberOfOcts);
  octLv = ivector(0, maxNumberOfOcts);
  octPrCell = ivector(0, maxNumberOfOcts);
  octNb = imatrix(0, nbNumberOfOct - 1, 0, maxNumberOfOcts);
  /* ftt cell variable */
  cellChOct = ivector(0, maxNumberOfCells);
  cellFlag = ivector(0, maxNumberOfCells);
  cellMark = ivector(0, maxNumberOfCells);
  cellType = ivector(0, maxNumberOfCells);
  cellHilb = ivector(0, maxNumberOfCells);
  cellNb = imatrix(0, nbNumberOfOct - 1, 0, maxNumberOfCells);

  /* geometric quantities of cell */
  xCell = dvector(0, maxNumberOfCells);
  yCell = dvector(0, maxNumberOfCells);

  mxCell = dvector(0, maxNumberOfCells);
  mzCell = dvector(0, maxNumberOfCells);
  alphaCell = dvector(0, maxNumberOfCells);
#if (ocTree) /* 3D */
  zCell = dvector(0, maxNumberOfCells);
#endif

  /* cell related physical quantities */
  u = dvector(0, maxNumberOfCells);
  v = dvector(0, maxNumberOfCells);
  p = dvector(0, maxNumberOfCells);
  U = dvector(0, maxNumberOfCells);
  V = dvector(0, maxNumberOfCells);
  dive = dvector(0, maxNumberOfCells);
  vof = dvector(0, maxNumberOfCells);
  // work1 = dvector(0, maxNumberOfCells);
  // work2 = dvector(0, maxNumberOfCells);
  // work3 = dvector(0, maxNumberOfCells);

  temp_vof = dvector(0, maxNumberOfCells);

#if (ocTree) /* 3D */
  w = dvector(0, maxNumberOfCells);
  W = dvector(0, maxNumberOfCells);
#endif

  return;
}
