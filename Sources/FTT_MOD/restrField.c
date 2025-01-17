#include <stdio.h>
#include <stdlib.h>
#include "ftt.h"
#include "nrutil.h"

void restrField(Real1D field)
{
  int it;

  // return;
  for (it = maxLevel - 1; it >= 0; it--)
  {
    restrFieldAtLevel(field, it);
  }
}

// calculates val of field for parent by averaging child values
void restrFieldAtLevel(Real1D field, int lev)
{
  int iCell, iLv, iOct, iFlag, chCell;
  Real value;

  for (iCell = 0; iCell < numberOfCells; iCell++)
  {
    iOct = iCell / cellNumberInOct;
    iLv = octLv[iOct];
    chCell = 4 * cellChOct[iCell];
    if (iLv == lev && chCell != 0)
    {
      value = field[chCell];
      chCell++;
      value += field[chCell];
      chCell++;
      value += field[chCell];
      chCell++;
      value += field[chCell];
      field[iCell] = 0.25 * value;
    }
  }
}
