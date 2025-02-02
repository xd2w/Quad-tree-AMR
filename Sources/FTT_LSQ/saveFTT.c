#include <stdio.h>
#include <stdlib.h>
#include "ftt.h"
#include "nrutil.h"

void saveFTT(FILE *fp)
{
  int iCell, iOct, i, i1, i2, i3;

  fprintf(fp, "maxLevel %d \nnumberOfCells %d  \nnumberOfOcts %d\n", 
              maxLevel, numberOfCells, numberOfOcts);
  fprintf(fp, "domainSize %g %g\n", Lx, Ly);

  fprintf(fp, "cellVariables_Children_Type %d\n", numberOfCells);
  for(iCell=0; iCell<numberOfCells; iCell++)
  {
    fprintf(fp, "%d %d\n", cellChOct[iCell], cellType[iCell]);
  }
  fprintf(fp, "octVariables_Parent_Level %d\n", numberOfOcts);
  for(iOct=0; iOct<numberOfOcts; iOct++)
  {
    fprintf(fp, "%d %d\n", octPrCell[iOct], octLv[iOct]);
  }
  printf("save ftt mesh: numberOfCells %d\n",numberOfCells);
}
