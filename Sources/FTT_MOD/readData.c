#include <stdio.h>
#include <stdlib.h>
#include "ftt.h"
#include "pfplib.h"

void readData(void)
{
  // ftt tree
  if (ifetch("maxLevel", &maxLevel) != 1)
  {
    printf("missing ftt maxLevel\n");
    exit(1);
  }
  if (ifetch("minLevel", &minLevel) != 1)
  {
    printf("missing ftt minLevel\n");
    exit(1);
  }
  if (ifetch("maxNumberOfOcts", &maxNumberOfOcts) != 1)
  {
    printf("missing ftt maxNumberOfOcts\n");
    exit(1);
  }
  if (ifetch("maxNumberOfCirclePoints", &maxNumberOfCirclePoints) != 1)
  {
    printf("missing ftt maxNumberOfCirclePoints\n");
    exit(1);
  }

  maxNumberOfCells = maxNumberOfOcts * cellNumberInOct;
  printf("Max Number Of Octs in FTT: %d\n", maxNumberOfOcts);
  printf("Max Number Of Cells in FTT: %d\n", maxNumberOfCells);
  printf("An Oct is composed of %d cells\n", cellNumberInOct);
  printf("Max Level of FTT: %d\n", maxLevel);
  printf("Max Number of Circles resolution : %d\n", maxNumberOfCirclePoints);

  // physical domain
  if (dfetch("Lx", &Lx) != 1)
  {
    printf("missing Lx\n");
    exit(1);
  }
  if (dfetch("Ly", &Ly) != 1)
  {
    printf("missing Ly\n");
    exit(1);
  }
  printf("physical domain: Lx = %g, Ly = %g\n", Lx, Ly);

  return;
}
