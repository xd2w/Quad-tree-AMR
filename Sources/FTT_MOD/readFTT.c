#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "ftt.h"
#include "nrutil.h"

void readFTT(FILE *fp)
{
  int iCell, iOct, iLv, child, parent, level, type, nb;
  char label[128];
  
  fscanf(fp,"%s %d\n", label, &maxLevel);
  if(strcmp(label,"maxLevel") != 0)
  {
    printf("error in mesh reading: maxLevel != %s\n", label);
    exit(1);
  }
  fscanf(fp,"%s %d\n", label, &numberOfCells);
  if(strcmp(label,"numberOfCells") != 0)
  {
    printf("error in mesh reading: numberOfCells != %s\n", label);
    exit(1);
  }
  fscanf(fp,"%s %d\n", label, &numberOfOcts);
  if(strcmp(label,"numberOfOcts") != 0)
  {
    printf("error in mesh reading: numberOfOcts != %s\n", label);
    exit(1);
  }
  fscanf(fp, "%d %d %d\n", &maxLevel, &numberOfCells, &numberOfOcts);
  fscanf(fp, "%s %le %le\n", label, &Lx, &Ly);
  if(strcmp(label,"domainSize") != 0)
  {
    printf("error in mesh reading: domainSize\n");
    exit(1);
  }

  /* initialise the memory */
  initMemory();
/* initialization of FTT */
  for(iCell=0; iCell<maxNumberOfCells; iCell++)
  {
    cellChOct[iCell] = 0;
  }


  fscanf(fp, "%s %d\n", label, &nb);
  if(strcmp(label,"cellVariables_Children_Type") != 0)
  {
    printf("error in mesh reading: cellVariables_Children_Type\n");
    exit(1);
  }
  for(iCell=0; iCell<numberOfCells; iCell++)
  {
    fscanf(fp, "%d %d\n", &child, &type);
    cellChOct[iCell] = child;
    cellType[iCell] = type;
  }
  fscanf(fp, "%s %d\n",label, &nb);
  if(strcmp(label,"octVariables_Parent_Level") != 0)
  {
    printf("error in mesh reading: octVariables_Parent_Level\n");
    exit(1);
  }
  for(iOct=0; iOct<numberOfOcts; iOct++)
  {
    fscanf(fp, "%d %d\n", &parent, &level);
    octPrCell[iOct] = parent;
    octLv[iOct] = level;
  }
  computeCoord();
  printf("read ftt mesh: numberOfCells %d\n",numberOfCells);
}
