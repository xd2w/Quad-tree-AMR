#include <stdio.h>
#include <stdlib.h>
#include "ftt.h"
#include "nrutil.h"
#include "box.h"

void fttStatistics(void)
{
  int count, iCell;
  printf("Statistics on a 2D Fully Threaded Tree DATA Sturcture\n");
  printf("Max Level %d\n", maxLevel);
  printf("Number Of Cells = %d Number Of Octs = %d \n", 
                   numberOfCells, numberOfOcts);
  count = 0;
  for(iCell=0; iCell<numberOfCells; iCell++)
  {
    if(cellFlag[iCell]) count++;
  } 
  printf("Flaged cell number = %d, flaged oct number = %d\n", count, count/4);

  return;
}
