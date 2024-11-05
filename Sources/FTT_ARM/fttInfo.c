#include <stdio.h>
#include <stdlib.h>
#include "ftt.h"
#include "nrutil.h"

void fttInfo(void)
{
  printf("FTT Mesh Information \n");
  printf("Number Of Cells = %d \n Number Of Octs = %d \n", 
          numberOfCells, numberOfOcts);
}
