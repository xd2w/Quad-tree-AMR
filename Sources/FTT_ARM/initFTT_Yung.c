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

void initFTT(void)
{
  int iCell, iLv;
  Real dx, dy;

/* initialization of FTT */
  for(iCell=0; iCell<maxNumberOfCells; iCell++)
  {
    cellChOct[iCell] = 0;
  }
  dx = 2.0*Lx; dy = 2.0*Ly;
  for(iLv=0; iLv<=maxLevel; iLv++)
  {
    dxCell[iLv] = dx;
    dyCell[iLv] = dy;
    dx *= 0.5; dy *= 0.5;
    printf("At level  %d, cell size %g %g\n", iLv, dxCell[iLv], dyCell[iLv]);
  }
/* computing domain [0:Lx][0:Ly] located as the north-east (third) cell of 
  the first Oct, 7th cell */

  /* start with zero-th oct, cells 0, 1, 2, and 3 */
  octLv[0] = 0;   /* start with oct level 0 */
  numberOfOcts = 1; numberOfCells = numberOfOcts*cellNumberInOct;
  cellChOct[0] = 0;
  cellChOct[1] = 0;
  cellChOct[2] = 0;
  cellChOct[3] = 0;
  cellType[0] = 1;
  cellType[1] = 1;
  cellType[2] = 1;
  cellType[3] = 1;
  xCell[0] = -Lx;
  yCell[0] = -Ly;
  xCell[1] = Lx;
  yCell[1] = -Ly;
  xCell[2] = -Lx;
  yCell[2] = Ly;
  xCell[3] = Lx;
  yCell[3] = Ly;
  // split cell 0
  splitCell(0);
  splitCell(1);
  splitCell(2);
  splitCell(3);
  cellType[7] = 0;
//  cellType[13] = 0;
  FILE* fp;
  fp = fopen("computingBox","w");
  plotFTTCell(7, fp);
 // plotFTTCell(13, fp);
  fclose(fp);
/*
  for(iCell=0; iCell<numberOfCells; iCell++)
  {
    printf("iCell %d type %d\n", iCell, cellType[iCell]);
  }
  exit(1);
*/
  return;
}
