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

void computeCoord(void)
{
  int iCell, iLv;
  Real dx, dy;

  dx = 2.0*Lx; dy = 2.0*Ly;
  for(iLv=0; iLv<=maxLevel; iLv++)
  {
    dxCell[iLv] = dx;
    dyCell[iLv] = dy;
    dx *= 0.5; dy *= 0.5;
    printf("At level  %d, cell size %g %g\n", iLv, dxCell[iLv], dyCell[iLv]);
  }
  
/* first cell [-Lx:Lx][-Ly:Ly]: located  as south-west cell of the first Oct */
  xCell[0] = -Lx;
  yCell[0] = -Ly;
  xCell[1] = Lx;
  yCell[1] = -Ly;
  xCell[2] = -Lx;
  yCell[2] = Ly;
  xCell[3] = Lx;
  yCell[3] = Ly;

  for(iLv=1; iLv <= maxLevel; iLv++)
  {
    computeCoordAtLevel(iLv);
  }  
}
void computeCoordAtLevel(int Lv)
{
  int iOct, iLv, octCell, prCell;
  Real dx, dy, xc, yc;

  dx = dxCell[Lv]; 
  dy = dyCell[Lv]; 
  for(iOct=1; iOct<numberOfOcts; iOct++)
  {
    iLv = octLv[iOct];
    if(iLv == Lv)
    {
      octCell = 4*iOct;
      prCell = octPrCell[iOct];
      xc = xCell[prCell];
      yc = yCell[prCell];
      xCell[octCell  ] = xc   ;
      yCell[octCell  ] = yc   ;
      xCell[octCell+1] = xc+dx;
      yCell[octCell+1] = yc   ;
      xCell[octCell+2] = xc   ;
      yCell[octCell+2] = yc+dy;
      xCell[octCell+3] = xc+dx;
      yCell[octCell+3] = yc+dy;
    }
  } 
}
