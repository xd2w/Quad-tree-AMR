#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "ftt.h"
#include "nrutil.h"
#include "box.h"
#include "pfplib.h"

void refineFTT(void)
{
  int cLv, iOct, iCell, chCell;
  Real xPoint, yPoint, theta, delta;
  Real dx, dy, Lx, Ly, radius, xc, yc, xShock;
  struct point *pList;
  struct box rectangle;
  int n, nPoint, nCircle, nShock;
  FILE *fPoints; 

   xc = 0.5;
   yc = 0.5;
   Lx = 1.0;
   Ly = 1.0;
   radius = 0.2;
   xShock = 0.2;
   dfetch("xc", &xc); 
   dfetch("yc", &yc); 
   dfetch("radius", &radius); 
   dfetch("xShock", &xShock); 
   dx = dxCell[maxLevel];
   dy = dyCell[maxLevel];
   dfetch("Lx", &Lx); 
   dfetch("Ly", &Ly); 
   nShock = 2*Ly/dy;
   nShock = 0;
   nCircle = 4*3.14159*radius/dx;
   delta = 2*3.14159/nCircle;
   nPoint = nCircle + nShock; 
// read a list of points 
  pList = (struct point *) malloc(nPoint*sizeof(struct point));
  for(n=0; n<nCircle; n++)
  {
    theta=(n+.5)*delta;
    pList[n].x = xc+radius*cos(theta); pList[n].y = yc+radius*sin(theta);
  }
/*  delta = Ly/nShock;
  for(n=nCircle; n<nPoint; n++)
  {
    pList[n].x = xShock; pList[n].y = (n+.5-nCircle)*delta;
  }
*/
    fPoints = fopen("meshSeeds", "w");
   for(n=0; n<nCircle; n++)
   {
     fprintf(fPoints, "%g %g\n", pList[n].x, pList[n].y);
   }
     fprintf(fPoints, "%g %g\n", pList[0].x, pList[0].y);
   fprintf(fPoints, "\n");
   fclose(fPoints);
   /*
   for(n=nCircle; n<nPoint; n++)
   {
     fprintf(fPoints, "%g %g\n", pList[n].x, pList[n].y);
   }
   */
  for(iCell=0; iCell<numberOfCells; iCell++)
  {
//    printf("split cell %d, numberOfCells %d\n", iCell, numberOfCells);
//    getchar();

    iOct = iCell/cellNumberInOct;
    cLv = octLv[iOct];
    if(cLv >= maxLevel) return;

    dx = dxCell[cLv]; 
    dy = dyCell[cLv]; 
    rectangle.pt1.x = xCell[iCell];
    rectangle.pt1.y = yCell[iCell];
    rectangle.pt2.x = xCell[iCell] + dx;
    rectangle.pt2.y = yCell[iCell] + dy;
    /* point within cell, split the cell */
    //if(ptListIntersectBox(pList, nPoint, rectangle) > 1)
    if(ptListIntersectBox(pList, nPoint, rectangle))
    {
      // split cell 
      if(cellChOct[iCell] == 0 && cellType[iCell] == 0) splitCell(iCell);
    }    
  }
  return;
}
