#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "ftt.h"
#include "nrutil.h"
#include "box.h"

void refineFTT(void)
{
  int cLv, iOct, iCell, chCell;
  Real xPoint, yPoint, deltaX, deltaY;
  Real dx, dy, Lx, Ly, xc, yc, deltaHeight;
  struct point *pList;
  struct box rectangle;
  int n, nPoint, nbHeight, nbEdge;
  FILE *fPoints; 

   xc = 0.5;
   yc = 0.5;
   deltaHeight=0.0625;
   Lx = 1.0;
   Ly = 1.0;
   dfetch("xc", &xc); 
   dfetch("yc", &yc); 
   dfetch("deltaHeight", &deltaHeight); 
   dx = dxCell[maxLevel];
   dy = dyCell[maxLevel];
   dfetch("Lx", &Lx); 
   dfetch("Ly", &Ly); 
   deltaHeight *= Lx;
   nbHeight = 2.0*deltaHeight/dx;
   deltaX = deltaHeight/nbHeight;
   nbEdge = 4.0*deltaHeight*Ly/Lx/dy;
   deltaY = 2.0*deltaHeight*Ly/Lx/nbEdge;

   nPoint = 2*nbHeight+nbEdge+1; 
// read a list of points 
  pList = (struct point *) malloc(nPoint*sizeof(struct point));
  for(n=0; n<nbHeight; n++)
  {
    pList[n].x = xc+n*deltaX; 
    pList[n].y = yc+n*deltaX*Ly/Lx;
  }
  for(n=0; n<nbHeight; n++)
  {
    pList[n+nbHeight].x = xc+n*deltaX; 
    pList[n+nbHeight].y = yc-n*deltaX*Ly/Lx;
  }
  xc += deltaHeight;
  yc -= deltaHeight*Ly/Lx;
  for(n=0; n<=nbEdge; n++)
  {
    pList[n+2*nbHeight].x = xc; 
    pList[n+2*nbHeight].y = yc+n*deltaY;
  }

    fPoints = fopen("meshSeeds", "w");
   for(n=0; n<nbHeight; n++)
   {
     fprintf(fPoints, "%g %g\n", pList[n].x, pList[n].y);
   }
   fprintf(fPoints, "\n");
   for(n=nbHeight; n<2*nbHeight; n++)
   {
     fprintf(fPoints, "%g %g\n", pList[n].x, pList[n].y);
   }
   fprintf(fPoints, "\n");

   for(n=2*nbHeight; n<nPoint; n++)
   {
     fprintf(fPoints, "%g %g\n", pList[n].x, pList[n].y);
   }

  for(iCell=0; iCell<numberOfCells; iCell++)
  {
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
}
