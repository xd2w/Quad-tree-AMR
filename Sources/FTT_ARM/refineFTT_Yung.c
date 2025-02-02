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

Real equation_analytical_curvature(Real x, Real y, Real xc, Real yc);
Real equation_val(Real x, Real y, Real xc, Real yc, Real radius);
int devide_the_cell(struct box b, Real xc, Real yc, Real radius);

void refineFTT_ellipse(void)
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
    if(devide_the_cell(rectangle, xc, yc, radius))
    {
      // split cell 
      printf("%f", equation_analytical_curvature(xCell[iCell] + dx/2, yCell[iCell] + dy/2, xc, yc));
      if(cellChOct[iCell] == 0 && cellType[iCell] == 0) splitCell(iCell);
    }    
  }
  return;
}

int devide_the_cell(struct box b, Real xc, Real yc, Real radius){
  // needs to be not completely in or out
    if (equation_val(b.pt1.x, b.pt1.y, xc, yc, radius) < 0)
      if (equation_val(b.pt2.x, b.pt1.y, xc, yc, radius) < 0)
        if (equation_val(b.pt1.x, b.pt2.y, xc, yc, radius) < 0)
          if (equation_val(b.pt2.x, b.pt2.y, xc, yc, radius) < 0)
            return 0;

    if (equation_val(b.pt1.x, b.pt1.y, xc, yc, radius) > 0)
      if (equation_val(b.pt2.x, b.pt1.y, xc, yc, radius) > 0)
        if (equation_val(b.pt1.x, b.pt2.y, xc, yc, radius) > 0)
          if (equation_val(b.pt2.x, b.pt2.y, xc, yc, radius) > 0)
            return 0;

    return 1;
}

Real equation_val(Real x, Real y, Real xc, Real yc, Real radius){
  return init_VOF_coefs[0] * (x - xc)*(x - xc) \
          + init_VOF_coefs[1] * (y - yc)*(x - xc)\
          + init_VOF_coefs[2] * (y - yc)*(y - yc)\
          - radius*radius;
}

Real equation_analytical_curvature(Real x, Real y, Real xc, Real yc){
  Real a, b, c;
  a = init_VOF_coefs[0];
  b = init_VOF_coefs[1];
  c = init_VOF_coefs[2];

  x -= xc;
  y -= yc;

  return fabs(
        2*a*(b * x + 2 * c * y)*(b * x + 2 * c * y)\
        + 2*c*(2 * a * x + b * y) *(2 * a * x + b * y)\
        - 2 * b * (b * x + 2 * c * y) * (2 * a * x + b * y)\
    ) / pow((2 * a * x + b * y)*(2 * a * x + b * y) + (2 * c * y + b * x)*(2 * c * y + b * x), 3 / 2);
}
