#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "ftt.h"
#include "nrutil.h"
#include "box.h"
#include "pfplib.h"

Real equation_val(Real x, Real y, Real xc, Real yc, Real radius)
{
  return init_VOF_coefs[0] * (x - xc) * (x - xc) + init_VOF_coefs[1] * (y - yc) * (x - xc) + init_VOF_coefs[2] * (y - yc) * (y - yc) - radius * radius;
}

int devide_the_cell(struct box b, Real xc, Real yc, Real radius)
{
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

Real equation_analytical_curvature(Real x, Real y, Real xc, Real yc)
{
  Real a, b, c, radius, A, B, phi, theta;
  a = init_VOF_coefs[0];
  b = init_VOF_coefs[1];
  c = init_VOF_coefs[2];

  x -= xc;
  y -= yc;

  dfetch("radius", &radius);

  phi = 0.5 * atan2(b, (a - c));
  // theta = atan2(y - yc, x - xc);
  theta = atan2(y, x);

  A = radius / (sqrt((pow(a * cos(phi), 2) + b * cos(phi) * sin(phi) + c * pow(sin(phi), 2))) + 1e-50);
  B = radius / (sqrt((pow(a * sin(phi), 2) - b * cos(phi) * sin(phi) + c * pow(cos(phi), 2))) + 1e-50);

  return (A * B) / (pow(A * A * sin(theta - phi) * sin(theta - phi) + B * B * cos(theta - phi) * cos(theta - phi), 1.5) + 1e-50);
}

void _refineFTT_ellipse(void)
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
  nShock = 2 * Ly / dy;
  nShock = 0;

  for (iCell = 0; iCell < numberOfCells; iCell++)
  {
    iOct = iCell / cellNumberInOct;
    cLv = octLv[iOct];
    if (cLv >= maxLevel)
      return;

    dx = dxCell[cLv];
    dy = dyCell[cLv];
    rectangle.pt1.x = xCell[iCell];
    rectangle.pt1.y = yCell[iCell];
    rectangle.pt2.x = xCell[iCell] + dx;
    rectangle.pt2.y = yCell[iCell] + dy;
    /* point within cell, split the cell */

    // if(cLv < 3){
    //   if(cellChOct[iCell] == 0 && cellType[iCell] == 0) splitCell(iCell);
    //   continue;
    // }

    if (devide_the_cell(rectangle, xc, yc, radius) == 1)
    {
      // split cell
      // printf("%f/n", equation_analytical_curvature(xCell[iCell] + dx/2, yCell[iCell] + dy/2, xc, yc));
      printf("in the intermidiate cell \n");
      if (cellChOct[iCell] == 0 && cellType[iCell] == 0)
        splitCell(iCell);
    }
  }
  return;
}

int ptListIntersectBox_extended(struct point *pList, int nPoint, struct box b, int *index)
{
  int ipt, true;
  true = 0;
  for (ipt = 0; ipt < nPoint; ipt++)
  {
    if (ptWithinBox(pList[ipt], b))
    {
      true ++;
      *index = ipt;
      return 1;
    }
  }
  return true;
}

void refineFTT(void)
{
  int cLv, iOct, iCell, chCell;
  Real xPoint, yPoint, theta, delta;
  Real dx, dy, Lx, Ly, radius, xc, yc, xShock;
  struct point *pList;
  struct box rectangle;
  int n, nPoint, nCircle, nShock;
  FILE *fPoints;

  Real A, B, C, a, b, phi, ave_kappa, refine_th;
  int index;

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
  nShock = 2 * Ly / dy;
  nShock = 0;
  nCircle = 4 * 4 * 3.14159 * radius / dx;
  delta = 2 * 3.14159 / nCircle;
  nPoint = nCircle + nShock;

  A = init_VOF_coefs[0];
  B = init_VOF_coefs[1];
  C = init_VOF_coefs[2];

  refine_th = 0.5;
  dfetch("refine_threshold", &refine_th);

  phi = 0.5 * atan2(B, A - C);
  a = radius / sqrt(A * cos(phi) * cos(phi) + B * cos(phi) * sin(phi) + C * sin(phi) * sin(phi));
  b = radius / sqrt(A * sin(phi) * sin(phi) - B * cos(phi) * sin(phi) + C * cos(phi) * cos(phi));

  // read a list of points
  pList = (struct point *)malloc(nPoint * sizeof(struct point));
  for (n = 0; n < nCircle; n++)
  {
    theta = (n + .5) * delta;
    pList[n].x = xc + a * cos(theta) * cos(phi) - b * sin(theta) * sin(phi);
    pList[n].y = yc + a * cos(theta) * sin(phi) + b * sin(theta) * cos(phi);
  }
  /*  delta = Ly/nShock;
    for(n=nCircle; n<nPoint; n++)
    {
      pList[n].x = xShock; pList[n].y = (n+.5-nCircle)*delta;
    }
  */
  fPoints = fopen("meshSeeds", "w");
  for (n = 0; n < nCircle; n++)
  {
    fprintf(fPoints, "%g %g %g\n", pList[n].x, pList[n].y, equation_analytical_curvature(pList[n].x, pList[n].y, xc, yc));
  }
  fprintf(fPoints, "%g %g %g\n", pList[n].x, pList[n].y, equation_analytical_curvature(pList[n].x, pList[n].y, xc, yc));
  fprintf(fPoints, "\n");
  fclose(fPoints);
  /*
  for(n=nCircle; n<nPoint; n++)
  {
    fprintf(fPoints, "%g %g\n", pList[n].x, pList[n].y);
  }
  */
  for (iCell = 0; iCell < numberOfCells; iCell++)
  {
    //    printf("split cell %d, numberOfCells %d\n", iCell, numberOfCells);
    //    getchar();

    iOct = iCell / cellNumberInOct;
    cLv = octLv[iOct];
    if (cLv >= maxLevel)
      return;

    if (cLv < minLevel)
    {
      if (cellChOct[iCell] == 0 && cellType[iCell] == 0)
        splitCell(iCell);
    }

    dx = dxCell[cLv];
    dy = dyCell[cLv];
    rectangle.pt1.x = xCell[iCell];
    rectangle.pt1.y = yCell[iCell];
    rectangle.pt2.x = xCell[iCell] + dx;
    rectangle.pt2.y = yCell[iCell] + dy;
    /* point within cell, split the cell */
    // if(ptListIntersectBox(pList, nPoint, rectangle) > 1)
    index = -1;
    if (ptListIntersectBox_extended(pList, nPoint, rectangle, &index))
    {
      // split cell
      // printf("kappa : %f\n", equation_analytical_curvature(xCell[iCell] + dx/2, yCell[iCell] + dy/2, xc, yc));
      // if(cellChOct[iCell] == 0 && cellType[iCell] == 0) splitCell(iCell);
      // if (cLv < minLevel)
      // {
      //   if (cellChOct[iCell] == 0 && cellType[iCell] == 0)
      //     splitCell(iCell);
      // }
      // else
      {
        // uncomment the following to initialise to curvature
        ave_kappa = equation_analytical_curvature(pList[index].x, pList[index].y, xc, yc);
        if (log(ave_kappa + 1) > refine_th * cLv)
        {
          if (cellChOct[iCell] == 0 && cellType[iCell] == 0)
            splitCell(iCell);
        }
      }
    }
  }
  return;
}
