#include <stdio.h>
#include <stdlib.h>
#include <math.h>
// #include "pfplib.h"
#include "ftt.h"
#include "nrutil.h"
#include "box.h"
#include "pfplib.h"

#define pi 3.14159265

void computePotential(int iCell, int itNb)
{
    int iOct, cLv, i, count;
    Real xc, yc, radius;
    Real dx, dy, xcoor[5], ycoor[5];
    Real xList[8], yList[8], color[5], lambda, fraction;
    Real x, y, circulation, Lx, Ly;

    xc = 0.75;
    yc = 0.5;
    radius = 0.1;
    circulation = 1.0;
    dfetch("xc", &xc);
    dfetch("yc", &yc);
    dfetch("circulation", &circulation);
    dfetch("radius", &radius);
    dfetch("Lx", &Lx);
    dfetch("Ly", &Ly);
    /* position of the circle at itNb interations */ //???
    // xc += .0075 * itNb;
    // yc += .005 * itNb;

    // printf("\ncircle center at (%g,%g), radius = %g itNb %d\n", xc, yc, radius, itNb);

    iOct = iCell / cellNumberInOct;
    cLv = octLv[iOct];

    dx = dxCell[cLv];
    dy = dyCell[cLv];
    // printf("dx %g dy %g\n", dx, dy);

    x = xCell[iCell] + dx * 0.5;
    y = yCell[iCell] + dy * 0.5;

    vx[iCell] = computeVX(x, y);
    vy[iCell] = computeVY(x, y);

    // for image method to create a rough boundary around the domain
    // vx[iCell] += computeVX(x, y - 2 * Ly, circulation);
    // vx[iCell] += computeVX(x, y + 2 * Ly, circulation);
    // vx[iCell] += computeVX(x - 2 * Lx, y, circulation);
    // vx[iCell] += computeVX(x + 2 * Lx, y, circulation);

    // vy[iCell] += computeVY(x, y - 2 * Ly, circulation);
    // vy[iCell] += computeVY(x, y + 2 * Ly, circulation);
    // vy[iCell] += computeVY(x - 2 * Lx, y, circulation);
    // vy[iCell] += computeVY(x + 2 * Lx, y, circulation);

    // printf("calculating potential at (%g, %g), with (vx, vy) = (%f, %f)\n", x, y, vx[iCell], vy[iCell]);
}

Real computeVX(Real x, Real y)
{
    return 2 * sin(pi * x) * sin(pi * x) * sin(pi * y) * cos(pi * y);
}

Real computeVY(Real x, Real y)
{
    return 2 * cos(pi * x) * sin(pi * x) * cos(pi * y) * cos(pi * y);
}

Real polyArea(Real1D xList, Real1D yList, int count)
{
    int i;
    Real x0, x1, x2;
    Real y0, y1, y2, area;
    area = 0.0;
    x0 = xList[0];
    y0 = yList[0];
    //  printf("count %d \n", count);
    for (i = 1; i < count - 1; i++)
    {
        //     printf("i %d\n", i);
        x1 = xList[i];
        y1 = yList[i];
        x2 = xList[i + 1];
        y2 = yList[i + 1];
        area += 0.5 * ((x1 - x0) * (y2 - y0) + (y1 - y0) * (x0 - x2));
    }
    return area;
}