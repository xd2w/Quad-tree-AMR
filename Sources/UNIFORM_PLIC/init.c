#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "uniform.h"
#include "nrutil.h"

#include "checkMemory.h"

Real f(Real xc, Real yc, Real r, Real x, Real y)
{
    return (x - xc) * (x - xc) + (y - yc) * (y - yc) - r * r;
}

Real polyArea(Real1D xList, Real1D yList, int count)
{
    int i;
    Real x0, x1, x2;
    Real y0, y1, y2, area;
    area = 0.0;
    x0 = xList[0];
    y0 = yList[0];

    for (i = 1; i < count - 1; i++)
    {
        x1 = xList[i];
        y1 = yList[i];
        x2 = xList[i + 1];
        y2 = yList[i + 1];
        area += 0.5 * ((x1 - x0) * (y2 - y0) + (y1 - y0) * (x0 - x2));
    }
    return area;
}

void init_vof(Real xc, Real yc, Real r)
{
    printf("initing vof\n");
    int i, j, count;
    Real cc[5];
    Real xList[5], yList[5];
    Real xcoor[5], ycoor[5];
    Real lambda;

    for (i = 0; i < nx; i++)
    {
        for (j = 0; j < ny; j++)
        {
            // printf("i %d j %d\n", i, j);
            temp_vof[i][j] = 0.0;

            xcoor[0] = x[i];
            xcoor[1] = x[i + 1];
            xcoor[2] = x[i + 1];
            xcoor[3] = x[i];
            xcoor[4] = x[i];

            ycoor[0] = y[j];
            ycoor[1] = y[j];
            ycoor[2] = y[j + 1];
            ycoor[3] = y[j + 1];
            ycoor[4] = y[j];

            cc[0] = f(xc, yc, r, x[i], y[j]);
            cc[1] = f(xc, yc, r, x[i + 1], y[j]);
            cc[2] = f(xc, yc, r, x[i + 1], y[j + 1]);
            cc[3] = f(xc, yc, r, x[i], y[j + 1]);
            cc[4] = cc[0];

            if (cc[0] < 0)
                if (cc[1] < 0)
                    if (cc[2] < 0)
                        if (cc[3] < 0)
                        {
                            vof[i][j] = 1.0;
                            continue;
                        }

            if (cc[0] > 0)
                if (cc[1] > 0)
                    if (cc[2] > 0)
                        if (cc[3] > 0)
                        {
                            vof[i][j] = 0.0;
                            continue;
                        }

            count = 0;
            for (int k = 0; k < 4; k++)
            {
                if (cc[k] <= 0.0)
                {
                    xList[count] = xcoor[k];
                    yList[count] = ycoor[k];
                    count++;
                }
                if (cc[k] * cc[k + 1] < 0.0)
                {
                    lambda = fabs(cc[k]) / (fabs(cc[k]) + fabs(cc[k + 1]) + 1e-50);
                    xList[count] = (1.0 - lambda) * xcoor[k] + lambda * xcoor[k + 1];
                    yList[count] = (1.0 - lambda) * ycoor[k] + lambda * ycoor[k + 1];
                    count++;
                }
            }

            if (count)
            {
                xList[count] = xList[0];
                yList[count] = yList[0];
            }

            vof[i][j] = polyArea(xList, yList, count) / (dx * dy);
        }
    }
}

void init_uniform_mesh(void)
{
    vof = dmatrix(0, nx, 0, ny);
    temp_vof = dmatrix(0, nx, 0, ny);
    ux = dmatrix(0, nx + 1, 0, ny + 1);
    uy = dmatrix(0, nx + 1, 0, ny + 1);

    x = dvector(0, nx + 1);
    y = dvector(0, ny + 1);

    dx = Lx / (double)nx;
    dy = Ly / (double)ny;

    for (int i = 0; i <= nx; i++)
        x[i] = i * dx;

    for (int j = 0; j <= ny; j++)
        y[j] = j * dy;

    printf("init vof\n");
    init_vof(0.5 * Lx, 0.5 * Ly, 0.2 * Lx);
    init_solver_vars();
}