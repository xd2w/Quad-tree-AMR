#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "uniform.h"
#include "nrutil.h"

Real VOL2(Real mx, Real mz, Real alpha, Real b)
{
    Real area, mb, alphaa;

    /* 0<mx or 0<mz */
    mb = mx * b;
    if (alpha <= 0.0 || b <= 0.0)
        area = 0.0;
    else if (alpha >= (mb + mz))
        area = b;
    else
    {
        if (mb <= mz)
        {
            if (alpha < mb)
                area = 0.5 * alpha * alpha / (mx * mz);
            else if (alpha > mz)
            {
                alphaa = mb + mz - alpha;
                area = b - 0.5 * alphaa * alphaa / (mx * mz);
            }
            else
                area = b * (alpha - 0.5 * mb) / mz;
        }
        else
        {
            if (alpha < mz)
                area = 0.5 * alpha * alpha / (mx * mz);
            else if (alpha > mb)
            {
                alphaa = mb + mz - alpha;
                area = b - 0.5 * alphaa * alphaa / (mx * mz);
            }
            else
                area = (alpha - 0.5 * mz) / mx;
        }
    }

    return area;
}

void copy_vof(Real2D from, Real2D to)
{
    int i, j;
    for (i = 0; i < nx; i++)
        for (j = 0; j < ny; j++)
            to[i][j] = from[i][j];
}

void clear_vof(Real2D from)
{
    int i, j;
    for (i = 0; i < nx; i++)
        for (j = 0; j < ny; j++)
            from[i][j] = 0.0;
}

void correct_vof(Real2D from)
{
    int i, j;
    Real tol = 1.0e-16;
    for (i = 0; i < nx; i++)
        for (j = 0; j < ny; j++)
        {
            if (from[i][j] > 1.0 - tol)
                from[i][j] = 1.0;

            else if (from[i][j] < tol)
                from[i][j] = 0.0;
        }
}