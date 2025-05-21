#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "uniform.h"
#include "nrutil.h"

void set_IC(void)
{
    // printf("setting IC\n");
    // for (int i = 0; i < nx + 1; i++)
    // {
    //     for (int j = 0; j < ny + 1; j++)
    //     {
    //         ux[i][j] = 0.0;
    //         uy[i][j] = 0.0;
    //     }
    // }
}

void init_solver_vars(void)
{
}

Real computeVX(Real x, Real y)
{
    // return x;
    // return 0.5 * y;
    if (t_total < 2)
    {
        return -2 * sin(M_PI * x) * sin(M_PI * x) * sin(M_PI * (y + 0.5)) * cos(M_PI * (y + 0.5));
    }
    if (t_total < 4)
    {
        // printf("flipped\n");
        return 2 * sin(M_PI * x) * sin(M_PI * x) * sin(M_PI * (y + 0.5)) * cos(M_PI * (y + 0.5));
    }
    return 0;
}

Real computeVY(Real x, Real y)
{
    // return 0;
    // return 0.25 * x;
    if (t_total < 2)
    {
        return -2 * cos(M_PI * x) * sin(M_PI * x) * cos(M_PI * (y + 0.5)) * cos(M_PI * (y + 0.5));
    }
    if (t_total < 4)
    {
        // printf("flipped\n");
        return 2 * cos(M_PI * x) * sin(M_PI * x) * cos(M_PI * (y + 0.5)) * cos(M_PI * (y + 0.5));
    }
    return 0;
}

void solveV(void)
{
    for (int i = 0; i < nx + 1; i++)
    {
        for (int j = 0; j < ny + 1; j++)
        {
            ux[i][j] = computeVX(x[i], y[j]);
            uy[i][j] = computeVY(x[i], y[j]);
        }
    }
}