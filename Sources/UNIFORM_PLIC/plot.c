#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "uniform.h"
#include "nrutil.h"
#define MAX(x, y) ((x) > (y) ? (x) : (y))
#define MIN(x, y) ((x) < (y) ? (x) : (y))

void plot(int ndata)
{
    // printf("plotting\n");
    int i, j, i1, i2, i3;
    char fname[] = "DATA/intf.000";
    FILE *fp;

    i = ndata;
    i1 = i % 10;
    i /= 10;
    i2 = i % 10;
    i /= 10;
    i3 = i % 10;

    fname[10] = '0' + i3;
    fname[11] = '0' + i2;
    fname[12] = '0' + i1;

    fp = fopen(fname, "w");
    if (fp == NULL)
    {
        printf("Error opening file!\n");
        exit(1);
    }

    for (i = 0; i < nx; i++)
    {
        for (j = 0; j < ny; j++)
        {
            if (0 < vof[i][j] && vof[i][j] < 1.0)
            {
                plot_one(fp, i, j);
                fprintf(fp, "\n");
            }
        }
    }

    fclose(fp);
    // printf("plotting done\n");
}

void plot_vof(int ndata)
{
    // printf("plotting vof\n");
    int i, j, i1, i2, i3;
    char fname[] = "DATA/vof.000";
    FILE *fp;

    i = ndata;
    i1 = i % 10;
    i /= 10;
    i2 = i % 10;
    i /= 10;
    i3 = i % 10;

    fname[9] = '0' + i3;
    fname[10] = '0' + i2;
    fname[11] = '0' + i1;

    fp = fopen(fname, "w");
    if (fp == NULL)
    {
        printf("Error opening file!\n");
        exit(1);
    }

    for (i = 0; i < nx; i++)
    {
        for (j = 0; j < ny; j++)
        {
            fprintf(fp, "%f %f %f\n", x[i], y[j], vof[i][j]);
        }
        fprintf(fp, "\n");
    }

    fclose(fp);
}

void plot_one(FILE *fp, int i, int j)
{
    Real x1, x2, z1, z2;
    Real mm1, mm2, mx, mz;
    Real alpha, V1, V2;
    int invx, invz;

    mm1 = vof[i - 1][j - 1] + 2.0 * vof[i - 1][j] + vof[i - 1][j + 1];
    mm2 = vof[i + 1][j - 1] + 2.0 * vof[i + 1][j] + vof[i + 1][j + 1];
    mx = mm1 - mm2;
    mm1 = vof[i - 1][j - 1] + 2.0 * vof[i][j - 1] + vof[i + 1][j - 1];
    mm2 = vof[i - 1][j + 1] + 2.0 * vof[i][j + 1] + vof[i + 1][j + 1];
    mz = mm1 - mm2;

    // symmetry
    invx = 0;
    if (mx < 0)
    {
        mx = -mx;
        invx = 1;
    }
    invz = 0;
    if (mz < 0)
    {
        mz = -mz;
        invz = 1;
    }
    mx += 1e-50;
    mz += 1e-50;
    mm2 = MAX(mx, mz);
    mx /= mm2;
    mz /= mm2;

    mm1 = MIN(mx, mz);
    V1 = 0.5 * mm1;
    V2 = 1 - V1;
    // determine alpha from volume fraction
    if (vof[i][j] < V1)
    {
        alpha = sqrt(2.0 * vof[i][j] * mm1);
        x1 = alpha / mx;
        z1 = 0.0;
        x2 = 0.0;
        z2 = alpha / mz;
        if (x1 > 1.0)
            x1 = 1.0;
        if (z2 > 1.0)
            z2 = 1.0;
    }
    else if (vof[i][j] < V2)
    {
        alpha = vof[i][j] + 0.5 * mm1;
        if (mx > mz)
        {
            x1 = alpha / mx;
            x2 = (alpha - mz) / mx;
            z1 = 0.0;
            z2 = 1.0;
        }
        else
        {
            x1 = 0.0;
            x2 = 1.0;
            z1 = alpha / mz;
            z2 = (alpha - mx) / mz;
        }
    }
    else
    {
        alpha = sqrt(2.0 * (1.0 - vof[i][j]) * mm1);
        x1 = 1.0 - alpha / mx;
        z1 = 1.0;
        x2 = 1.0;
        z2 = 1.0 - alpha / mz;
        if (x1 < 0.0)
            x1 = 0.0;
        if (z2 < 0.0)
            z2 = 0.0;
    }
    // use symmetry
    if (invx != 0)
    {
        x1 = 1.0 - x1;
        x2 = 1.0 - x2;
    }
    if (invz != 0)
    {
        z1 = 1.0 - z1;
        z2 = 1.0 - z2;
    }
    x1 = x[i] + dx * x1;
    z1 = y[j] + dy * z1;
    x2 = x[i] + dx * x2;
    z2 = y[j] + dy * z2;
    fprintf(fp, "%f %f \n", x1, z1);
    fprintf(fp, "%f %f \n", x2, z2);
    fprintf(fp, "\n");
}