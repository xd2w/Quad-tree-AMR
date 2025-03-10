#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "ftt.h"
#include "nrutil.h"
#include "pfplib.h"

void printcc3(Real cc[3][3])
{
    int n = 3;
    int i, j;
    printf("\n\n");
    for (j = 0; j < n; j++)
    {
        for (i = 0; i < n; i++)
        {
            printf("\t%f", cc[i][n - 1 - j]);
        }
        printf("\n");
    }
    printf("\n\n");
}

void printcc6(Real cc[6][6])
{
    int n = 6;
    int i, j;
    printf("\n\n");
    for (j = 0; j < n; j++)
    {
        for (i = 0; i < n; i++)
        {
            printf("\t%f", cc[i][n - 1 - j]);
        }
        printf("\n");
    }
    printf("\n\n");
}

void checkSplitVOF(Real cc6[6][6], int iCell)
{
    int i, j, u, v;
    Real refcc[3][3], cc3[3][3];
    Real tol = 1e-5;
    int flag = 0;
    getCellNgbVOF(octPrCell[iCell / 4], refcc);
    for (i = 0; i < 3; i++)
    {
        for (j = 0; j < 3; j++)
        {
            cc3[i][j] = 0.25 * (cc6[(2 * i)][(2 * j)] + cc6[(2 * i) + 1][(2 * j)] + cc6[(2 * i)][(2 * j) + 1] + cc6[(2 * i) + 1][(2 * j) + 1]);
            if (fabs(cc3[i][j] - refcc[i][j]) > tol)
            {
                flag = 1;
            }
        }
    }

    if (flag)
    {
        printf("wrong VOF splitting \n");
        printf("at iCell=%d \n", iCell);

        printf("cc6 : \n");
        printcc6(cc6);

        printf("summed to : \n");
        printcc3(cc3);

        printf("ref : \n");
        printcc3(refcc);
        // exit(1);
    }
    // printf("%d fine\n", iCell);
}

void plotCurvatureAtLevel(int ndata, int level)
{
    Real fraction, kappa, theta;

    Real cc[6][6];

    int ip[] = {0, 1, 0, 1};
    int jp[] = {0, 0, 1, 1};

    int iLv, iCell;
    char fname[] = "DATA/kappa.000";

    fname[13] = '0' + ndata % 10;
    fname[12] = '0' + (ndata / 10) % 10;
    fname[11] = '0' + (ndata / 100) % 10;

    FILE *fp = fopen(fname, "w");

    Real xc, yc;
    xc = 0.5;
    yc = 0.5;

    dfetch("xc", &xc);
    dfetch("yc", &yc);

    for (iCell = 0; iCell < numberOfCells; iCell += 4)
    {
        // if (cellChOct[iCell] == 0)
        if (octLv[iCell / 4] == level)
        {
            fraction = vof[cellChOct[iCell / 4]];
            if (fraction > 0.0 && fraction < 1.0)
            {
                // calc stuff
                iLv = octLv[iCell / 4];
                getCellNgbVOF_6x6(iCell / 4, cc);
                for (int k = 0; k < 4; k++)
                {
                    fraction = vof[iCell + k];
                    if (fraction > 0.0 && fraction < 1.0)
                    {
                        // kappa = curvature_5x5(cc, ip[k], jp[k], dxCell[iLv], dyCell[iLv]);
                        kappa = kappaBarickALELike(iCell, cc);
                        theta = atan2(yCell[iCell + k] - yc, xCell[iCell + k] - xc);
                        fprintf(fp, "%f %f \n", theta, kappa);
                    }
                }
            }
        }
    }
    fprintf(fp, "\n");
    fclose(fp);
}

void plotCurvatureAtLeafCells(int ndata)
{
    Real fraction, kappa, theta, th_kappa;

    Real cc[6][6], cc3[3][3];

    int ip[] = {0, 1, 0, 1};
    int jp[] = {0, 0, 1, 1};

    int iLv, iCell;
    char fname[] = "DATA/kappa.000";

    fname[13] = '0' + ndata % 10;
    fname[12] = '0' + (ndata / 10) % 10;
    fname[11] = '0' + (ndata / 100) % 10;

    FILE *fp = fopen(fname, "w");

    Real xc, yc;
    xc = 0.5;
    yc = 0.5;

    dfetch("xc", &xc);
    dfetch("yc", &yc);

    for (iCell = 0; iCell < numberOfCells; iCell += 4)
    {
        if (cellChOct[iCell] == 0)
        {
            // fraction = vof[cellChOct[iCell / 4]];
            // if (fraction > 0.0 && fraction < 1.0)
            {
                iLv = octLv[iCell / 4];
                getCellNgbVOF_6x6(iCell / 4, cc);
                // checkSplitVOF(cc, iCell);

                for (int k = 0; k < 4; k++)
                {
                    fraction = vof[iCell + k];
                    if (fraction > 0.0 && fraction < 1.0)
                    {
                        // kappa = curvature_5x5(cc, ip[k], jp[k], dxCell[iLv], dyCell[iLv]);
                        kappa = kappaBarickALELike(iCell + k, cc);
                        // kappa = kappaBarickALELike_wider(iCell + k, cc) / 10;
                        // printf("kappa = %g\n", kappa);
                        // exit(0);
                        // kappa = kappaMeier(iCell + k);
                        theta = atan2(yCell[iCell + k] - yc, xCell[iCell + k] - xc + 1e-50);
                        // th_kappa = 1 / (sqrt(pow(yCell[iCell + k] - yc, 2) + pow(xCell[iCell + k] - xc, 2)) + 1e-50);
                        th_kappa = equation_analytical_curvature(xCell[iCell + k], yCell[iCell + k], xc, yc);
                        fprintf(fp, "%f %f %f \n", theta, kappa, th_kappa);
                    }
                }
            }
        }
    }
    fprintf(fp, "\n");
    fclose(fp);
}