#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "ftt.h"
#include "pfplib.h"
#include "nrutil.h"

void setTimeStep(void)
{
    int iCell, intrinsicMaxLevel;
    Real intrinsicMaxU2;

    intrinsicMaxLevel = 0;
    intrinsicMaxU2 = 0;

    for (iCell = 0; iCell < numberOfCells; iCell++)
    {
        // if (ux[iCell] * ux[iCell] + uy[iCell] * uy[iCell] > intrinsicMaxU2)
        // {
        //     intrinsicMaxU2 = ux[iCell] * ux[iCell] + uy[iCell] * uy[iCell];
        // }
        if (octLv[iCell / cellNumberInOct] > intrinsicMaxLevel)
        {
            intrinsicMaxLevel = octLv[iCell / cellNumberInOct];
        }
    }

    intrinsicMaxU2 = 1;

    if (intrinsicMaxU2 == 0)
    {
        printf("****************************************************\n");
        printf("Error in miscUtils.c: setTimeStep: invalid ux and uy\n\n");
        exit(1);
    }
    global_dt = CFL * (DMIN(dxCell[intrinsicMaxLevel + 1], dyCell[intrinsicMaxLevel + 1])) / sqrt(intrinsicMaxU2);
}

void show6x6VofGrid(int iCell)
{
    // iCell = 196; // octPrCell[196 / 4];
    printf("iCell = %d\n", iCell);
    printf("vof = %g\n", vof[iCell]);

    for (int k = 0; k < 4; k++)
    {
        printf("cell neighbour[%d] = %d\n", k, cellNb[k][iCell]);
        for (int w = 0; w < 4; w++)
        {
            printf("Nb[%d]'s neighbour[%d] = %d\n", cellNb[k][iCell], w, cellNb[w][cellNb[k][iCell]]);
        }
    }

    for (int k = 0; k < 4; k++)
    {
        printf("cell neighbour[%d] oct size = %d\n", k, octLv[cellNb[k][iCell] / 4]);
    }

    for (int k = 0; k < 4; k++)
    {
        printf("Oct neighbour[%d] = %d\n", k, octNb[k][iCell / 4]);
    }

    Real cc[6][6];
    getCellNgbVOF_6x6(iCell / 4, cc);

    int i, j;
    printf("\n\n");
    for (j = 0; j < 6; j++)
    {
        for (i = 0; i < 6; i++)
        {
            printf("\t%f", cc[i][5 - j]);
        }
        printf("\n");
    }
    printf("\n\n");

    double kappas[4];
    curvature_6x6(cc, kappas, dxCell[octLv[iCell / 4]], dyCell[octLv[iCell / 4]]);

    for (j = 0; j < 4; j++)
    {
        printf("\t%g", kappas[i]);
    }
    printf("\n");
}