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
    global_dt = CFL * 0.5 * (DMIN(dxCell[intrinsicMaxLevel], dyCell[intrinsicMaxLevel])) / sqrt(intrinsicMaxU2);
    printf("current Max Level = %d\n", intrinsicMaxLevel);
}

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