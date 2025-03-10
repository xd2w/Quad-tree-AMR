#include <stdio.h>
#include <stdlib.h>
#include "ftt.h"
#include "nrutil.h"
#include <stdbool.h>
#include <math.h>

void plotVOFCross(FILE *fp, int iCell)
{
    int cLv;
    Real xLeft, yLeft, dx, dy;
    cLv = octLv[iCell / 4];
    xLeft = xCell[iCell];
    yLeft = yCell[iCell];
    dx = dxCell[cLv];
    dy = dyCell[cLv];

    fprintf(fp, "%g %g %g\n", xLeft, yLeft, vof[iCell]);
    fprintf(fp, "%g %g %g\n\n", xLeft + 2.0 * dx, yLeft + 2.0 * dy, vof[iCell]);
    fprintf(fp, "%g %g %g\n", xLeft + 2.0 * dx, yLeft, vof[iCell]);
    fprintf(fp, "%g %g %g\n\n", xLeft, yLeft + 2.0 * dy, vof[iCell]);
}

// plot VoF at the leaf cells
void plotVOF(int ndata)
{
    int iCell, i, i1, i2, i3;
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
    //  iCell = 339;
    //  plotCellInterf(iCell, finterf);
    for (iCell = 0; iCell < numberOfCells; iCell++)
    {
        if (cellChOct[iCell] == 0)
        {
            plotVOFCross(fp, iCell);
        }
    }

    fclose(fp);
    return;
}

void plotVOFAtLevel(int ndata, int level)
{
    int iOct, iCell, i, i1, i2, i3;
    char fname[] = "DATA/vof_lev0.000";
    FILE *fp;

    i = ndata;
    i1 = i % 10;
    i /= 10;
    i2 = i % 10;
    i /= 10;
    i3 = i % 10;
    fname[12] = '0' + level;
    fname[14] = '0' + i3;
    fname[15] = '0' + i2;
    fname[16] = '0' + i1;

    fp = fopen(fname, "w");
    //  iCell = 339;
    //  plotCellInterf(iCell, finterf);
    for (iOct = 0; iOct < numberOfOcts; iOct++)
    {
        if (octLv[iOct] == level)
        {
            for (i = 0; i < cellNumberInOct; i++)
            {
                iCell = cellNumberInOct * iOct + i;
                plotVOFCross(fp, iCell);
            }
        }
    }

    fclose(fp);
    return;
}