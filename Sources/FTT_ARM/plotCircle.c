#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "ftt.h"

void plotCircle(int ndata)
{
    int iCell, iOct, iLv, i, i1, i2, i3, index;
    Real fraction, left, right, top, bottom, x, y;
    char fname[] = "DATA/intf.000";
    FILE *finterf;

    i = ndata;
    i1 = i % 10;
    i /= 10;
    i2 = i % 10;
    i /= 10;
    i3 = i % 10;
    fname[10] = '0' + i3;
    fname[11] = '0' + i2;
    fname[12] = '0' + i1;

    finterf = fopen(fname, "w");
    //  iCell = 339;
    //  plotCellInterf(iCell, finterf);
    for (int i = 0; i < numberOfCirclePoints; i++)
    {
        fprintf(finterf, "%f %f \n", xCircle[i], yCircle[i]);
    }
    fprintf(finterf, "%f %f \n", xCircle[0], yCircle[0]);

    fclose(finterf);
    return;
}