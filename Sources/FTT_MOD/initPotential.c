#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "ftt.h"
#include "nrutil.h"
#include "box.h"
#include "pfplib.h"

void initPotential(int itNb)
{
    int iCell, iLv, iOct, iFlag, cellNumber;

    // flag the leave cells
    cellNumber = numberOfCells;
    for (iCell = 0; iCell < cellNumber; iCell++)
    {
        cellFlag[iCell] = 0;
    }
    // printf("Computing potential init _ got this far");
    for (iCell = 0; iCell < cellNumber; iCell++)
    {
        if (cellChOct[iCell] == 0)
        {
            // vof[iCell] = computeVOF(iCell, itNb);
            computePotential(iCell, itNb);
            // printf("Computing potential");
        }
    }
    // restrField();

    return;
}

void initCircle()
{
    Real xc, yc, radius, pi, theta;
    xc = 0.5;
    yc = 0.5;
    radius = 0.2;
    numberOfCirclePoints = 20;

    pi = 3.14159265;
    theta = 0.02;

    dfetch("xc", &xc);
    dfetch("yc", &yc);
    dfetch("radius", &radius);
    ifetch("numberOfCirclePoints", &numberOfCirclePoints);
    if (numberOfCirclePoints > maxNumberOfCirclePoints)
    {
        numberOfCirclePoints = maxNumberOfCirclePoints;
        printf("numberOfCirclePoints too high, capped to max");
    }

    for (int i = 0; i < numberOfCirclePoints; i++)
    {
        xCircle[i] = xc + radius * cos(theta);
        yCircle[i] = yc + radius * sin(theta);
        theta += 2 * pi / (numberOfCirclePoints);
    }
}