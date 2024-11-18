#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "ftt.h"
#include "nrutil.h"
#include "box.h"
#include "pfplib.h"
#include "interpolation.h"

void copyArray(Real1D from, double *to, int length);
void propagateCirclesTruePosition(float dt);

// calculated VOF of flagged cells in both x and y dirs
void stream(void)
{
    // printf("streaming\n");
    // initPotential(0);
    // propagateCircles(0.001);
    propagateCirclesTruePosition(0.01);
}

// flag cells at the interface of 2 fluids
// void flagInterfCells(void)
// {
//     printf("Flagging cells\n");
//     int iCell, iOct, iLv, index;
//     Real fraction, left, right, top, bottom, x, y;
//     for (iCell = 0; iCell < numberOfCells; iCell++)
//     {
//         cellFlag[iCell] = 0;
//         iOct = iCell / cellNumberInOct;
//         iLv = octLv[iOct];
//         if (iLv == maxLevel)
//         {
//             left = xCell[iCell];
//             bottom = yCell[iCell];
//             right = left + dxCell[iLv];
//             top = bottom + dyCell[iLv];

//             // printf("%f, %f, %f, %f\n", left, right, top, bottom);
//             // printf("%f, %f\n", dxCell[iLv], dyCell[iLv]);
//             // exit(1);

//             for (index = 0; index < numberOfCirclePoints; index++)
//             {
//                 x = xCircle[index];
//                 y = yCircle[index];

//                 if (((left <= x) && (x < right)) && ((bottom <= y) && (y < top)))
//                 {
//                     cellFlag[iCell] = 1;
//                     printf("cell %i flagged\n", iCell);
//                 }
//             }
//         }
//     }
// }

void flagInterfCells(void)
{
    printf("Flagging cells\n");
    int iCell, iOct, iLv, index;
    Real fraction, left, right, top, bottom, x, y;
    for (iCell = 0; iCell < numberOfCells; iCell++)
    {
        cellFlag[iCell] = 0;
        iOct = iCell / cellNumberInOct;
        iLv = octLv[iOct];
        // if (iLv == maxLevel)
        // {
        left = xCell[iCell];
        bottom = yCell[iCell];
        right = left + dxCell[iLv];
        top = bottom + dyCell[iLv];

        // printf("%f, %f, %f, %f\n", left, right, top, bottom);
        // printf("%f, %f\n", dxCell[iLv], dyCell[iLv]);
        // exit(1);

        for (index = 0; index < numberOfCirclePoints; index++)
        {
            x = xCircle[index];
            y = yCircle[index];

            if (((left <= x) && (x < right)) && ((bottom <= y) && (y < top)))
            {
                cellFlag[iCell] = 1;
                // printf("cell %i flagged\n", iCell);
            }
        }
        // }
    }
}

void propagateCircles(float dt)
{
    int index, iCell, iOct, iLv;
    float x, y, left, right, top, bottom;
    for (iCell = 0; iCell < numberOfCells; iCell++)
    {

        iOct = iCell / cellNumberInOct;
        iLv = octLv[iOct];

        left = xCell[iCell];
        bottom = yCell[iCell];
        right = left + dxCell[iLv];
        top = bottom + dyCell[iLv];

        for (index = 0; index < numberOfCirclePoints; index++)
        {
            x = xCircle[index];
            y = yCircle[index];

            if (((left < x) && (x < right)) && ((bottom < y) && (y < top)))
            {
                cellFlag[iCell] = 1;
                xCircle[index] += dt * vx[iCell];
                yCircle[index] += dt * vy[iCell];
            }
        }
    }
}

void propagateCirclesTruePosition(float dt)
{
    double x, y, vx, vy, xprev, yprev, newx, newy;
    ffetch("dt", &dt);
    double *bufferx = dvector(1, numberOfCirclePoints);
    double *buffery = dvector(1, numberOfCirclePoints);
    double *y2 = dvector(1, numberOfCirclePoints);

    int deg = 4;

    double *xslice = dvector(1, deg);
    double *yslice = dvector(1, deg);

    int og_length = numberOfCirclePoints;
    float tol = 5 * sqrt(2) * sqrt(Lx * Ly) * pow(2.0, -maxLevel);
    double dydx = (yCircle[og_length - 1] - yCircle[0]) / (xCircle[og_length - 1] - xCircle[0]);

    copyArray(xCircle, bufferx, og_length);
    copyArray(yCircle, buffery, og_length);

    // printf("spline started\n");
    spline(buffery, buffery, og_length, dydx, dydx, y2);
    printf("spline is done\n");

    xprev = xCircle[og_length - 1];
    yprev = yCircle[og_length - 1];

    int counter = 0;
    for (int index = 0; index < og_length; index++)
    {
        x = bufferx[index];
        y = buffery[index];
        if (numberOfCirclePoints < maxNumberOfCirclePoints)
        {
            // printf("so for so good");
            if (sqrt((x - xprev) * (x - xprev) + (y - yprev) * (y - yprev)) > tol)
            {
                for (int k = 0; k < deg; k++)
                {
                    xslice[k] = bufferx[(og_length - 2 + k + index) % og_length];
                    xslice[k] = bufferx[(og_length - 2 + k + index) % og_length];
                }
                // printf("%d, (%d) : %f \n", index, (og_length + index - 1) % og_length, sqrt((x - xprev) * (x - xprev) + (y - yprev) * (y - yprev)));
                newx = 0.5 * (x + xprev);
                // splint(xslice, yslice, y2, og_length, newx, &newy); // interpolates y for given x
                newy = 0.5 * (y + yprev);
                printf("(%f %f) (%f %f) (%f %f)\n", x, y, newx, newy, xprev, yprev);

                vx = computeVX(newx, newy);
                vy = computeVY(newx, newy);

                xCircle[counter] = newx + (dt * vx);
                yCircle[counter] = newy + (dt * vy);
                counter++;
                numberOfCirclePoints++;
            }
        }

        vx = computeVX(x, y);
        vy = computeVY(x, y);

        xCircle[counter] = x + (dt * vx);
        yCircle[counter] = y + (dt * vy);
        counter++;

        xprev = x;
        yprev = y;
    }
    // free_dvector(bufferx, 1, maxNumberOfCirclePoints);
    // free_dvector(buffery, 1, maxNumberOfCirclePoints);
    // free_dvector(y2, 1, maxNumberOfCirclePoints);
}

// copy from one array to another
void copyCellInt1D(Int1D from, Int1D to)
{
    int iCell;
    for (iCell = 0; iCell < numberOfCells; iCell++)
    {
        to[iCell] = from[iCell];
    }
}

// copy from one array to another
void copyArray(Real1D from, double *to, int length)
{
    int i;
    for (i = 0; i < length; i++)
    {
        to[i] = from[i];
    }
}

/* propagate flag in direction dir
(flags the neighboring cells of the flagged cells in direction x or y)
*/
void propagateFlag(int dir)
{
    int iCell, ngbCell;

    copyCellInt1D(cellFlag, cellMark);
    for (iCell = 0; iCell < numberOfCells; iCell++)
    {
        if (cellMark[iCell]) // cellMark is Copy of cellFlag
        {
            if (dir == 0)
            {
                // left 0  right 1
                ngbCell = cellNb[0][iCell];
                cellFlag[ngbCell] = 1;
                ngbCell = cellNb[1][iCell];
                cellFlag[ngbCell] = 1;
            }
            else
            {
                // down 2  up 3
                ngbCell = cellNb[2][iCell];
                cellFlag[ngbCell] = 1;
                ngbCell = cellNb[3][iCell];
                cellFlag[ngbCell] = 1;
            }
        }
    }
}