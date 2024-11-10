#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "ftt.h"
#include "nrutil.h"
#include "box.h"
#include "pfplib.h"

// calculated VOF of flagged cells in both x and y dirs
void stream(void)
{
    printf("streaming\n");
    initPotential(0);
    // propagateCircles(0.001);
    propagateCirclesTruePosition(0.01);
    // flagInterfCells();

    // flagInterfCells();
    // propagateFlag(0);
    // propagateFlag(0);
    // propagateFlag(1);
    // propagateFlag(1);
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
    float x, y, gam, vx, vy;

    gam = 5;
    ffetch("circulation", &gam);
    ffetch("dt", &dt);

    for (int index = 0; index < numberOfCirclePoints; index++)
    {
        x = xCircle[index];
        y = yCircle[index];

        vx = computeVX(x, y, gam);
        vy = computeVY(x, y, gam);

        // for image method to create a rough boundary around the domain
        vx += computeVX(x, y - 2 * Ly, gam);
        vx += computeVX(x, y + 2 * Ly, gam);
        vx += computeVX(x - 2 * Lx, y, gam);
        vx += computeVX(x + 2 * Lx, y, gam);

        vy += computeVY(x, y - 2 * Ly, gam);
        vy += computeVY(x, y + 2 * Ly, gam);
        vy += computeVY(x - 2 * Lx, y, gam);
        vy += computeVY(x + 2 * Lx, y, gam);

        xCircle[index] += dt * vx;
        yCircle[index] += dt * vy;
    }
}
// void propagateCirclesTruePosition(float dt)
// {
//     float x, y, xp, yp, gam, vx, vy;
//     int queue[maxNumberOfCirclePoints];
//     int queueStart = 0;
//     int queueEnd = 0;

//     gam = 5;
//     ffetch("circulation", &gam);
//     ffetch("dt", &dt);

//     for (int index = 0; index + 1 < numberOfCirclePoints + (queueEnd - queueStart); index++)
//     {
//         x = xCircle[index];
//         y = yCircle[index];

//         xp = xCircle[index + 1];
//         yp = yCircle[index + 1];

//         if ((x - xp) * (x - xp) + (y - xp) * (y - yp) > sqrt(2.0 * Ly * Lx) / (1 << (maxLevel)))
//         {
//             x, y = cubicSpine()
//         }

//         vx = computeVX(x, y, gam);
//         vy = computeVY(x, y, gam);

//         // for image method to create a rough boundary around the domain
//         vx += computeVX(x, y - 2 * Ly, gam);
//         vx += computeVX(x, y + 2 * Ly, gam);
//         vx += computeVX(x - 2 * Lx, y, gam);
//         vx += computeVX(x + 2 * Lx, y, gam);

//         vy += computeVY(x, y - 2 * Ly, gam);
//         vy += computeVY(x, y + 2 * Ly, gam);
//         vy += computeVY(x - 2 * Lx, y, gam);
//         vy += computeVY(x + 2 * Lx, y, gam);

//         xCircle[index] += dt * vx;
//         yCircle[index] += dt * vy;
//     }
// }

// copy from one array to another
void copyCellInt1D(Int1D from, Int1D to)
{
    int iCell;
    for (iCell = 0; iCell < numberOfCells; iCell++)
    {
        to[iCell] = from[iCell];
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

// Real cubicSpine(Real x0, Real y0, Real x1, Real y1, Real x2, Real y2, Real x3, Real y3)
// {

// }