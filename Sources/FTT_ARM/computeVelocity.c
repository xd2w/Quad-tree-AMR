#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "pfplib.h"
#include "ftt.h"
#include "nrutil.h"
#include "box.h"
#include "pfplib.h"

#define pi 3.14159265

// void computeVelocity(void)
// {
// }

Real computeVX(Real x, Real y)
{
    // return x;
    // return 0.5 * y;
    if (runtime == 0)
        return -2 * sin(pi * x) * sin(pi * x) * sin(pi * (y + 0.5)) * cos(pi * (y + 0.5));

    if (t_total < runtime / 2)
    {
        return -2 * sin(pi * x) * sin(pi * x) * sin(pi * (y + 0.5)) * cos(pi * (y + 0.5));
    }
    if (t_total < runtime)
    {
        return 2 * sin(pi * x) * sin(pi * x) * sin(pi * (y + 0.5)) * cos(pi * (y + 0.5));
    }
    return 0;
}

Real computeVY(Real x, Real y)
{
    // return 0;
    // return 0.25 * x;
    if (runtime == 0)
        return -2 * cos(pi * x) * sin(pi * x) * cos(pi * (y + 0.5)) * cos(pi * (y + 0.5));

    if (t_total < runtime / 2)
    {
        return -2 * cos(pi * x) * sin(pi * x) * cos(pi * (y + 0.5)) * cos(pi * (y + 0.5));
    }
    if (t_total < runtime)
    {
        return 2 * cos(pi * x) * sin(pi * x) * cos(pi * (y + .5)) * cos(pi * (y + .5));
    }
    return 0;
}

void computeVelocityAtLeaves(void)
{
    float px[] = {0, 1, 0, 1};
    float py[] = {0, 0, 1, 1};
    Real dx, dy, x, y;
    int iCell, iLv, k;

    for (iCell = 0; iCell < numberOfCells; iCell++)
    {
        if (cellChOct[iCell] == 0)
        {
            iLv = octLv[iCell / 4];
            dx = dxCell[iLv];
            dy = dyCell[iLv];
            x = xCell[iCell];
            y = yCell[iCell];

            for (k = 0; k < 4; k++)
            {
                u[k][iCell] = computeVX(x + px[k] * dx, y + py[k] * dy);
                v[k][iCell] = computeVY(x + px[k] * dx, y + py[k] * dy);
            }
        }
    }
}
