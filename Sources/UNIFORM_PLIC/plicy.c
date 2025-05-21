#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "uniform.h"
#include "nrutil.h"
#define MAX(x, y) ((x) > (y) ? (x) : (y))
#define MIN(x, y) ((x) < (y) ? (x) : (y))

void plicy(Real2D from, Real2D to)
{
    int i, j, inv;
    Real mm1, mm2, mx, mz, alpha;
    Real V1, V2, work1, work2, work3;
    Real s1, s2;

    for (i = 1; i < nx - 1; i++)
    {
        for (j = 1; j < ny - 1; j++)
        {

            s1 = 0.5 * (dt / dy) * (uy[i][j] + uy[i + 1][j]);
            s2 = 0.5 * (dt / dy) * (uy[i][j + 1] + uy[i + 1][j + 1]);

            // s1 = 0.5 * (dt / dy) * (computeVY(x[i], y[j]) + computeVY(x[i + 1], y[j]));
            // s2 = 0.5 * (dt / dy) * (computeVY(x[i], y[j + 1]) + computeVY(x[i + 1], y[j + 1]));

            // s1 = -0.1;
            // s2 = -0.1;

            if (from[i][j] == 1)
            {
                work1 = DMAX(-s1, 0.0);
                work2 = 1.0 - DMAX(s1, 0.0) + DMIN(s2, 0.0);
                work3 = DMAX(s2, 0.0);
            }

            if (from[i][j] == 0)
            {
                work1 = 0.0;
                work2 = 0.0;
                work3 = 0.0;
            }

            if (0 < from[i][j] && from[i][j] < 1.0)
            {
                mm1 = from[i - 1][j - 1] + 2. * from[i - 1][j] + from[i - 1][j + 1];
                mm2 = from[i + 1][j - 1] + 2. * from[i + 1][j] + from[i + 1][j + 1];
                mx = mm1 - mm2;
                mm1 = from[i - 1][j - 1] + 2. * from[i][j - 1] + from[i + 1][j - 1];
                mm2 = from[i - 1][j + 1] + 2. * from[i][j + 1] + from[i + 1][j + 1];
                mz = mm1 - mm2;

                inv = 0;
                if (mz < 0.)
                {
                    mm1 = -s1;
                    s1 = -s2;
                    s2 = mm1;
                    mz = -mz;
                    inv = 1;
                }
                mx = fabs(mx) + 1.e-50;
                mz = mz + 1.e-50;
                mm2 = DMAX(mx, mz);
                mx = mx / mm2;
                mz = mz / mm2;

                /* get alpha to determine the equation of the interface */
                mm1 = DMIN(mx, mz);
                V1 = 0.5 * mm1;
                V2 = 1.0 - V1;
                if (from[i][j] <= V1)
                    alpha = sqrt(2.0 * from[i][j] * mm1);
                else if (from[i][j] <= V2)
                    alpha = from[i][j] + 0.5 * mm1;
                else
                    alpha = mm1 + 1.0 - sqrt(2.0 * (1.0 - from[i][j]) * mm1);

                /* the new equation of the interface after advection */
                mz = mz / (1.0 - s1 + s2);
                alpha = alpha + mz * s1;

                /* calculate works */
                mm1 = DMAX(-s1, 0.0);
                mm2 = alpha + mz * mm1;
                work1 = VOL2(mz, mx, mm2, mm1);

                mm1 = 1.0 - DMAX(s1, 0.0) + DMIN(s2, 0.0);
                mm2 = alpha - mz * DMAX(s1, 0.0);
                work2 = VOL2(mz, mx, mm2, mm1);

                mm1 = DMAX(s2, 0.0);
                mm2 = alpha - mz;
                work3 = VOL2(mz, mx, mm2, mm1);

                /* symmetry conditions */
                if (inv == 1)
                {
                    mm1 = work1;
                    work1 = work3;
                    work3 = mm1;
                }
            }
            to[i][j - 1] += work1;
            to[i][j] += work2;
            to[i][j + 1] += work3;
        }
    }
}