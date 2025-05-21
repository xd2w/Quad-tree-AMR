#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "uniform.h"
#include "nrutil.h"

void init_solver_vars(void)
{
    ux_star = dmatrix(0, nx, 0, ny);
    uy_star = dmatrix(0, nx, 0, ny);
    phi = dmatrix(0, nx, 0, ny);
}

void solveV(void)
{
    Real duxdx, duydy, duxdy, duydx, temp, delta, maxDelta;
    Real d2uxdx2, d2uxdy2, d2uydy2, d2uydx2;

    Real nu = 0.1; // viscosity

    // solve the possion equation
    for (int i = 1; i < nx - 1; i++)
    {
        for (int j = 1; j < ny - 1; j++)
        {
            duxdx = (ux[i + 1][j] - ux[i - 1][j]) / (2 * dx);
            duxdy = (ux[i][j + 1] - ux[i][j - 1]) / (2 * dy);
            duydx = (uy[i + 1][j] - uy[i - 1][j]) / (2 * dx);
            duydy = (uy[i][j + 1] - uy[i][j - 1]) / (2 * dy);

            d2uxdx2 = (ux[i + 1][j] - 2 * ux[i][j] + ux[i - 1][j]) / (dx * dx);
            d2uxdy2 = (ux[i][j + 1] - 2 * ux[i][j] + ux[i][j - 1]) / (dy * dy);
            d2uydx2 = (uy[i + 1][j] - 2 * uy[i][j] + uy[i - 1][j]) / (dx * dx);
            d2uydy2 = (uy[i][j + 1] - 2 * uy[i][j] + uy[i][j - 1]) / (dy * dy);

            ux_star[i][j] = ux[i][j] - dt * (ux[i][j] * duxdx + uy[i][j] * duydy + nu * (d2uxdx2 + d2uxdy2));
            uy_star[i][j] = uy[i][j] - dt * (ux[i][j] * duydx + uy[i][j] * duydy + nu * (d2uydx2 + d2uydy2));
        }
    }

    // gauss seidel iteration
    Real omega = 0.9; // relaxation factor
    maxDelta = 1;
    while (maxDelta > 1e-6)
    {
        maxDelta = 0;
        for (int i = 1; i < nx - 1; i++)
        {
            for (int j = 1; j < ny - 1; j++)
            {
                duxdx = (ux_star[i + 1][j] - ux_star[i - 1][j]) / (2 * dx);
                duydy = (uy_star[i][j + 1] - uy_star[i][j - 1]) / (2 * dy);

                temp = (phi[i + 1][j] + phi[i - 1][j] + phi[i][j + 1] + phi[i][j - 1]);
                delta = 0.25 * (temp - vof[i][j] * dx * dy * dt * (duxdx + duydy));

                if (fabs(delta) > maxDelta)
                    maxDelta = fabs(delta);

                phi[i][j] = (1 - omega) * phi[i][j] + omega * delta;
            }
        }
    }

    for (int i = 1; i < nx - 1; i++)
    {
        for (int j = 1; j < ny - 1; j++)
        {
            duxdx = (phi[i + 1][j] - phi[i - 1][j]) / (2 * dx);
            duydy = (phi[i][j + 1] - phi[i][j - 1]) / (2 * dy);

            ux[i][j] = ux_star[i][j] - (duxdx + duydy);
            uy[i][j] = uy_star[i][j] - (duxdx + duydy);
        }
    }
}