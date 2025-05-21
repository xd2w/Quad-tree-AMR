#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "uniform.h"
#include "nrutil.h"

int main()
{
    Lx = 1.0;
    Ly = 1.0;
    nx = 128;
    ny = 128;

    printf("init uniform mesh\n");
    init_uniform_mesh();
    printf("init uniform mesh done\n");
    plot(0);
    plot_vof(0);

    Real cfl = 0.7;
    Real Vmax = 1.0;
    dt = cfl * fmin(dx, dy) / Vmax;

    t_total = 0.0;

    int t_max = 800;
    int t_plot = 8;

    // plic();

    // plot(1);
    // plot_vof(1);

    for (int t = 1; t < t_max; t++)
    {
        // set_IC();
        // printf("set IC done\n");
        solveV();
        // printf("solveV done\n");
        plic();
        t_total += dt;
        if (t % t_plot == 0)
        {
            printf("t = \t%g\n", t_total);
            plot(t / t_plot);
            plot_vof(t / t_plot);
        }
    }
    return 0;
}