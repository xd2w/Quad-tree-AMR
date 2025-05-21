#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "uniform.h"
#include "nrutil.h"

// extern Real w[9];
// extern Real e_x[9];
// extern Real e_y[9];
// extern Real2D f[9];
// extern Real2D f_2[9];
// extern Real2D feq[9];
// extern Real2D rho;
// extern Real tau;

static double w[9] = {4.0 / 9.0, 1.0 / 9.0, 1.0 / 9.0, 1.0 / 9.0, 1.0 / 9.0, 1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0};
static double e_x[9] = {0, 0, 1, 0, -1, -1, 1, 1, -1};
static double e_y[9] = {0, 1, 0, -1, 0, 1, 1, -1, -1};

void init_solver_vars(void)
{
    // D2Q9 model
    //  5 1 6
    //  4 0 2
    //  8 3 7

    f0 = dmatrix(0, nx, 0, ny);
    f1 = dmatrix(0, nx, 0, ny);
    f2 = dmatrix(0, nx, 0, ny);
    f3 = dmatrix(0, nx, 0, ny);
    f4 = dmatrix(0, nx, 0, ny);
    f5 = dmatrix(0, nx, 0, ny);
    f6 = dmatrix(0, nx, 0, ny);
    f7 = dmatrix(0, nx, 0, ny);
    f8 = dmatrix(0, nx, 0, ny);

    f0_2 = dmatrix(0, nx, 0, ny);
    f1_2 = dmatrix(0, nx, 0, ny);
    f2_2 = dmatrix(0, nx, 0, ny);
    f3_2 = dmatrix(0, nx, 0, ny);
    f4_2 = dmatrix(0, nx, 0, ny);
    f5_2 = dmatrix(0, nx, 0, ny);
    f6_2 = dmatrix(0, nx, 0, ny);
    f7_2 = dmatrix(0, nx, 0, ny);
    f8_2 = dmatrix(0, nx, 0, ny);

    rho = dmatrix(0, nx, 0, ny);

    tau = 0.6; // relaxation time
    printf("init solver vars\n");
}

void set_BC(void)
{
    for (int i = 0; i < nx; i++)
    {
        ux[i][0] = ux[i][1];
        ux[i][ny] = ux[i][ny - 1];
    }
    for (int j = 0; j < ny + 1; j++)
    {
        uy[0][j] = uy[1][j];
        uy[nx][j] = uy[nx - 1][j];
    }
}

void set_IC(void)
{
    int i, j, k;
    for (i = 0; i < nx; i++)
    {
        // f0[i][ny - 1] = 1.0;
        // f1[i][ny - 1] = 1.0;
        // f2[i][ny - 1] = 2.0;
        // f3[i][ny - 1] = 1.0;
        // f4[i][ny - 1] = 1.0;
        // f5[i][ny - 1] = 1.0;
        // f6[i][ny - 1] = 1.0;
        // f7[i][ny - 1] = 1.0;
        // f8[i][ny - 1] = 1.0;
        for (j = 0; j < ny; j++)
        {
            // printf("i %d j %d\n", i, j);
            f0[i][j] = 1.0;
            f1[i][j] = 1.0;
            f2[i][j] = 1.0;
            f3[i][j] = 1.0;
            f4[i][j] = 1.0;
            f5[i][j] = 1.0;
            f6[i][j] = 1.0;
            f7[i][j] = 1.0;
            f8[i][j] = 1.0;
        }
    }
}

void BC(void)
{
    for (int i = 0; i < nx; i++)
    {

        f1_2[i][0] = f3_2[i][0];
        f5_2[i][0] = f7_2[i][0];
        f6_2[i][0] = f8_2[i][0];
    }
}

Real f_eq(int i, int j, int k)
{
    Real eu = e_x[k] * ux[i][j] + e_y[k] * uy[i][j];
    Real uu = ux[i][j] * ux[i][j] + uy[i][j] * uy[i][j];
    return rho[i][j] * w[k] * (1 + 3 * (eu) + 9 * (eu) * (eu) / 2 - 3 * (uu) / 2);
}

Real F_term(int i, int j, int k)
{

    // return w[k] * (1 - (1 / (2 * tau))) * ((e_x[k] - ux[i][j]) + (e_y[k] - uy[i][j]) + (e_x[k] * (e_x[k] - ux[i][j]) + e_y[k] * (e_y[k] - uy[i][j])))*F;
}

void stream_and_BC(void)
{
    // boundary condition
    for (int j = 0; j < ny; j++)
    {
        f5[nx - 1][j] = f5[nx - 2][j];
        f4[nx - 1][j] = f4[nx - 2][j];
        f8[nx - 1][j] = f8[nx - 2][j];

        f6[0][j] = f6[1][j];
        f2[0][j] = f2[1][j];
        f7[0][j] = f7[1][j];
    }

    for (int i = 0; i < nx; i++)
    {

        f1[i][0] = f3[i][0];
        f5[i][0] = f7[i][0];
        f6[i][0] = f8[i][0];

        f3[i][0] = f1[i][0];
        f7[i][0] = f5[i][0];
        f8[i][0] = f6[i][0];
    }

    // stream step
    for (int i = 0; i < nx; i++)
    {
        for (int j = 0; j < ny; j++)
        {
            f0_2[i][j] = f0[i][j];                                 // self
            f1_2[i][(j + ny - 1) % ny] = f1[i][j];                 // up
            f2_2[(i + 1) % nx][j] = f2[i][j];                      // ->
            f3_2[i][(j + 1) % ny] = f3[i][j];                      // down
            f4_2[(i + nx - 1) % nx][j] = f4[i][j];                 // <-
            f5_2[(i + nx - 1) % nx][(j + ny - 1) % ny] = f5[i][j]; // top left
            f6_2[(i + 1) % nx][(j + ny - 1) % ny] = f6[i][j];      // top right
            f7_2[(i + 1) % nx][(j + 1) % ny] = f7[i][j];           // bot right
            f8_2[(i + nx - 1) % nx][(j + 1) % ny] = f8[i][j];      // bot left
        }
    }

    // for (int i = 0; i < nx; i++)
    // {

    //     f1_2[i][0] = f3_2[i][0];
    //     f5_2[i][0] = f7_2[i][0];
    //     f6_2[i][0] = f8_2[i][0];
    // }

    // calculating ux and uy
    for (int i = 0; i < nx; i++)
    {
        for (int j = 0; j < ny; j++)
        {

            rho[i][j] = 0;
            rho[i][j] += f0_2[i][j];
            rho[i][j] += f1_2[i][j];
            rho[i][j] += f2_2[i][j];
            rho[i][j] += f3_2[i][j];
            rho[i][j] += f4_2[i][j];
            rho[i][j] += f5_2[i][j];
            rho[i][j] += f6_2[i][j];
            rho[i][j] += f7_2[i][j];
            rho[i][j] += f8_2[i][j];

            ux[i][j] = 0;
            ux[i][j] += e_x[0] * f0_2[i][j];
            ux[i][j] += e_x[1] * f1_2[i][j];
            ux[i][j] += e_x[2] * f2_2[i][j];
            ux[i][j] += e_x[3] * f3_2[i][j];
            ux[i][j] += e_x[4] * f4_2[i][j];
            ux[i][j] += e_x[5] * f5_2[i][j];
            ux[i][j] += e_x[6] * f6_2[i][j];
            ux[i][j] += e_x[7] * f7_2[i][j];
            ux[i][j] += e_x[8] * f8_2[i][j];
            ux[i][j] /= rho[i][j];

            uy[i][j] = 0;
            uy[i][j] += e_y[0] * f0_2[i][j];
            uy[i][j] += e_y[1] * f1_2[i][j];
            uy[i][j] += e_y[2] * f2_2[i][j];
            uy[i][j] += e_y[3] * f3_2[i][j];
            uy[i][j] += e_y[4] * f4_2[i][j];
            uy[i][j] += e_y[5] * f5_2[i][j];
            uy[i][j] += e_y[6] * f6_2[i][j];
            uy[i][j] += e_y[7] * f7_2[i][j];
            uy[i][j] += e_y[8] * f8_2[i][j];
            uy[i][j] /= rho[i][j];
            // printf("%f\n", ux[i][j]);
        }
    }

    // // collision
    // int i, j;
    // float temp1, temp2, temp3, temp4;
    // for (int k = 0; k < 517; k++)
    // {
    //     i = objx[k];
    //     j = objy[k];
    //     temp1 = f_2[i][j][1];
    //     temp2 = f_2[i][j][5];
    //     temp3 = f_2[i][j][8];
    //     temp4 = f_2[i][j][4];

    //     f_2[i][j][1] = f_2[i][j][3];
    //     f_2[i][j][5] = f_2[i][j][7];
    //     f_2[i][j][8] = f_2[i][j][6];
    //     f_2[i][j][4] = f_2[i][j][2];

    //     f_2[i][j][3] = temp1;
    //     f_2[i][j][7] = temp2;
    //     f_2[i][j][6] = temp3;
    //     f_2[i][j][2] = temp4;

    //     ux[i][j] = 0.0;
    //     uy[i][j] = 0.0;
    // }

    // calculating f_eq (the relaxation)
    float temp_f;
    int k;
    for (int i = 0; i < nx; i++)
    {
        for (int j = 0; j < ny; j++)
        {
            for (int k = 0; k < 9; k++)
            {
                // calculating f_eq (switching between f ad f_2 for convinience)
                // and calculating f for next step
                f0[i][j] = f0_2[i][j] - ((1 / tau) * (f0_2[i][j] - f_eq(i, j, 0)));
                f1[i][j] = f1_2[i][j] - ((1 / tau) * (f0_2[i][j] - f_eq(i, j, 1)));
                f2[i][j] = f2_2[i][j] - ((1 / tau) * (f0_2[i][j] - f_eq(i, j, 2)));
                f3[i][j] = f3_2[i][j] - ((1 / tau) * (f0_2[i][j] - f_eq(i, j, 3))) + 1000 * vof[i][j];
                f4[i][j] = f4_2[i][j] - ((1 / tau) * (f0_2[i][j] - f_eq(i, j, 4)));
                f5[i][j] = f5_2[i][j] - ((1 / tau) * (f0_2[i][j] - f_eq(i, j, 5)));
                f6[i][j] = f6_2[i][j] - ((1 / tau) * (f0_2[i][j] - f_eq(i, j, 6)));
                f7[i][j] = f7_2[i][j] - ((1 / tau) * (f0_2[i][j] - f_eq(i, j, 7)));
                f8[i][j] = f8_2[i][j] - ((1 / tau) * (f0_2[i][j] - f_eq(i, j, 8)));
            }
        }
    }
}

void solveV(void)
{
    stream_and_BC();
}