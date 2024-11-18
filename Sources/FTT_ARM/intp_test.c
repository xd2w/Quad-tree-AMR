#include <stdio.h>
#include <math.h>

#include "interpolation.h"
#include "nrutil.h"

#define N 20
#define PI 3.1415926

int main()
{
    int i;
    double yp1, ypn, *x, *y, *y2;

    x = dvector(1, N);
    y = dvector(1, N);
    y2 = dvector(1, N);
    printf("\nsecond-derivatives for sin(x) from 0 to pi\n");
    /* Generate array for interpolation */
    for (i = 1; i <= 20; i++)
    {
        x[i] = i * PI / N;
        y[i] = sin(x[i]);
    }
    /* calculate 2nd derivative with spline */
    yp1 = cos(x[1]);
    ypn = cos(x[N]);
    spline(x, y, N, yp1, yp1, y2);
    /* test result */
    printf("%23s %16s\n", "spline", "actual");
    printf("%11s %14s %16s\n", "angle", "2nd deriv", "2nd deriv");
    for (i = 1; i <= N; i++)
        printf("%10.2f %16.6f %16.6f\n", x[i], y2[i], -sin(x[i]));
    free_dvector(y2, 1, N);
    free_dvector(y, 1, N);
    free_dvector(x, 1, N);
    return 0;
}