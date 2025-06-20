/*
This code is extracted form Chapter 3 of

Numerical Recipes in C - The Art of Scientific Computing
(Second Edition) published in 1997

Written by:

William H. Press
    Harvard-Smithsonian Center for Astrophysics
Saul A. Teukolsky
    Department of Physics, Cornell University
William T. Vetterling
    Polaroid Corporation
Brian P. Flannery
    EXXON Research and Engineering Company
*/
#include <stdio.h>
#include "nrutil.h"
#include "interpolation.h"
/*
Given arrays x[1..n] and y[1..n] containing a tabulated function, i.e.,
yi = f (xi), with x1 < x2 < . . . < xN , and given values yp1 and ypn for the
first derivative of the interpolating function at points 1 and n, respectively,
this routine returns an array y2[1..n] that contains the second derivatives of
the interpolating function at the tabulated points xi. If yp1 and/or ypn are
equal to 1 × 1030 or larger, the routine is signaled to set the corresponding
boundary condition for a natural spline, with zero second derivative on that
boundary
 */
void spline(double x[], double y[], int n, double yp1, double ypn, double y2[])
{
    printf("0");
    int i, k;
    float p, qn, sig, un, *u;

    u = vector(1, n - 1);
    if (yp1 > 0.99e30)
        // The lower boundary condition is set either to be “natural”
        y2[1] = u[1] = 0.0;
    else
    {
        // or else to have a specified first derivative.
        y2[1] = -0.5;
        u[1] = (3.0 / (x[2] - x[1])) * ((y[2] - y[1]) / (x[2] - x[1]) - yp1);
    }
    // printf("1");

    // This is the decomposition loop of the tridiagonal algorithm.
    // y2 and u are used for temporary storage of the decomposed factors.
    for (i = 2; i <= n - 1; i++)
    {
        sig = (x[i] - x[i - 1]) / (x[i + 1] - x[i - 1]);
        p = sig * y2[i - 1] + 2.0;
        y2[i] = (sig - 1.0) / p;
        u[i] = (y[i + 1] - y[i]) / (x[i + 1] - x[i]) - (y[i] - y[i - 1]) / (x[i] - x[i - 1]);
        u[i] = (6.0 * u[i] / (x[i + 1] - x[i - 1]) - sig * u[i - 1]) / p;
    }
    // printf("2");
    if (ypn > 0.99e30)
        // The upper boundary condition is set either to be “natural”
        qn = un = 0.0;
    else
    {
        // or else to have a specified first derivative.
        qn = 0.5;
        un = (3.0 / (x[n] - x[n - 1])) * (ypn - (y[n] - y[n - 1]) / (x[n] - x[n - 1]));
    }
    y2[n] = (un - qn * u[n - 1]) / (qn * y2[n - 1] + 1.0);

    // This is the backsubstitution loop of the tridiagonalalgorithm.
    for (k = n - 1; k >= 1; k--)
    {
        y2[k] = y2[k] * y2[k + 1] + u[k];
    }
    free_vector(u, 1, n - 1);
}

/*
Given the arrays xa[1..n] and ya[1..n], which tabulate a function (with the xai’s in order),
and given the array y2a[1..n], which is the output from spline above, and given a value of
x, this routine returns a cubic-spline interpolated value y
*/
void splint(double xa[], double ya[], double y2a[], int n, float x, double *y)
{
    void nrerror(char error_text[]);
    int klo, khi, k;
    float h, b, a;

    klo = 1;
    khi = n;
    while (khi - klo > 1)
    {
        k = (khi + klo) >> 1;
        if (xa[k] > x)
        {
            khi = k;
        }
        else
        {
            klo = k;
        }
    }
    h = xa[khi] - xa[klo];
    if (h == 0.0)
    {
        nrerror("Bad xa input to routine splint");
    }
    a = (xa[khi] - x) / h;
    b = (x - xa[klo]) / h;
    *y = a * ya[klo] + b * ya[khi] + ((a * a * a - a) * y2a[klo] + (b * b * b - b) * y2a[khi]) * (h * h) / 6.0;
}