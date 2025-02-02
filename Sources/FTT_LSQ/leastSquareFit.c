/*
This code is extracted form Chapter 15 of

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

#include <math.h>
#include "nrutil.h"

#define SWAP(a, b)  \
    {               \
        swap = (a); \
        (a) = (b);  \
        (b) = swap; \
    }
#define TOL 1.0e-5

/*
Given a set of data points x[1..ndat], y[1..ndat] with individual standard deviations
sig[1..ndat], use χ2 minimization to fit for some or all of the coefficients a[1..ma] of
a function that depends linearly on a, y = ∑
i ai × afunci(x). The input array ia[1..ma]
indicates by nonzero entries those components of a that should be fitted for, and by zero entries
those components that should be held fixed at their input values. The program returns values
for a[1..ma], χ2 = chisq, and the covariance matrix covar[1..ma][1..ma]. (Parameters
held fixed will return zero covariances.) The user supplies a routine funcs(x,afunc,ma)
that
    returns the ma basis functions evaluated at x = x in the array afunc[1..ma].*/
void lfit(float x[], float y[], float sig[], int ndat, float a[], int ia[],
          int ma, float **covar, float *chisq, void (*funcs)(float, float[], int))
{
    void covsrt(float **covar, int ma, int ia[], int mfit);
    void gaussj(float **a, int n, float **b, int m);
    int i, j, k, l, m, mfit = 0;
    float ym, wt, sum, sig2i, **beta, *afunc;
    beta = matrix(1, ma, 1, 1);
    afunc = vector(1, ma);
    for (j = 1; j <= ma; j++)
        if (ia[j])
            mfit++;
    if (mfit == 0)
        nrerror("lfit: no parameters to be fitted");
    for (j = 1; j <= mfit; j++)
    {
        // Initialize the(symmetric) matrix.
        for (k = 1; k <= mfit; k++)
            covar[j][k] = 0.0;
        beta[j][1] = 0.0;
    }
    for (i = 1; i <= ndat; i++)
    {
        (*funcs)(x[i], afunc, ma);
        ym = y[i];
        if (mfit < ma)
        {
            // Subtract off dependences on known pieces of the fitting function.
            for (j = 1; j <= ma; j++)
                if (!ia[j])
                    ym -= a[j] * afunc[j];
        }
        sig2i = 1.0 / SQR(sig[i]);
        for (j = 0, l = 1; l <= ma; l++)
        {
            if (ia[l])
            {
                wt = afunc[l] * sig2i;
                for (j++, k = 0, m = 1; m <= l; m++)
                    if (ia[m])
                        covar[j][++k] += wt * afunc[m];
                beta[j][1] += ym * wt;
            }
        }
    }
    for (j = 2; j <= mfit; j++)
        // Fill in above the diagonal from symmetry.
        for (k = 1; k < j; k++)
            covar[k][j] = covar[j][k];
    gaussj(covar, mfit, beta, 1);
    // Matrix solution.
    for (j = 0, l = 1; l <= ma; l++)
        if (ia[l])
            a[l] = beta[++j][1];
    // Partition solution to appropriate coefficients
    // a.
    *chisq = 0.0;
    for (i = 1; i <= ndat; i++)
    {
        // Evaluate χ2 of the fit.
        (*funcs)(x[i], afunc, ma);
        for (sum = 0.0, j = 1; j <= ma; j++)
            sum += a[j] * afunc[j];
        *chisq += SQR((y[i] - sum) / sig[i]);
    }
    covsrt(covar, ma, ia, mfit);
    // Sort covariance matrix to true order of fitting coefficients.
    free_vector(afunc, 1, ma);
    free_matrix(beta, 1, ma, 1, 1);
}

/* Expand in storage the covariance matrix covar,
    so as to take into account parameters that are
    being held fixed.(For the latter, return zero covariances.)*/
void covsrt(float **covar, int ma, int ia[], int mfit)
{
    int i, j, k;
    float swap;
    for (i = mfit + 1; i <= ma; i++)
        for (j = 1; j <= i; j++)
            covar[i][j] = covar[j][i] = 0.0;
    k = mfit;
    for (j = ma; j >= 1; j--)
    {
        if (ia[j])
        {
            for (i = 1; i <= ma; i++)
                SWAP(covar[i][k], covar[i][j])
            for (i = 1; i <= ma; i++)
                SWAP(covar[k][i], covar[j][i])
            k--;
        }
    }
}

/*Given a set of data points x[1..ndata],
y[1..ndata] with individual standard deviations
sig[1..ndata],
use χ2 minimization to determine the coefficients a[1..ma] of the fitting function
y = ∑ i ai × afunci(x).
Here we solve the fitting equations using singular value decomposition of the ndata
by ma matrix, as in §2.6. Arrays u[1..ndata][1..ma], v[1..ma][1..ma], and w[1..ma]
rovide workspace on input; on output they define the singular value decomposition,
and can be used to obtain the covariance matrix. The program returns values for the
ma fit parameters a, and χ2, chisq. The user supplies a routine funcs(x,afunc,ma)
that returns the ma basis functions evaluated at x = x in the array
                                                     afunc[1..ma]*/
void svdfit(float x[], float y[], float sig[], int ndata, float a[], int ma,
            float **u, float **v, float w[], float *chisq,
            void (*funcs)(float, float[], int))
{
    void svbksb(float **u, float w[], float **v, int m, int n, float b[],
                float x[]);
    void svdcmp(float **a, int m, int n, float w[], float **v);
    int j, i;
    float wmax, tmp, thresh, sum, *b, *afunc;
    b = vector(1, ndata);
    afunc = vector(1, ma);
    for (i = 1; i <= ndata; i++)
    {
        // Accumulate coefficients of the fitting matrix.
        (*funcs)(x[i], afunc, ma);
        tmp = 1.0 / sig[i];
        for (j = 1; j <= ma; j++)
            u[i][j] = afunc[j] * tmp;
        b[i] = y[i] * tmp;
    }
    svdcmp(u, ndata, ma, w, v);
    // Singular value decomposition.
    wmax = 0.0;
    //     Edit the singular values, given TOL from the
    // #define statement , between here...
    for (j = 1; j <= ma; j++)
        if (w[j] > wmax)
            wmax = w[j];
    thresh = TOL * wmax;
    for (j = 1; j <= ma; j++)
        if (w[j] < thresh)
            w[j] = 0.0;
    // ... and here.
    svbksb(u, w, v, ndata, ma, b, a);
    *chisq = 0.0;
    // Evaluate chi-square.
    for (i = 1; i <= ndata; i++)
    {
        (*funcs)(x[i], afunc, ma);
        for (sum = 0.0, j = 1; j <= ma; j++)
            sum += a[j] * afunc[j];
        *chisq += (tmp = (y[i] - sum) / sig[i], tmp * tmp);
    }
    free_vector(afunc, 1, ma);
    free_vector(b, 1, ndata);
}

/*To evaluate the covariance matrix cvm[1..ma][1..ma] of the fit for ma parameters obtained
by svdfit, call this routine with matrices v[1..ma][1..ma], w[1..ma] as returned from
svdfit.*/
void svdvar(float **v, int ma, float w[], float **cvm)
{
    int k, j, i;
    float sum, *wti;
    wti = vector(1, ma);
    for (i = 1; i <= ma; i++)
    {
        wti[i] = 0.0;
        if (w[i])
            wti[i] = 1.0 / (w[i] * w[i]);
    }
    for (i = 1; i <= ma; i++)
    {
        // Sum contributions to covariance matrix.
        for (j = 1; j <= i; j++)
        {
            for (sum = 0.0, k = 1; k <= ma; k++)
                sum += v[i][k] * v[j][k] * wti[k];
            cvm[j][i] = cvm[i][j] = sum;
        }
    }
    free_vector(wti, 1, ma);
}