/*
This code is extracted form Chapter 2 of

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
        temp = (a); \
        (a) = (b);  \
        (b) = temp; \
    }

/*Linear equation solution by Gauss - Jordan elimination,
    equation(2.1.1) above.a[1..n][1..n] is the input matrix.b[1..n][1..m] is input
    containing the m right - hand side vectors.
    On output, a is replaced by its matrix inverse, and b is replaced by the
    corresponding set of solution vectors.*/
void gaussj(float **a, int n, float **b, int m)
{
    int *indxc, *indxr, *ipiv;
    int i, icol, irow, j, k, l, ll;
    float big, dum, pivinv, temp;
    indxc = ivector(1, n);
    // The integer arrays ipiv, indxr, and indxc are used for bookkeeping on the pivoting.
    indxr = ivector(1, n);
    ipiv = ivector(1, n);
    for (j = 1; j <= n; j++)
        ipiv[j] = 0;
    for (i = 1; i <= n; i++)
    {
        // This is the main loop over the columns to be reduced.
        big = 0.0;
        for (j = 1; j <= n; j++)
            // This is the outer loop of the search for a pivot element.
            if (ipiv[j] != 1)
                for (k = 1; k <= n; k++)
                {
                    if (ipiv[k] == 0)
                    {
                        if (fabs(a[j][k]) >= big)
                        {
                            big = fabs(a[j][k]);
                            irow = j;
                            icol = k;
                        }
                    }
                    else if (ipiv[k] > 1)
                        nrerror("gaussj: Singular Matrix-1");
                }
        ++(ipiv[icol]);
        // We now have the pivot element, so we interchange rows, if needed,
        //  to put the pivot element on the diagonal.The columns are not physically
        //  interchanged, only relabeled : indxc[i], the column of the ith pivot element,
        //  is the ith column that is reduced, while indxr[i] is the row in which that
        //  pivot element was originally located.If indxr[i]  = indxc[i] there is an
        //  implied column interchange.With this form of bookkeeping, the solution b’s
        //  will end up in the correct order, and the inverse matrix will be scrambled by columns.

        if (irow != icol)
        {
            for (l = 1; l <= n; l++)
                SWAP(a[irow][l], a[icol][l]);
            for (l = 1; l <= m; l++)
                SWAP(b[irow][l], b[icol][l]);
        }
        indxr[i] = irow;
        // We are now ready to divide the pivot row by the pivot element, located at irow and icol.
        indxc[i] = icol;
        if (a[icol][icol] == 0.0)
            nrerror("gaussj: Singular Matrix-2");
        pivinv = 1.0 / a[icol][icol];
        a[icol][icol] = 1.0;
        for (l = 1; l <= n; l++)
            a[icol][l] *= pivinv;
        for (l = 1; l <= m; l++)
            b[icol][l] *= pivinv;

        for (ll = 1; ll <= n; ll++)
            // Next, we reduce the rows...
            if (ll != icol)
            { // ...except for the pivot one, of course.
                dum = a[ll][icol];
                a[ll][icol] = 0.0;
                for (l = 1; l <= n; l++)
                    a[ll][l] -= a[icol][l] * dum;
                for (l = 1; l <= m; l++)
                    b[ll][l] -= b[icol][l] * dum;
            }
    }
    // This is the end of the main loop over columns of the reduction.It only remains to unscram -
    //     ble the solution in view of the column interchanges.We do this by interchanging pairs of
    //         columns in the reverse order that the permutation was built up.
    for (l = n; l >= 1; l--)
    {
        if (indxr[l] != indxc[l])
            for (k = 1; k <= n; k++)
                SWAP(a[k][indxr[l]], a[k][indxc[l]]);
    }
    // And we are done.
    free_ivector(ipiv, 1, n);
    free_ivector(indxr, 1, n);
    free_ivector(indxc, 1, n);
}

/*Solves A·X = B for a vector X, where A is specified by the arrays u[1..m][1..n], w[1..n],
v[1..n][1..n] as returned by svdcmp. m and n are the dimensions of a, and will be equal for
square matrices. b[1..m] is the input right-hand side. x[1..n] is the output solution vector.
No input quantities are destroyed, so the routine may be called sequentially with different b’s.*/
void svbksb(float **u, float w[], float **v, int m, int n, float b[], float x[])
{
    int jj, j, i;
    float s, *tmp;
    tmp = vector(1, n);
    for (j = 1; j <= n; j++)
    {
        // Calculate U T B.
        s = 0.0;
        if (w[j])
        {
            // Nonzero result only if wj is nonzero.
            for (i = 1; i <= m; i++)
                s += u[i][j] * b[i];
            s /= w[j];
            // This is the divide by wj.
        }
        tmp[j] = s;
    }
    for (j = 1; j <= n; j++)
    {
        // Matrix multiply by V to get answer.
        s = 0.0;
        for (jj = 1; jj <= n; jj++)
            s += v[j][jj] * tmp[jj];
        x[j] = s;
    }
    free_vector(tmp, 1, n);
}

/*Given a matrix a[1..m][1..n],
    this routine computes its singular value decomposition, A = U ·W ·V T.
    The matrix U replaces a on output.The diagonal matrix of singular values W
    is out - put as a vector w[1..n].The matrix V(not the transpose V T) is output
    as v[1..n][1..n]. */
void svdcmp(float **a, int m, int n, float w[], float **v)
{
    float pythag(float a, float b);
    int flag, i, its, j, jj, k, l, nm;
    float anorm, c, f, g, h, s, scale, x, y, z, *rv1;
    rv1 = vector(1, n);
    g = scale = anorm = 0.0;
    // Householder reduction to bidiagonal form.
    for (i = 1; i <= n; i++)
    {
        l = i + 1;
        rv1[i] = scale * g;
        g = s = scale = 0.0;
        if (i <= m)
        {
            for (k = i; k <= m; k++)
                scale += fabs(a[k][i]);
            if (scale)
            {
                for (k = i; k <= m; k++)
                {
                    a[k][i] /= scale;
                    s += a[k][i] * a[k][i];
                }
                f = a[i][i];
                g = -SIGN(sqrt(s), f);
                h = f * g - s;
                a[i][i] = f - g;
                for (j = l; j <= n; j++)
                {
                    for (s = 0.0, k = i; k <= m; k++)
                        s += a[k][i] * a[k][j];
                    f = s / h;
                    for (k = i; k <= m; k++)
                        a[k][j] += f * a[k][i];
                }
                for (k = i; k <= m; k++)
                    a[k][i] *= scale;
            }
        }
        w[i] = scale * g;
        g = s = scale = 0.0;
        if (i <= m && i != n)
        {
            for (k = l; k <= n; k++)
                scale += fabs(a[i][k]);
            if (scale)
            {

                for (k = l; k <= n; k++)
                {
                    a[i][k] /= scale;
                    s += a[i][k] * a[i][k];
                }
                f = a[i][l];
                g = -SIGN(sqrt(s), f);
                h = f * g - s;
                a[i][l] = f - g;
                for (k = l; k <= n; k++)
                    rv1[k] = a[i][k] / h;
                for (j = l; j <= m; j++)
                {
                    for (s = 0.0, k = l; k <= n; k++)
                        s += a[j][k] * a[i][k];
                    for (k = l; k <= n; k++)
                        a[j][k] += s * rv1[k];
                }
                for (k = l; k <= n; k++)
                    a[i][k] *= scale;
            }
        }
        anorm = FMAX(anorm, (fabs(w[i]) + fabs(rv1[i])));
    }
    for (i = n; i >= 1; i--)
    {
        // Accumulation of right hand transformations.
        if (i < n)
        {
            if (g)
            {
                // Double division to avoid possible underflow.
                for (j = l; j <= n; j++)
                    v[j][i] = (a[i][j] / a[i][l]) / g;

                for (j = l; j <= n; j++)
                {
                    for (s = 0.0, k = l; k <= n; k++)
                        s += a[i][k] * v[k][j];

                    for (k = l; k <= n; k++)
                        v[k][j] += s * v[k][i];
                }
            }
            for (j = l; j <= n; j++)
                v[i][j] = v[j][i] = 0.0;
        }
        v[i][i] = 1.0;
        g = rv1[i];
        l = i;
    }
    for (i = IMIN(m, n); i >= 1; i--)
    {
        // Accumulation of left hand transformations.
        l = i + 1;
        g = w[i];
        for (j = l; j <= n; j++)
            a[i][j] = 0.0;
        if (g)
        {
            g = 1.0 / g;
            for (j = l; j <= n; j++)
            {
                for (s = 0.0, k = l; k <= m; k++)
                    s += a[k][i] * a[k][j];
                f = (s / a[i][i]) * g;
                for (k = i; k <= m; k++)
                    a[k][j] += f * a[k][i];
            }
            for (j = i; j <= m; j++)
                a[j][i] *= g;
        }
        else
            for (j = i; j <= m; j++)
                a[j][i] = 0.0;
        ++a[i][i];
    }
    for (k = n; k >= 1; k--)
    {
        // Diagonalization of the bidiagonal form : Loop over singular values,
        //      and over allowed iterations.
        for (its = 1; its <= 30; its++)
        {
            flag = 1;
            for (l = k; l >= 1; l--)
            {
                // Test for splitting.
                nm = l - 1;
                // Note that rv1[1] is always zero.
                if ((float)(fabs(rv1[l]) + anorm) == anorm)
                {
                    flag = 0;
                    break;
                }
                if ((float)(fabs(w[nm]) + anorm) == anorm)
                    break;
            }
            if (flag)
            {
                c = 0.0;
                // Cancellation of rv1[l], if l > 1.
                s = 1.0;
                for (i = l; i <= k; i++)
                {
                    f = s * rv1[i];
                    rv1[i] = c * rv1[i];
                    if ((float)(fabs(f) + anorm) == anorm)
                        break;
                    g = w[i];
                    h = pythag(f, g);
                    w[i] = h;
                    h = 1.0 / h;
                    c = g * h;
                    s = -f * h;
                    for (j = 1; j <= m; j++)
                    {
                        y = a[j][nm];
                        z = a[j][i];
                        a[j][nm] = y * c + z * s;
                        a[j][i] = z * c - y * s;
                    }
                }
            }
            z = w[k];
            if (l == k)
            {
                // Convergence.
                if (z < 0.0)
                {
                    // Singular value is made nonnegative.
                    w[k] = -z;
                    for (j = 1; j <= n; j++)
                        v[j][k] = -v[j][k];
                }
                break;
            }
            if (its == 30)
                nrerror("no convergence in 30 svdcmp iterations");
            x = w[l];
            // Shift from bottom 2 - by - 2 minor.
            nm = k - 1;
            y = w[nm];
            g = rv1[nm];
            h = rv1[k];
            f = ((y - z) * (y + z) + (g - h) * (g + h)) / (2.0 * h * y);
            g = pythag(f, 1.0);
            f = ((x - z) * (x + z) + h * ((y / (f + SIGN(g, f))) - h)) / x;
            c = s = 1.0;
            // Next QR transformation :
            for (j = l; j <= nm; j++)
            {
                i = j + 1;
                g = rv1[i];
                y = w[i];
                h = s * g;
                g = c * g;
                z = pythag(f, h);
                rv1[j] = z;
                c = f / z;
                s = h / z;
                f = x * c + g * s;
                g = g * c - x * s;
                h = y * s;
                y *= c;
                for (jj = 1; jj <= n; jj++)
                {
                    x = v[jj][j];
                    z = v[jj][i];
                    v[jj][j] = x * c + z * s;
                    v[jj][i] = z * c - x * s;
                }
                z = pythag(f, h);
                w[j] = z;
                // Rotation can be arbitrary if z = 0.
                if (z)
                {
                    z = 1.0 / z;
                    c = f * z;
                    s = h * z;
                }
                f = c * g + s * y;
                x = c * y - s * g;

                for (jj = 1; jj <= m; jj++)
                {
                    y = a[jj][j];
                    z = a[jj][i];
                    a[jj][j] = y * c + z * s;
                    a[jj][i] = z * c - y * s;
                }
            }
            rv1[l] = 0.0;
            rv1[k] = f;
            w[k] = x;
        }
    }
    free_vector(rv1, 1, n);
}

/*Computes(a2 + b2) 1 /2 without destructive underflow or overflow.*/
float pythag(float a, float b)
{
    float absa, absb;
    absa = fabs(a);
    absb = fabs(b);
    if (absa > absb)
        return absa * sqrt(1.0 + SQR(absb / absa));
    else
        return (absb == 0.0 ? 0.0 : absb * sqrt(1.0 + SQR(absa / absb)));
}