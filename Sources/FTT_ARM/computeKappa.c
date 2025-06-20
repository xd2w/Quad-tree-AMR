#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "pfplib.h"
#include "ftt.h"
#include "nrutil.h"
#include "box.h"
#include "pfplib.h"

Real kappaBarickALELike(int iCell, Real cc[][6])
{
    int i, j, ip, jp, iOct, iLv;
    Real ccs[4][4], nxvert[2][2], nyvert[2][2], nmagvert[2][2]; //  cc[6][6],
    Real nmagcell, nxcell, nycell, dx, dy, nmagdx, nmagdy, nxdx, nydy, mag;
    // getCellNgbVOF_6x6(iCell / 4, cc);
    // int px[] = {1, 2, 1, 2};
    // int py[] = {1, 1, 2, 2};

    int px[] = {0, 1, 0, 1};
    int py[] = {0, 0, 1, 1};

    // smooth2ndOrd(cc, ccs);
    // copycc(cc, ccs);
    // smooth2ndOrd(cc, ccs);

    smoothUnifrom(cc, ccs);
    // copycc(cc, ccs);
    // smoothUnifrom(cc, ccs);

    // getCellNgbVOF_unifrom(iCell, ccs);

    iOct = iCell / cellNumberInOct;
    iLv = octLv[iOct];
    dx = dxCell[iLv];
    dy = dyCell[iLv];

    // for (jp = 0; jp < 4; jp++)
    // {
    //     for (ip = 0; ip < 4; ip++)
    //     {
    //         printf("\t%f, ", ccs[ip][jp]);
    //     }
    //     printf("\n");
    // }
    // printf("\n");

    for (ip = 0; ip < 2; ip++)
    {
        for (jp = 0; jp < 2; jp++)
        {
            i = ip + px[iCell % nbNumberOfOct];
            j = jp + py[iCell % nbNumberOfOct];
            nxvert[ip][jp] = (ccs[i + 1][j] + ccs[i + 1][j + 1] - ccs[i][j] - ccs[i][j + 1]) / (2 * dx + 1e-50);
            nyvert[ip][jp] = (ccs[i][j + 1] + ccs[i + 1][j + 1] - ccs[i][j] - ccs[i + 1][j]) / (2 * dy + 1e-50);
            nmagvert[ip][jp] = sqrt(nxvert[ip][jp] * nxvert[ip][jp] + nyvert[ip][jp] * nyvert[ip][jp]);
            // printf("nx = %f, ny = %f, |n| = %f\n", nxvert[ip][jp], nyvert[ip][jp], nmagvert[ip][jp]);
        }
    }
    i = 0;
    j = 0;
    nxcell = 0.25 * (nxvert[0][0] + nxvert[1][0] + nxvert[0][1] + nxvert[1][1]);
    nycell = 0.25 * (nyvert[0][0] + nyvert[0][1] + nyvert[1][0] + nyvert[1][1]);
    nmagcell = 0.25 * (nmagvert[0][0] + nmagvert[1][0] + nmagvert[0][1] + nmagvert[1][1]);
    // nmagcell = sqrt(nxcell * nxcell + nycell * nycell);
    mag = 1 / (nmagcell + 1e-50);

    nmagdx = (nmagvert[i + 1][j] + nmagvert[i + 1][j + 1] - nmagvert[i][j] - nmagvert[i][j + 1]) / (2 * dx + 1e-50);
    nmagdy = (nmagvert[i][j + 1] + nmagvert[i + 1][j + 1] - nmagvert[i][j] - nmagvert[i + 1][j]) / (2 * dy + 1e-50);

    nxdx = (nxvert[i + 1][j] + nxvert[i + 1][j + 1] - nxvert[i][j] - nxvert[i][j + 1]) / (2 * dx + 1e-50);
    nydy = (nyvert[i][j + 1] + nyvert[i + 1][j + 1] - nyvert[i][j] - nyvert[i + 1][j]) / (2 * dy + 1e-50);

    return fabs(mag * ((nxcell * nmagdx + nycell * nmagdy) * mag - (nxdx + nydy)));
    // return fabs((1 / (nmagcell*nmagcell)) * ((nxcell / nmagcell) * nmagdx + (nycell / nmagcell) * nmagdy) - (nxdx + nydy));
}

Real kappaBarickALELike_wider(int iCell, Real cc[][6])
{
    int i, j, ip, jp, iOct, iLv;
    Real ccs[4][4], nxvert[3][3], nyvert[3][3], nmagvert[3][3]; //  cc[6][6],
    Real nmagcell, nxcell, nycell, dx, dy, nmagdx, nmagdy, nxdx, nydy, mag;

    getCellNgbVOF_6x6(iCell / 4, cc);
    int px[] = {1, 2, 1, 2};
    int py[] = {1, 1, 2, 2};

    // smooth2ndOrd(cc, ccs);
    // copycc(cc, ccs);
    // smooth2ndOrd(cc, ccs);

    // smoothUnifrom(cc, ccs);
    // copycc(cc, ccs);
    // smoothUnifrom(cc, ccs);

    iOct = iCell / cellNumberInOct;
    iLv = octLv[iOct];
    dx = dxCell[iLv];
    dy = dyCell[iLv];

    // for (jp = 0; jp < 4; jp++)
    // {
    //     for (ip = 0; ip < 4; ip++)
    //     {
    //         printf("\t%f, ", ccs[ip][jp]);
    //     }
    //     printf("\n");
    // }
    // printf("\n");
    Real mm1, mm2;

    for (ip = 0; ip < 3; ip++)
    {
        for (jp = 0; jp < 3; jp++)
        {
            i = ip + px[iCell % nbNumberOfOct];
            j = jp + py[iCell % nbNumberOfOct];
            mm1 = cc[i - 1][j - 1] + 2. * cc[i - 1][j] + cc[i - 1][j + 1];
            mm2 = cc[i + 1][j - 1] + 2. * cc[i + 1][j] + cc[i + 1][j + 1];
            nxvert[ip][jp] = 0.5 * (mm1 - mm2) / dx;
            mm1 = cc[i - 1][j - 1] + 2. * cc[i][j - 1] + cc[i + 1][j - 1];
            mm2 = cc[i - 1][j + 1] + 2. * cc[i][j + 1] + cc[i + 1][j + 1];
            nyvert[ip][jp] = 0.5 * (mm1 - mm2) / dy;
            nmagvert[ip][jp] = sqrt(nxvert[ip][jp] * nxvert[ip][jp] + nyvert[ip][jp] * nyvert[ip][jp]);
            // printf("nx = %f, ny = %f, |n| = %f\n", nxvert[ip][jp], nyvert[ip][jp], nmagvert[ip][jp]);
        }
    }
    i = 1;
    j = 1;

    // nxcell = nxvert[i][j];
    // nycell = nyvert[i][j];
    // nmagcell = sqrt(nxvert[i][j] * nxvert[i][j] + nyvert[i][j] * nyvert[i][j])

    nxcell = (nxvert[0][0] + nxvert[1][0] + nxvert[2][0] +
              nxvert[0][1] + nxvert[1][1] + nxvert[2][1] +
              nxvert[0][2] + nxvert[1][2] + nxvert[2][2]) /
             9;

    nycell = (nyvert[0][0] + nyvert[1][0] + nyvert[2][0] +
              nyvert[0][1] + nyvert[1][1] + nyvert[2][1] +
              nyvert[0][2] + nyvert[1][2] + nyvert[2][2]) /
             9;

    nmagcell = (nmagvert[0][0] + nmagvert[1][0] + nmagvert[2][0] +
                nmagvert[0][1] + nmagvert[1][1] + nmagvert[2][1] +
                nmagvert[0][2] + nmagvert[1][2] + nmagvert[2][2]) /
               9;

    // nmagcell = sqrt(nxcell * nxcell + nycell * nycell);
    mag = 1 / (nmagcell + 1e-50);

    mm1 = nmagvert[i - 1][j - 1] + 2. * nmagvert[i - 1][j] + nmagvert[i - 1][j + 1];
    mm2 = nmagvert[i + 1][j - 1] + 2. * nmagvert[i + 1][j] + nmagvert[i + 1][j + 1];
    nmagdx = (mm1 - mm2) / (2 * dx + 1e-50);

    mm1 = nmagvert[i - 1][j - 1] + 2. * nmagvert[i][j - 1] + nmagvert[i + 1][j - 1];
    mm2 = nmagvert[i - 1][j + 1] + 2. * nmagvert[i][j + 1] + nmagvert[i + 1][j + 1];
    nmagdy = (mm1 - mm2) / (2 * dy + 1e-50);

    // del.n
    mm1 = nxvert[i - 1][j - 1] + 2. * nxvert[i - 1][j] + nxvert[i - 1][j + 1];
    mm2 = nxvert[i + 1][j - 1] + 2. * nxvert[i + 1][j] + nxvert[i + 1][j + 1];
    nxdx = (mm1 - mm2) / (2 * dx + 1e-50);

    mm1 = nyvert[i - 1][j - 1] + 2. * nyvert[i][j - 1] + nyvert[i + 1][j - 1];
    mm2 = nyvert[i - 1][j + 1] + 2. * nyvert[i][j + 1] + nyvert[i + 1][j + 1];
    nydy = (mm1 - mm2) / (2 * dy + 1e-50);

    return fabs(mag * ((nxcell * nmagdx + nycell * nmagdy) * mag - (nxdx + nydy)));
    // return fabs((1 / (nmagcell*nmagcell)) * ((nxcell / nmagcell) * nmagdx + (nycell / nmagcell) * nmagdy) - (nxdx + nydy));
}

Real kappaMeier(int iCell)
{
    int i, j, invx, invz, iLv;
    Real mx, mz, alpha, s1, s2, mm1, mm2, V1, V2, dx;
    Real cc[3][3], A, sigma, xi, omega, kappa1, kappa2, kappa3;
    getCellNgbVOF(iCell, cc);

    iLv = octLv[iCell / nbNumberOfOct];
    dx = dxCell[iLv];

    i = 1;
    j = 1;

    mm1 = cc[i - 1][j - 1] + 2. * cc[i - 1][j] + cc[i - 1][j + 1];
    mm2 = cc[i + 1][j - 1] + 2. * cc[i + 1][j] + cc[i + 1][j + 1];
    mx = mm1 - mm2;
    mm1 = cc[i - 1][j - 1] + 2. * cc[i][j - 1] + cc[i + 1][j - 1];
    mm2 = cc[i - 1][j + 1] + 2. * cc[i][j + 1] + cc[i + 1][j + 1];
    mz = mm1 - mm2;

    mx = fabs(mx) + 1.e-50;
    mz = fabs(mz) + 1.e-50;
    mm2 = DMAX(mx, mz);
    mm1 = DMIN(mx, mz);
    mx = mx / mm2;
    mz = mz / mm2;

    omega = atan2(mm1, mm2);

    /* get alpha to determine the equation of the interface */
    mm1 = DMIN(mx, mz);
    /* compute the two critical volume fraction */
    V1 = 0.5 * mm1;
    V2 = 1.0 - V1;
    if (cc[i][j] <= V1)
        alpha = sqrt(2.0 * cc[i][j] * mm1);
    else if (cc[i][j] <= V2)
        alpha = cc[i][j] + 0.5 * mm1;
    else
        alpha = mm1 + 1.0 - sqrt(2.0 * (1.0 - cc[i][j]) * mm1);

    mm2 = (alpha + mx + mz) / 3;
    A = 9 * VOL2(mx, mz, mm2, 1);
    sigma = 0;
    for (i = 0; i < 3; i++)
    {
        for (j = 0; j < 3; j++)
        {
            sigma += cc[i][j];
        }
    }
    sigma = sigma - A;

    kappa1 = (est1[0] + est1[1] * xi + est1[2] * omega + est1[3] * sigma);
    kappa1 += (est1[4] * xi * xi + est1[5] * omega * omega + est1[6] * sigma * sigma);
    kappa1 += (est1[7] * xi * omega + est1[8] * xi * sigma + est1[9] * omega * sigma);
    kappa1 += (est1[10] * xi * xi * xi + est1[11] * omega * omega * omega + est1[12] * sigma * sigma * sigma);
    kappa1 += (est1[13] * xi * xi * omega + est1[14] * xi * xi * sigma);
    kappa1 += (est1[15] * omega * omega * xi + est1[16] * omega * omega * sigma);
    kappa1 += (est1[17] * sigma * sigma * xi + est1[18] * sigma * sigma * omega);
    kappa1 += (est1[19] * xi * omega * sigma);
    kappa1 = fabs(kappa1);

    if (kappa1 > 0.4)
    {
        // printf("check 1\n");
        return kappa1;
    }

    kappa2 = (est2[0] + est2[1] * xi + est2[2] * omega + est2[3] * sigma);
    kappa2 += (est2[4] * xi * xi + est2[5] * omega * omega + est2[6] * sigma * sigma);
    kappa2 += (est2[7] * xi * omega + est2[8] * xi * sigma + est2[9] * omega * sigma);
    kappa2 += (est2[10] * xi * xi * xi + est2[11] * omega * omega * omega + est2[12] * sigma * sigma * sigma);
    kappa2 += (est2[13] * xi * xi * omega + est2[14] * xi * xi * sigma);
    kappa2 += (est2[15] * omega * omega * xi + est2[16] * omega * omega * sigma);
    kappa2 += (est2[17] * sigma * sigma * xi + est2[18] * sigma * sigma * omega);
    kappa2 += (est2[19] * xi * omega * sigma);
    kappa2 = fabs(kappa2);

    if (kappa2 > 0.1)
    {
        // printf("check 2\n");
        return kappa2;
    }

    kappa3 = (est3[0] + est3[1] * xi + est3[2] * omega + est3[3] * sigma);
    kappa3 += (est3[4] * xi * xi + est2[5] * omega * omega + est3[6] * sigma * sigma);
    kappa3 += (est3[7] * xi * omega + est3[8] * xi * sigma + est3[9] * omega * sigma);
    kappa3 += (est3[10] * xi * xi * xi + est3[11] * omega * omega * omega + est3[12] * sigma * sigma * sigma);
    kappa3 += (est3[13] * xi * xi * omega + est3[14] * xi * xi * sigma);
    kappa3 += (est3[15] * omega * omega * xi + est3[16] * omega * omega * sigma);
    kappa3 += (est3[17] * sigma * sigma * xi + est3[18] * sigma * sigma * omega);
    kappa3 += (est3[19] * xi * omega * sigma);
    kappa3 = fabs(kappa3);
    // printf("check 3\n");

    return kappa3;

    // printf("***********************\n");
    // printf("Error kappa est. not working\n\n");
}

Real kappaHF(int iCell, Real cc[][6])
{
    int i, iOct, iLocal, iLv;
    Real mx, mz, hx, hxx, delta, dx, dy, kappa1, kappa2, kappa;
    mx = mxCell[iCell];
    mz = mzCell[iCell];
    iOct = iCell / 4;
    iLocal = iCell % 4;
    iLv = octLv[iOct];

    dx = dxCell[iLv];
    dy = dyCell[iLv];

    Real h[3] = {0, 0, 0};

    if (fabs(mx) < fabs(mz))
    { // add up vertically
        for (i = 0; i < 6; i++)
        {
            // 1 or 3   % 2 = 1
            // 0 or 2   % 2 = 0
            h[0] += cc[1 + (iLocal % 2)][i];
            h[1] += cc[2 + (iLocal % 2)][i];
            h[2] += cc[3 + (iLocal % 2)][i];
        }
        delta = dy / (dx * dx + 1e-50);
        // hx = (h[2] - h[0]) / 2;
        // hxx = (h[2] - 2 * h[1] + h[0]);

        // kappa1 = fabs(hxx / (pow(1 + hx * hx, 1.5) + 1e-50)) / (delta + 1e-50);
    }
    else
    { // add up horizontally
        for (i = 0; i < 6; i++)
        {
            // 2 or 3   / 2 = 1
            // 0 or 1   / 2 = 0
            h[0] += cc[i][1 + (iLocal / 2)];
            h[1] += cc[i][2 + (iLocal / 2)];
            h[2] += cc[i][3 + (iLocal / 2)];
        }
        delta = dx / (dy * dy + 1e-50);
    }
    hx = (h[2] - h[0]) / 2;
    hxx = (h[2] - 2 * h[1] + h[0]);

    kappa = fabs(hxx / (pow(1 + hx * hx, 1.5) + 1e-50)) * delta;

    // kappa = DMIN(kappa1, kappa2);

    // if (h[0] == 0 || h[1] == 0 || h[2] == 0)
    // {
    //     kappa = -1;
    //     printf("0\n");
    //     printf("VOF = %f\n", vof[iCell]);
    //     printf("mx = %f, mz = %f \n", mx, mz);
    //     printf("h0 = %f, h1 = %f, h2 = %f", h[0], h[1], h[2]);
    //     printcc6(cc);
    //     exit(1);
    // }

    // if (h[0] == 6 || h[1] == 6 || h[2] == 6)
    // {
    //     kappa = -1;
    //     printf("6\n");
    //     printf("VOF = %f\n", vof[iCell]);
    //     printf("mx = %f, mz = %f \n", mx, mz);
    //     printf("h0 = %f, h1 = %f, h2 = %f", h[0], h[1], h[2]);
    //     printcc6(cc);
    //     exit(1);
    // }

    // printcc6(cc);

    // printf("\nh[0] = %f, h[1] = %f, h[2] = %f, d = %f\n", h[0], h[1], h[2], delta);
    // printf("mx = %f, mz = %f\n", mx, mz);
    // printf("hx = %f, hxx = %f\n", hx, hxx);
    // printf("kappa = %f\n", kappa);

    return kappa;
}

Real kappaHF3x3(int iCell, Real cc[][6])
{
    int i, iOct, iLocal, iLv;
    Real mx, mz, hx, hxx, delta, dx, dy, kappa, ccs[3][3];
    mx = mxCell[iCell];
    mz = mzCell[iCell];
    iOct = iCell / 4;
    iLocal = iCell % 4;
    iLv = octLv[iOct];

    dx = dxCell[iLv];
    dy = dyCell[iLv];

    getCellNgbVOF(iCell, ccs);

    Real h[3] = {0, 0, 0};

    if (fabs(mx) < fabs(mz))
    { // vertical
        for (i = 0; i < 3; i++)
        {
            // 1 or 3   % 2 = 1
            // 0 or 2   % 2 = 0
            h[0] += ccs[0][i];
            h[1] += ccs[1][i];
            h[2] += ccs[2][i];
        }
        delta = dx;
    }
    else
    { // horizontal
        for (i = 0; i < 3; i++)
        {
            // 2 or 3   / 2 = 1
            // 0 or 1   / 2 = 0
            h[0] += ccs[i][0];
            h[1] += ccs[i][1];
            h[2] += ccs[i][2];
        }
        delta = dy;
    }

    hx = (h[2] - h[0]) / (2 * delta);
    hxx = (h[2] - 2 * h[1] + h[0]) / (delta * delta);

    kappa = fabs(hxx / (pow(1 + hx * hx, 1.5) + 1e-50));

    printcc3(ccs);

    printf("\nh[0] = %f, h[1] = %f, h[2] = %f, d = %f\n", h[0], h[1], h[2], delta);
    printf("mx = %f, mz = %f\n", mx, mz);
    printf("hx = %f, hxx = %f\n", hx, hxx);
    printf("kappa = %f\n", kappa);

    return kappa;
}

void smooth2ndOrd8(Real from[][6], Real to[][4])
{
    int i, j;
    for (i = 1; i < 5; i++)
    {
        for (j = 1; j < 5; j++)
        {
            to[i - 1][j - 1] = (256. * from[i][j] +
                                81. * (from[i + 1][j] + from[i - 1][j] + from[i][j + 1] + from[i][j - 1]) +
                                1. * (from[i + 1][j + 1] + from[i - 1][j + 1] + from[i + 1][j + 1] + from[i + 1][j - 1])) /
                               36.;
            // printf("\t%f,", to[i - 1][j - 1]);
        }
        // printf("\n");
    }
    // printf("\n");
}

void smooth2ndOrd(Real from[][6], Real to[][4])
{
    int i, j;
    for (i = 1; i < 5; i++)
    {
        for (j = 1; j < 5; j++)
        {
            to[i - 1][j - 1] = (16. * from[i][j] +
                                4. * (from[i + 1][j] + from[i - 1][j] + from[i][j + 1] + from[i][j - 1]) +
                                1. * (from[i + 1][j + 1] + from[i - 1][j + 1] + from[i + 1][j + 1] + from[i + 1][j - 1])) /
                               36.;
            // printf("\t%f,", to[i - 1][j - 1]);
        }
        // printf("\n");
    }
    // printf("\n");
}

void smoothUnifrom(Real from[][6], Real to[][4])
{
    int i, j;
    for (i = 1; i < 5; i++)
    {
        for (j = 1; j < 5; j++)
        {
            to[i - 1][j - 1] = (from[i - 1][j - 1] + from[i][j - 1] + from[i + 1][j - 1] +
                                from[i - 1][j] + from[i][j] + from[i + 1][j] +
                                from[i - 1][j + 1] + from[i][j + 1] + from[i + 1][j + 1]) /
                               9;
        }
    }
}

void copycc(Real cc[][6], Real ccs[][4])
{
    int i, j;
    Real ave = 0;
    for (i = 0; i < 4; i++)
    {
        for (j = 0; j < 4; j++)
        {
            // ave += ccs[i][j];
            cc[i + 1][j + 1] = ccs[i][j];
        }
    }
    // printf("ave=%f\n", ave);
    ave /= 16;
    for (i = 0; i < 6; i++)
    {
        // cc[i][0] = ave;
        // cc[i][5] = ave;

        cc[i][0] = cc[i][1];
        cc[i][5] = cc[i][4];
    }
    for (i = 1; i < 5; i++)
    {
        cc[5][i] = cc[4][i];
        cc[0][i] = cc[1][i];
    }

    cc[0][0] = ccs[0][0];
    cc[5][0] = ccs[3][0];
    cc[0][5] = ccs[0][3];
    cc[5][5] = ccs[3][3];
}

void copyCellReal1D(Real1D from, Real1D to)
{
    for (int iCell = 0; iCell < maxNumberOfCells; iCell++)
    {
        to[iCell] = from[iCell];
    }
}

void smooth1D(void)
{
    int iCell, i, j;
    Real cc[3][3];
    i = 1;
    j = 1;
    for (iCell = 0; iCell < maxNumberOfCells; iCell++)
    {
        if (cellChOct[iCell] == 0)
        {
            getCellNgbTempVOF(iCell, cc);
            vof[iCell] = (16. * cc[i][j] +
                          4. * (cc[i + 1][j] + cc[i - 1][j] + cc[i][j + 1] + cc[i][j - 1]) +
                          1. * (cc[i + 1][j + 1] + cc[i - 1][j + 1] + cc[i + 1][j + 1] + cc[i + 1][j - 1])) /
                         36.;
        }
    }
}

void smooth1DUnifrom(void)
{
    int iCell, i, j;
    Real cc[3][3];
    i = 1;
    j = 1;
    for (iCell = 0; iCell < maxNumberOfCells; iCell++)
    {
        if (cellChOct[iCell] == 0)
        {
            getCellNgbTempVOF(iCell, cc);
            temp_vof[iCell] = (cc[i][j] +
                               (cc[i + 1][j] + cc[i - 1][j] + cc[i][j + 1] + cc[i][j - 1]) +
                               (cc[i + 1][j + 1] + cc[i - 1][j + 1] + cc[i + 1][j + 1] + cc[i + 1][j - 1])) /
                              9;
        }
    }
}