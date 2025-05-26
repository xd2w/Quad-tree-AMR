#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "ftt.h"
#include "nrutil.h"

#include "pfplib.h"

#define MAX(x, y) ((x) > (y) ? (x) : (y))
#define MIN(x, y) ((x) < (y) ? (x) : (y))
#define SWAP(a, b) \
  {                \
    temp = a;      \
    a = b;         \
    b = temp;      \
  }

/* ------------------------------------------------------------------- */
/* v is the flux on the cell faces and cc volume fraction */
/* work1, work2, work3 temporal variables providing working allocations */

void octTreeYSwp(int iCell)
{
  int i, j, invx, invz, iLv;
  Real mx, mz, alpha, s1, s2, mm1, mm2, V1, V2;
  Real cc[3][3];

  getCellNgbVOF(iCell, cc);
  i = 1;
  j = 1;

  /*
      s1 = (v[i][j]/dy)*dt(cfl);
      s2 = (v[i+1][j]/dy)*dt(cfl);
  */
  // Real uy = 0.8; // change to actual u later
  // dfetch("uy", &uy);
  Real uy1, uy2, x, y, dx, dy;

  iLv = octLv[iCell / 4];

  x = xCell[iCell];
  y = yCell[iCell];
  dx = dxCell[iLv];
  dy = dyCell[iLv];

  uy1 = (computeVY(x, y) + computeVY(x + dx, y)) * 0.5;
  uy2 = (computeVY(x, y + dy) + computeVY(x + dx, y + dy)) * 0.5;

  s1 = (uy1 / dyCell[iLv]) * global_dt;
  s2 = (uy2 / dyCell[iLv]) * global_dt;

  int botNb = cellNb[2][iCell];
  int topNb = cellNb[3][iCell];
  if (botNb == 0 || topNb == 0)
  {
    printf("***************************************\n");
    printf("Error in calcWorksX: invalid neighbours\n\n");
    exit(1);
  }

  Real ux11, ux12, ux21, ux22;
  Real s11, s12, s21, s22;

  ux11 = (computeVY(x, y) + computeVY(x + 0.5 * dx, y)) * 0.5;
  ux12 = (computeVY(x + 0.5 * dx, y) + computeVY(x + dx, y)) * 0.5;
  ux21 = (computeVY(x, y + dy) + computeVY(x + 0.5 * dx, y + dy)) * 0.5;
  ux22 = (computeVY(x + 0.5 * dx, y + dy) + computeVY(x + dx, y + dy)) * 0.5;

  s11 = (ux11 / dyCell[iLv]) * global_dt;
  s12 = (ux12 / dyCell[iLv]) * global_dt;
  s21 = (ux21 / dyCell[iLv]) * global_dt;
  s22 = (ux22 / dyCell[iLv]) * global_dt;

  if (octLv[botNb / 4] == octLv[iCell / 4] && cellChOct[botNb] != 0)
  {
    s1 = 0.5 * (s11 + s12);
  }

  if (octLv[topNb / 4] == octLv[iCell / 4] && cellChOct[topNb] != 0)
  {
    s2 = 0.5 * (s21 + s22);
  }

  if (cc[i][j] == 0.0)
  {
    return;
  }
  else if (cc[i][j] == 1.0)
  {
    calcWorksYFull(iCell, s1, s2);
  }
  else
  {

    /* normal to the interface */
    mm1 = cc[i - 1][j - 1] + 2. * cc[i - 1][j] + cc[i - 1][j + 1];
    mm2 = cc[i + 1][j - 1] + 2. * cc[i + 1][j] + cc[i + 1][j + 1];
    mx = mm1 - mm2;
    mm1 = cc[i - 1][j - 1] + 2. * cc[i][j - 1] + cc[i + 1][j - 1];
    mm2 = cc[i - 1][j + 1] + 2. * cc[i][j + 1] + cc[i + 1][j + 1];
    mz = mm1 - mm2;

    invz = 0;
    if (mz < 0.)
    {
      mm1 = -s1;
      s1 = -s2;
      s2 = mm1;
      mz = -mz;
      invz = 1;
    }

    invx = 0;
    if (mx < 0.0)
    {
      invx = 1;
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
    if (cc[i][j] <= V1)
      alpha = sqrt(2.0 * cc[i][j] * mm1);
    else if (cc[i][j] <= V2)
      alpha = cc[i][j] + 0.5 * mm1;
    else
      alpha = mm1 + 1.0 - sqrt(2.0 * (1.0 - cc[i][j]) * mm1);

    // /* the new equation of the interface after advection */
    mz = mz / (1.0 - s1 + s2);
    alpha = alpha + mz * s1;

    // mx and mz swapped here
    calcWorksY(iCell, vof[iCell], alpha, mx, mz, invx, invz, s1, s2);
  }

  return;
}

void calcWorksY(int iCell, Real vofVal, Real alpha, Real mx, Real mz, int invx, int invz, Real s1, Real s2)
{
  int botNb, topNb, temp;
  Real V1, V2, V3, mm1, mm2, vol;
  int ltop, lbot, rtop, rbot;

  // -----------

  Real s11, s12, s21, s22;
  Real ux11, ux12, ux21, ux22;
  Real x, y, dx, dy;
  Real calpha, cmz;

  int iLv = octLv[iCell / 4];

  x = xCell[iCell];
  y = yCell[iCell];
  dx = dxCell[iLv];
  dy = dyCell[iLv];

  ux11 = (computeVY(x, y) + computeVY(x + 0.5 * dx, y)) * 0.5;
  ux12 = (computeVY(x + 0.5 * dx, y) + computeVY(x + dx, y)) * 0.5;
  ux21 = (computeVY(x, y + dy) + computeVY(x + 0.5 * dx, y + dy)) * 0.5;
  ux22 = (computeVY(x + 0.5 * dx, y + dy) + computeVY(x + dx, y + dy)) * 0.5;

  s11 = (ux11 / dyCell[iLv]) * global_dt;
  s12 = (ux12 / dyCell[iLv]) * global_dt;
  s21 = (ux21 / dyCell[iLv]) * global_dt;
  s22 = (ux22 / dyCell[iLv]) * global_dt;

  if (invz)
  {
    mm1 = -s11;
    s11 = -s21;
    s21 = mm1;
    mm1 = -s12;
    s12 = -s22;
    s22 = mm1;
  }

  if (invx)
  {
    mm1 = s11;
    s11 = s12;
    s12 = mm1;
    mm1 = s21;
    s21 = s22;
    s22 = mm1;
  }

  // s11 = s12 = s1;
  // s21 = s22 = s2;

  // 2   3
  // 0   1
  // lt rt
  // lb rb
  // 2   3
  // 0   1

  // x- l  r x+

  // the child position if the destination
  ltop = 0;
  lbot = 2;
  rtop = 1;
  rbot = 3;

  botNb = cellNb[2][iCell];
  topNb = cellNb[3][iCell];
  if (botNb == 0 || topNb == 0)
  {
    printf("***************************************\n");
    printf("Error in calcWorksX: invalid neighbours\n\n");
    exit(1);
  }

  if (invz)
  { // inverting top and bot
    SWAP(ltop, lbot)
    SWAP(rtop, rbot)
    SWAP(topNb, botNb)
  }

  if (invx)
  {
    SWAP(ltop, rtop)
    SWAP(lbot, rbot)
  }

  if (vofVal >= 1)
  {
    return;
  }

  // complicated case
  if (0 < vofVal && vofVal < 1)
  {
    mm1 = 1.0 - DMAX(s1, 0.0) + DMIN(s2, 0.0);
    mm2 = alpha - mz * DMAX(s1, 0.0);
    V2 = VOL2(mz, mx, mm2, mm1);

    mm1 = DMAX(-s1, 0.0);
    mm2 = alpha + mx * mm1;
    V1 = VOL2(mx, mz, mm2, mm1);

    mm1 = DMAX(s2, 0.0);
    mm2 = alpha - mx;
    V3 = VOL2(mx, mz, mm2, mm1);

    temp_vof[iCell] += V2;
    // V2 = vofVal;

    if (octLv[botNb / 4] == octLv[iCell / 4] && cellChOct[botNb] != 0)
    { // if neighbour is smaller

      // calpha = alpha + mz * s1;
      // cmz = mz * (1.0 - s1 + s2);

      // top
      mm1 = DMAX(-s12, 0.0);
      mm2 = 2 * (alpha + mz * mm1 - mx * 0.5);
      vol = VOL2(mz, mx, mm2, 2 * mm1);
      temp_vof[4 * cellChOct[botNb] + rbot] += vol;
      V1 -= 0.25 * vol;
      // bottom
      mm1 = DMAX(-s11, 0.0);
      mm2 = 2 * (alpha + mz * mm1);
      vol = VOL2(mz, mx, mm2, 2 * mm1);
      temp_vof[4 * cellChOct[botNb] + lbot] += vol;
      V1 -= 0.25 * vol;

      // correction to V2
      // temp_vof[iCell] += V1;

      // cs1 = 0.5 * (s11 + s12);
    }
    else
    {
      mm1 = DMAX(-s1, 0.0);
      mm2 = alpha + mz * mm1;
      if (octLv[botNb / 4] < octLv[iCell / 4])
      {
        vol = VOL2(mz, mx, mm2, mm1);
        temp_vof[botNb] += 0.25 * vol;
        // V2 -= vol;
      }
      else
      {
        vol = VOL2(mz, mx, mm2, mm1);
        temp_vof[botNb] += vol;
        // V2 -= vol;
      }
    }

    if (octLv[topNb / 4] == octLv[iCell / 4] && cellChOct[topNb] != 0)
    { // if neighbour is smaller
      mm1 = DMAX(s22, 0.0);
      mm2 = 2 * (alpha - mz - mx * 0.5);
      vol = VOL2(mz, mx, mm2, 2 * mm1);
      temp_vof[4 * cellChOct[topNb] + rtop] += vol;
      V3 -= 0.25 * vol;

      mm1 = DMAX(s21, 0.0);
      mm2 = 2 * (alpha - mz);
      vol = VOL2(mz, mx, mm2, 2 * mm1);
      temp_vof[4 * cellChOct[topNb] + ltop] += vol;
      V3 -= 0.25 * vol;

      // correction to V2
      // temp_vof[iCell] += V3;

      // s2 = 0.5 * (s21 + s22);
    }
    else
    {
      mm1 = DMAX(s2, 0.0);
      mm2 = (alpha - mz);
      if (octLv[topNb / 4] < octLv[iCell / 4])
      {
        vol = VOL2(mz, mx, mm2, mm1);
        temp_vof[topNb] += 0.25 * vol;
        V2 -= vol;
      }
      else
      {
        vol = VOL2(mz, mx, mm2, mm1);
        temp_vof[topNb] += vol;
        V2 -= vol;
      }
    }
    // if (V2 < 0)
    // {
    //   printf("\nV2 <0\n");
    //   printf("V2 = %f\n", V2);
    //   // exit(1);
    //   V2 = 0;
    // }

    // temp_vof[iCell] += V2

    return;
  }
}

void calcWorksYFull(int iCell, Real s1, Real s2)
{
  int botNb, topNb, temp;
  Real V1, V2, V3, mm1, mm2;
  int ltop, lbot, rtop, rbot;

  // -----------

  Real s11, s12, s21, s22, cs1, cs2;
  Real ux11, ux12, ux21, ux22;
  Real x, y, dx, dy;

  int iLv = octLv[iCell / 4];

  x = xCell[iCell];
  y = yCell[iCell];
  dx = dxCell[iLv];
  dy = dyCell[iLv];

  ux11 = (computeVY(x, y) + computeVY(x + 0.5 * dx, y)) * 0.5;
  ux12 = (computeVY(x + 0.5 * dx, y) + computeVY(x + dx, y)) * 0.5;
  ux21 = (computeVY(x, y + dy) + computeVY(x + 0.5 * dx, y + dy)) * 0.5;
  ux22 = (computeVY(x + 0.5 * dx, y + dy) + computeVY(x + dx, y + dy)) * 0.5;

  s11 = (ux11 / dyCell[iLv]) * global_dt;
  s12 = (ux12 / dyCell[iLv]) * global_dt;
  s21 = (ux21 / dyCell[iLv]) * global_dt;
  s22 = (ux22 / dyCell[iLv]) * global_dt;

  // s11 = s12 = s1;
  // s21 = s22 = s2;

  // 2   3
  // 0   1
  // lt rt
  // lb rb
  // 2   3
  // 0   1

  // the child position if the destination
  ltop = 0;
  lbot = 2;
  rtop = 1;
  rbot = 3;

  V1 = DMAX(-s1, 0.0);
  // temp_vof[iCell] += 1.0 - DMAX(s1, 0.0) + DMIN(s2, 0.0);
  // V2 = 1;
  cs1 = s1;
  cs2 = s2;
  V3 = DMAX(s2, 0.0);

  botNb = cellNb[2][iCell];
  topNb = cellNb[3][iCell];
  if (botNb == 0 || topNb == 0)
  {
    printf("***************************************\n");
    printf("Error in calcWorksX: invalid neighbours\n\n");
    exit(1);
  }

  // left
  if ((octLv[botNb / 4] == octLv[iCell / 4]) && (cellChOct[botNb] != 0))
  {
    // 2 instead of 4 because V3 = 0.5 should fill up both
    // TODO calculate this separately with different u
    temp_vof[4 * cellChOct[botNb] + rbot] += 2 * DMAX(-s12, 0.0);
    temp_vof[4 * cellChOct[botNb] + lbot] += 2 * DMAX(-s11, 0.0);
    // temp_vof[iCell] += 0.5 * ((fmax(s1, 0.0)) - (fmax(s11, 0.0)));
    // temp_vof[iCell] += 0.5 * ((fmax(s1, 0.0)) - (fmax(s12, 0.0)));
    // V2 -= 0.5 * DMAX(-s11, 0.0);
    // V2 -= 0.5 * DMAX(-s12, 0.0);
    cs1 = 0.5 * (s11 + s12);
  }
  else
  {
    if (octLv[botNb / 4] < octLv[iCell / 4])
    {
      temp_vof[botNb] += 0.25 * V1;
    }
    else
    {
      temp_vof[botNb] += V1;
    }
    // V2 -= V1;
  }

  // right
  if ((octLv[topNb / 4] == octLv[iCell / 4]) && (cellChOct[topNb] != 0))
  {
    // 2 instead of 4 because V3 = 0.5 should fill up both
    // TODO calculate this separately with different u
    temp_vof[4 * cellChOct[topNb] + rtop] += 2 * DMAX(s22, 0.0);
    temp_vof[4 * cellChOct[topNb] + ltop] += 2 * DMAX(s21, 0.0);
    // temp_vof[iCell] -= 0.5 * (fmin(s2, 0.0) - fmin(s21, 0.0));
    // temp_vof[iCell] -= 0.5 * (fmin(s2, 0.0) - fmin(s22, 0.0));
    // V2 -= 0.5 * DMAX(s21, 0.0);
    // V2 -= 0.5 * DMAX(s22, 0.0);
    cs2 = 0.5 * (s21 + s22);
  }
  else
  {
    if (octLv[topNb / 4] < octLv[iCell / 4])
    {
      temp_vof[topNb] += 0.25 * V3;
    }
    else
    {
      temp_vof[topNb] += V3;
    }
    // V2 -= V3;
  }
  // temp_vof[iCell] += V2;
  temp_vof[iCell] += 1.0 - DMAX(cs1, 0.0) + DMIN(cs2, 0.0);
}
