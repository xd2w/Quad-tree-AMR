#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "ftt.h"
#include "nrutil.h"

void plic(void)
{

  // x-sweep
  flagInterfLeaves();
  propagateFlag(0);
  // plotFlagFTT(100);
  propagateFlag(0);
  // plotFlagFTT(101);
  // exit(1);
  computeXVOF();
  restrField(vof);

  // y-sweep
  flagInterfLeaves();
  propagateFlag(1);
  // plotFlagFTT(200);
  propagateFlag(1);
  // plotFlagFTT(201);
  // exit(0);
  computeYVOF();
  restrField(vof);
}
void computeXVOF(void)
{
  int iCell, leftNgb, rightNgb, leftLv, rightLv, myLv;
  Real fraction, MinVof = 1.0e-16, MaxVof = 1.0 - MinVof;
  Real flux1, flux2, flux3;

  /* compute vof1, vof2, vof3 on flaged cells */
  for (iCell = 0; iCell < numberOfCells; iCell++)
  {
    if (cellFlag[iCell])
    {
      octTreeXSwp(iCell);
    }
  }
  /* add vof1, vof2, vof3 on marked cells which are one layer less than
     flaged cells */
  for (iCell = 0; iCell < numberOfCells; iCell++)
  {
    if (cellMark[iCell])
    {
      fraction = temp_vof[iCell];
      if (fraction > MaxVof)
        fraction = 1.0;
      if (fraction < MinVof)
        fraction = 0.0;
      vof[iCell] = fraction;
    }
    temp_vof[iCell] = 0;
  }
}
void computeYVOF(void)
{
  int iCell, leftNgb, rightNgb;
  Real fraction, MinVof = 1.0e-16, MaxVof = 1.0 - MinVof;

  /* compute vof1, vof2, vof3 on flaged cells */
  for (iCell = 0; iCell < numberOfCells; iCell++)
  {
    if (cellFlag[iCell])
    {
      octTreeYSwp(iCell);
    }
  }
  /* add vof1, vof2, vof3 on marked cells which are one layer less than
     flaged cells */
  for (iCell = 0; iCell < numberOfCells; iCell++)
  {
    if (cellMark[iCell])
    {
      fraction = temp_vof[iCell];
      if (fraction > MaxVof)
        fraction = 1.0;
      if (fraction < MinVof)
        fraction = 0.0;
      vof[iCell] = fraction;
    }
    temp_vof[iCell] = 0;
  }
}
void copyCellInt1D(Int1D from, Int1D to)
{
  int iCell;
  for (iCell = 0; iCell < numberOfCells; iCell++)
  {
    to[iCell] = from[iCell];
  }
}

void flagInterfCells(void)
{
  int iCell, iOct, iLv;
  Real fraction;
  for (iCell = 0; iCell < numberOfCells; iCell++)
  {
    cellFlag[iCell] = 0;
    iOct = iCell / cellNumberInOct;
    iLv = octLv[iOct];
    // if(iLv == maxLevel)
    {
      fraction = vof[iCell];
      if (fraction > 0.0 && fraction < 1.0)
      {
        cellFlag[iCell] = 1;
      }
    }
  }
}

void flagInterfLeaves(void)
{
  int iCell, iOct, iLv;
  Real fraction;
  for (iCell = 0; iCell < numberOfCells; iCell++)
  {
    cellFlag[iCell] = 0;
    if (cellChOct[iCell] == 0)
    {
      fraction = vof[iCell];
      if (fraction > 0.0 && fraction < 1.0)
      {
        cellFlag[iCell] = 1;
      }
    }
  }
}

/* propagate flag in direction dir */
void propagateFlag(int dir)
{
  int iCell, ngbCell;

  copyCellInt1D(cellFlag, cellMark);
  for (iCell = 0; iCell < numberOfCells; iCell++)
  {
    if (cellMark[iCell])
    {
      if (dir == 0)
      {
        ngbCell = cellNb[0][iCell];
        if (cellChOct[ngbCell] == 0)
        {
          cellFlag[ngbCell] = 1;
        }
        else
        {
          cellFlag[4 * cellChOct[ngbCell] + 1] = 1;
          cellFlag[4 * cellChOct[ngbCell] + 3] = 1;
        }

        ngbCell = cellNb[1][iCell];
        if (cellChOct[ngbCell] == 0)
        {
          cellFlag[ngbCell] = 1;
        }
        else
        {
          cellFlag[4 * cellChOct[ngbCell] + 0] = 1;
          cellFlag[4 * cellChOct[ngbCell] + 2] = 1;
        }
      }
      else
      {
        ngbCell = cellNb[2][iCell];
        if (cellChOct[ngbCell] == 0)
        {
          cellFlag[ngbCell] = 1;
        }
        else
        {
          cellFlag[4 * cellChOct[ngbCell] + 2] = 1;
          cellFlag[4 * cellChOct[ngbCell] + 3] = 1;
        }

        ngbCell = cellNb[3][iCell];
        if (cellChOct[ngbCell] == 0)
        {
          cellFlag[ngbCell] = 1;
        }
        else
        {
          cellFlag[4 * cellChOct[ngbCell] + 0] = 1;
          cellFlag[4 * cellChOct[ngbCell] + 1] = 1;
        }
      }
    }
  }
}

Real curvature_5x5(Real cc[][6], int ip, int jp, Real dx, Real dy)
{
  Real mm1, mm2, mx1, mz1, mx2, mz2, dndx, dndy;
  int i, j;

  i = 3 + ip;
  j = 2 + jp;
  mm1 = cc[i - 1][j - 1] + 2. * cc[i - 1][j] + cc[i - 1][j + 1];
  mm2 = cc[i + 1][j - 1] + 2. * cc[i + 1][j] + cc[i + 1][j + 1];
  mx1 = mm1 - mm2;

  mm1 = cc[i - 1][j - 1] + 2. * cc[i][j - 1] + cc[i + 1][j - 1];
  mm2 = cc[i - 1][j + 1] + 2. * cc[i][j + 1] + cc[i + 1][j + 1];
  mz1 = mm1 - mm2;

  mm1 = sqrt(mx1 * mx1 + mz1 * mz1) + 1e-50;
  mx1 = mx1 / mm1;
  // mx1 = mx1 / (3 * dx);

  i = 1 + ip;
  j = 2 + jp;
  mm1 = cc[i - 1][j - 1] + 2. * cc[i - 1][j] + cc[i - 1][j + 1];
  mm2 = cc[i + 1][j - 1] + 2. * cc[i + 1][j] + cc[i + 1][j + 1];
  mx2 = mm1 - mm2;

  mm1 = cc[i - 1][j - 1] + 2. * cc[i][j - 1] + cc[i + 1][j - 1];
  mm2 = cc[i - 1][j + 1] + 2. * cc[i][j + 1] + cc[i + 1][j + 1];
  mz2 = mm1 - mm2;

  mm1 = sqrt(mx2 * mx2 + mz2 * mz2) + 1e-50;
  mx2 = mx2 / mm1;
  // mx2 = mx2 / (3 * dx);

  dndx = (mx1 - mx2); // dx;
  // dndx = (mx1 - mx2) / dx; // dx;

  i = 2 + ip;
  j = 3 + jp;
  mm1 = cc[i - 1][j - 1] + 2. * cc[i - 1][j] + cc[i - 1][j + 1];
  mm2 = cc[i + 1][j - 1] + 2. * cc[i + 1][j] + cc[i + 1][j + 1];
  mx1 = mm1 - mm2;

  mm1 = cc[i - 1][j - 1] + 2. * cc[i][j - 1] + cc[i + 1][j - 1];
  mm2 = cc[i - 1][j + 1] + 2. * cc[i][j + 1] + cc[i + 1][j + 1];
  mz1 = mm1 - mm2;

  mm1 = sqrt(mx1 * mx1 + mz1 * mz1) + 1e-50;
  mz1 = mz1 / mm1;
  // mz1 = mz1 / (3 * dy);

  i = 2 + ip;
  j = 1 + jp;
  mm1 = cc[i - 1][j - 1] + 2. * cc[i - 1][j] + cc[i - 1][j + 1];
  mm2 = cc[i + 1][j - 1] + 2. * cc[i + 1][j] + cc[i + 1][j + 1];
  mx2 = mm1 - mm2;

  mm1 = cc[i - 1][j - 1] + 2. * cc[i][j - 1] + cc[i + 1][j - 1];
  mm2 = cc[i - 1][j + 1] + 2. * cc[i][j + 1] + cc[i + 1][j + 1];
  mz2 = mm1 - mm2;

  mm1 = sqrt(mx2 * mx2 + mz2 * mz2) + 1e-50;
  mz2 = mz2 / mm1;
  // mz2 = mz2 / (3 * dy);

  dndy = (mz1 - mz2); // dy;
  // dndy = (mz1 - mz2) / dy; // dy;

  // (mx1 - mx2) / (3 * dx * dx) + (mz1 - mz2) / (3 * dy * dy)
  // printf("\t%f\t %f\t %f\t %f\t %f\t %f\n", mx1, mx2, mz1, mz2, dx, dy);

  // printf("\t%f\t %f\t %f\t %f\n", dndx, dndy, dx, dy);
  // printf("%f \n", dndx + dndy);
  return dndx + dndy;
}

void curvature_6x6(Real cc[][6], double kappas[4], Real dx, Real dy)
{
  Real mm1, mm2, mx1, mz1, mx2, mz2, dndx, dndy;
  int i, j;

  int k = 0;

  for (int jp = 0; jp < 2; jp++)
  {
    for (int ip = 0; ip < 2; ip++)
    {
      kappas[k] = curvature_5x5(cc, ip, jp, dx, dy);
      printf("%f \n", kappas[k]);
      k += 1;
    }
  }
}

void get_m_alpha(int iCell, Real alpha, Real mx, Real mz)
{
  int i, j;
  Real mxp, mzp, mm1, mm2, V1, V2;
  Real cc[3][3];

  getCellNgbVOF(iCell, cc);
  i = 1;
  j = 1;

  /* normal to the interface */
  mm1 = cc[i - 1][j - 1] + 2. * cc[i - 1][j] + cc[i - 1][j + 1];
  mm2 = cc[i + 1][j - 1] + 2. * cc[i + 1][j] + cc[i + 1][j + 1];
  mx = mm1 - mm2;
  mm1 = cc[i - 1][j - 1] + 2. * cc[i][j - 1] + cc[i + 1][j - 1];
  mm2 = cc[i - 1][j + 1] + 2. * cc[i][j + 1] + cc[i + 1][j + 1];
  mz = mm1 - mm2;

  mxp = fabs(mx) + 1.e-50;
  mzp = fabs(mz) + 1.e-50;
  mm2 = DMAX(mxp, mzp);
  mxp = mxp / mm2;
  mzp = mzp / mm2;

  /* get alpha to determine the equation of the interface */
  mm1 = DMIN(mxp, mzp);
  /* compute the two critical volume fraction */
  V1 = 0.5 * mm1;
  V2 = 1.0 - V1;
  if (cc[i][j] <= V1)
    alpha = sqrt(2.0 * cc[i][j] * mm1);
  else if (cc[i][j] <= V2)
    alpha = cc[i][j] + 0.5 * mm1;
  else
    alpha = mm1 + 1.0 - sqrt(2.0 * (1.0 - cc[i][j]) * mm1);

  return;
}

void initialise_temp_vof_to(Real val)
{
  int iCell;
  for (iCell = 0; iCell < numberOfCells; iCell++)
  {
    temp_vof[iCell] = val;
  }
}