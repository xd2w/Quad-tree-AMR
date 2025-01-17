#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "ftt.h"
#include "nrutil.h"
#include "linalg.h"
#include "leastSquareFit.h"

#define TOL 1.0e-5

void plic(void)
{
  // x-sweep
  flagInterfCells();
  propagateFlag(0);
  propagateFlag(1);
  computeXVOF();
  // plotFlagFTT(); exit(1);
  // computeSurfaceNormal();

  // y-sweep
  flagInterfCells();
  propagateFlag(1);
  propagateFlag(1);
  computeYVOF();
}
void computeSurfaceNormal(int iNt)
{
  int iCell, leftNgb, rightNgb;
  Real fraction, MinVof = 1.0e-16, MaxVof = 1.0 - MinVof;
  int *outlist, *stack, list_size;
  char filename[20];
  sprintf(filename, "./DATA/nvec.%03d", iNt);
  FILE *f = fopen(filename, "w");
  float *nvec = vector(1, 2);
  float kappa;

  for (iCell = 0; iCell < numberOfCells; iCell++)
  {
    if ((0 < vof[iCell] && vof[iCell] < 1) && cellFlag[iCell])
    {
      computeN(iCell, nvec);
      kappa = computeKappa(iCell, nvec, -0.5);
      fprintf(f, "%f %f %f %f %f\n", xCell[iCell], yCell[iCell], nvec[1], nvec[2], kappa);
    }
  }
  free_vector(nvec, 1, 2);
  printf("done + closed\n");
  fprintf(f, "\n");
  fclose(f);
}

void computeXVOF(void)
{
  int iCell, leftNgb, rightNgb;
  Real fraction, MinVof = 1.0e-16, MaxVof = 1.0 - MinVof;

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
      leftNgb = cellNb[0][iCell];
      rightNgb = cellNb[1][iCell];
      fraction = work3[leftNgb] + work2[iCell] + work1[rightNgb];
      if (fraction > MaxVof)
        fraction = 1.0;
      if (fraction < MinVof)
        fraction = 0.0;
      vof[iCell] = fraction;
    }
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
      leftNgb = cellNb[2][iCell];
      rightNgb = cellNb[3][iCell];
      fraction = work3[leftNgb] + work2[iCell] + work1[rightNgb];
      if (fraction > MaxVof)
        fraction = 1.0;
      if (fraction < MinVof)
        fraction = 0.0;
      vof[iCell] = fraction;
    }
  }
}

void planefpoly(int iCell, float p[], int deg)
{
  p[1] = 1.0;
  int iLv = octLv[(int)iCell / 4];
  p[2] = xCell[(int)iCell] + 0.5 * dxCell[iLv];
  p[3] = yCell[(int)iCell] + 0.5 * dyCell[iLv];
  // for (int j = 1; j <= deg; j++)
  // {
  //   p[j + 1] = p[j] * x;
  // }
}

void computeN(int iCell, float *n)
{
  int list_size, i, iLv;
  float chisq, *y, *sig, *a, *w, **cvm, **u, **v;
  int *index;
  int deg = 3; // number of coefficients

  int *intfs = ivector(1, 32);
  int *stack = ivector(1, 33);
  list_size = 0;

  neighbouringInterf(iCell, intfs, &list_size, stack);
  // printf("finshed fine \n");
  intfs[list_size] = iCell;
  list_size += 1;

  index = ivector(1, list_size);
  y = vector(1, list_size);
  sig = vector(1, list_size);
  a = vector(1, deg);
  w = vector(1, deg);
  cvm = matrix(1, deg, 1, deg);
  u = matrix(1, list_size, 1, deg);
  v = matrix(1, deg, 1, deg);

  for (i = 0; i < list_size; i++)
  {
    iCell = intfs[i];
    iLv = octLv[(int)iCell / 4];
    index[i + 1] = iCell;
    y[i + 1] = vof[iCell];
    sig[i + 1] = sqrt(dxCell[iLv] * dyCell[iLv]);
  }
  svdfit_mod(index, y, sig, list_size, a, deg, u, v, w, &chisq, planefpoly);
  // svdvar(v, deg, w, cvm);
  float mag = hypot(a[2], a[3]);
  n[1] = -a[2] / mag;
  n[2] = -a[3] / mag;
  for (int k = 1; k <= deg; k++)
  {
    printf("a[%d] = %f\n", k, a[k]);
  }
  printf("\n");
  free_matrix(v, 1, deg, 1, deg);
  free_matrix(u, 1, list_size, 1, deg);
  free_matrix(cvm, 1, deg, 1, deg);
  free_vector(w, 1, deg);
  free_vector(a, 1, deg);
  free_vector(sig, 1, list_size);
  free_vector(y, 1, list_size);
  free_ivector(index, 1, list_size);
}

void fpoly(float x, float p[], int deg)
{
  p[1] = 1.0;
  for (int j = 1; j <= deg; j++)
  {
    p[j + 1] = p[j] * x;
  }
}

float computeKappa(int iCell, float *n, float f_th)
{
  int list_size, i, iLv;
  float chisq, *x, *y, *sig, *a, *w, **cvm, **u, **v;
  int deg = 3; // number of coefficients
  int ogCell = iCell;
  int ogLv = octLv[iCell / 4];

  int *intfs = ivector(1, 32);
  int *stack = ivector(1, 33);
  list_size = 0;

  neighbouringInterf(iCell, intfs, &list_size, stack);
  // printf("finshed fine \n");
  if (f_th > 0)
  {
    intfs[list_size] = iCell;
    list_size += 1;
  }

  if (list_size <= 2)
  {
    printf("insufficint data but probably need refinenign");
    return -1.0;
  }

  x = vector(1, list_size);
  y = vector(1, list_size);
  sig = vector(1, list_size);
  a = vector(1, deg);
  w = vector(1, deg);
  cvm = matrix(1, deg, 1, deg);
  u = matrix(1, list_size, 1, deg);
  v = matrix(1, deg, 1, deg);

  float xp, yp, sfrac;

  for (i = 0; i < list_size; i++)
  {
    iCell = intfs[i];
    iLv = octLv[(int)iCell / 4];
    if (f_th < 0)
    {
      xp = xCell[iCell] + 0.5 * dxCell[iLv];
      yp = yCell[iCell] + 0.5 * dyCell[iLv];
    }
    else
    {
      sfrac = (f_th - vof[ogCell]) / (vof[iCell] - vof[ogCell]);
      if (0 > sfrac || sfrac > 1)
      {
        continue;
      }
      xp = (1 - sfrac) * (xCell[ogCell] + dxCell[ogLv]) + sfrac * (xCell[iCell] + dxCell[iLv]);
      yp = (1 - sfrac) * (yCell[ogCell] + dyCell[ogLv]) + sfrac * (yCell[iCell] + dyCell[iLv]);
    }

    x[i + 1] = xp * n[1] + yp * n[2];
    y[i + 1] = xp * n[2] - yp * n[1];

    sig[i + 1] = sqrt(dxCell[iLv] * dyCell[iLv]);
  }
  svdfit(x, y, sig, list_size, a, deg, u, v, w, &chisq, fpoly);
  // svdvar(v, deg, w, cvm);
  for (int k = 1; k <= deg; k++)
  {
    printf("a[%d] = %f\n", k, a[k]);
  }
  float kappa = fabs(a[3]) / pow(1 + (a[2]) * (a[2]), 1.5);
  printf("\n");
  free_matrix(v, 1, deg, 1, deg);
  free_matrix(u, 1, list_size, 1, deg);
  free_matrix(cvm, 1, deg, 1, deg);
  free_vector(w, 1, deg);
  free_vector(a, 1, deg);
  free_vector(sig, 1, list_size);
  free_vector(y, 1, list_size);
  free_vector(x, 1, list_size);
  return kappa;
}

// modified version of the svdfit from leastSquareFit.h
void svdfit_mod(int index[], float y[], float sig[], int ndata, float a[], int ma,
                float **u, float **v, float w[], float *chisq,
                void (*funcs)(int, float[], int))
{
  void svbksb(float **u, float w[], float **v, int m, int n, float b[], float x[]);
  void svdcmp(float **a, int m, int n, float w[], float **v);
  int j, i;
  float wmax, tmp, thresh, sum, *b, *afunc;
  b = vector(1, ndata);
  afunc = vector(1, ma);
  for (i = 1; i <= ndata; i++)
  {
    // Accumulate coefficients of the fitting matrix.
    (*funcs)(index[i], afunc, ma);
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
    (*funcs)(index[i], afunc, ma);
    for (sum = 0.0, j = 1; j <= ma; j++)
      sum += a[j] * afunc[j];
    *chisq += (tmp = (y[i] - sum) / sig[i], tmp * tmp);
  }
  free_vector(afunc, 1, ma);
  free_vector(b, 1, ndata);
}

void getVOFnonUniformMesh(int iCell, Real1D x, Real1D y, Real1D Lvs, Real1D Vfracs, float f_th)
{
  Int1D stack = ivector(1, 32);
  int stack_size = 0;
  int points = 0;
  int ogCell = iCell;
  float ogVOF = vof[iCell];

  int iOct = iCell / 4;
  int i, iLv;
  double sfrac, l, r0, r1;

  // for (i=0; i<8; i++){
  //   stack[i] =
  // }
  int ngbCell, southCell, northCell;

  // center cell

  ngbCell = cellNb[0][iCell];
  stack[0] = ngbCell;

  southCell = cellNb[2][ngbCell];
  stack[1] = southCell;

  northCell = cellNb[3][ngbCell];
  stack[2] = northCell;

  ngbCell = cellNb[1][iCell];
  stack[3] = ngbCell;

  southCell = cellNb[2][ngbCell];
  stack[4] = southCell;

  northCell = cellNb[3][ngbCell];
  stack[5] = northCell;

  ngbCell = cellNb[2][iCell];
  stack[5] = ngbCell;

  ngbCell = cellNb[3][iCell];
  stack[7] = ngbCell;
  stack_size = 8;

  while (stack_size > 0)
  {
    iCell = stack[stack_size];
    iLv = octLv[(int)iCell / 4];

    if (cellChOct[iCell] == 0)
    {
      // approx for now adds data of cells that will contain the interface
      if (f_th < 0)
      {
        // centre of cell method
        if (0 < vof[iCell] && vof[iCell] < 1)
        {
          Vfracs[points] = vof[iCell];
          x[points] = xCell[iCell] + 0.5 * dxCell[iLv];
          y[points] = yCell[iCell] + 0.5 * dyCell[iLv];
          Lvs[points] = octLv[(int)iCell / 4];
          points += 1;
          stack_size -= 1;
        }
      }
      else
      {
        // thresholding method
        l = sqrt(dxCell[iLv] * dxCell[iLv] + dyCell[iLv] * dyCell[iLv]);
        sfrac = (f_th - ogVOF) / (vof[iCell] - ogVOF);
        if (0 < sfrac && sfrac < 1)
        { // TODO fix this
          x[points] = (sfrac / l) * (xCell[iCell] - xCell[ogCell]);
          y[points] = (sfrac / l) * (yCell[iCell] - yCell[ogCell]);
          Lvs[points] = l;
          points += 1;
        }
        stack_size -= 1;
      }
    }
    else
    {
      iOct = cellChOct[iCell];
      stack[stack_size] = 4 * iOct;
      stack[stack_size + 1] = 4 * iOct + 1;
      stack[stack_size + 2] = 4 * iOct + 2;
      stack[stack_size + 3] = 4 * iOct + 3;
      stack_size += 4;
      if (stack_size > 32)
      {
        printf("cell balance not enforced");
        exit(1);
      }
    }
  }
}

void neighbouringInterf(int iCell, int *outlist, int *list_size, int *stack)
{
  // printf("func reached\n\n");
  int i;
  double sfrac, l, r0, r1;
  int ngbCell, southCell, northCell;

  // Int1D stack = ivector(1, 32);
  int stack_size = 0;
  int ogCell = iCell;
  float ogVOF = vof[iCell];
  int iOct = iCell / 4;
  *list_size = 0;

  ngbCell = cellNb[0][iCell];
  stack[0] = ngbCell;

  southCell = cellNb[2][ngbCell];
  stack[1] = southCell;

  northCell = cellNb[3][ngbCell];
  stack[2] = northCell;

  ngbCell = cellNb[1][iCell];
  stack[3] = ngbCell;

  southCell = cellNb[2][ngbCell];
  stack[4] = southCell;

  northCell = cellNb[3][ngbCell];
  stack[5] = northCell;

  ngbCell = cellNb[2][iCell];
  stack[6] = ngbCell;

  ngbCell = cellNb[3][iCell];
  stack[7] = ngbCell;
  stack_size = 8;

  // printf("stack : %d\n", stack_size);

  while (stack_size > 0)
  {
    stack_size -= 1;
    iCell = stack[stack_size];
    // printf("stack : %d\n", stack_size);
    // printf("iCell : %d\n", iCell);
    // printf("cell child : %d\n", cellChOct[iCell]);
    if (cellChOct[iCell] == 0)
    {
      // approx for now adds data of cells that will contain the interface
      // centre of cell method
      // if (0 < vof[iCell] && vof[iCell] < 1)
      // {
      outlist[*list_size] = iCell;
      ++*list_size;
      // }
    }
    else
    {
      iOct = cellChOct[iCell];
      stack[stack_size] = iOct * 4;
      stack[stack_size + 1] = iOct * 4 + 1;
      stack[stack_size + 2] = iOct * 4 + 2;
      stack[stack_size + 3] = iOct * 4 + 3;
      stack_size += 4;
      if (stack_size > 32)
      {
        printf("cell balance not enforced\n");
        exit(1);
      }
    }
  }
  // free_ivector(stack, 1, 32);
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
    if (iLv == maxLevel)
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
        cellFlag[ngbCell] = 1;
        ngbCell = cellNb[1][iCell];
        cellFlag[ngbCell] = 1;
      }
      else
      {
        ngbCell = cellNb[2][iCell];
        cellFlag[ngbCell] = 1;
        ngbCell = cellNb[3][iCell];
        cellFlag[ngbCell] = 1;
      }
    }
  }
}
