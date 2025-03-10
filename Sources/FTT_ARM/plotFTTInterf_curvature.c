#include <stdio.h>
#include <stdlib.h>
#include "ftt.h"
#include "nrutil.h"
#include <stdbool.h>
#include <math.h>

void plotFTTInterf(int ndata)
{
   int iCell, i, i1, i2, i3;
   Real fraction;
   char fname[] = "DATA/intf.000";
   FILE *finterf;

   i = ndata;
   i1 = i % 10;
   i /= 10;
   i2 = i % 10;
   i /= 10;
   i3 = i % 10;
   fname[10] = '0' + i3;
   fname[11] = '0' + i2;
   fname[12] = '0' + i1;

   finterf = fopen(fname, "w");
   //  iCell = 339;
   //  plotCellInterf(iCell, finterf);
   for (iCell = 0; iCell < numberOfCells; iCell++)
   {
      if (cellChOct[iCell] == 0)
      {
         fraction = vof[iCell];
         if (fraction > 0.0 && fraction < 1.0)
         {
            plotCellInterf(iCell, finterf);
         }
      }
   }

   fclose(finterf);
   return;
}

void plotFTTInterfAtLevel(int ndata, int level)
{
   int iCell, i, i1, i2, i3;
   Real fraction;
   char fname[] = "DATA/intf_lev0.000";
   FILE *finterf;

   i = ndata;
   i1 = i % 10;
   i /= 10;
   i2 = i % 10;
   i /= 10;
   i3 = i % 10;
   fname[13] = '0' + level;
   fname[15] = '0' + i3;
   fname[16] = '0' + i2;
   fname[17] = '0' + i1;

   finterf = fopen(fname, "w");
   //  iCell = 339;
   //  plotCellInterf(iCell, finterf);
   for (iCell = 0; iCell < numberOfCells; iCell++)
   {
      if (octLv[iCell / 4] == level)
      {
         fraction = vof[iCell];
         if (fraction > 0.0 && fraction < 1.0)
         {
            plotCellInterf(iCell, finterf);
         }
      }
   }

   fclose(finterf);
   return;
}

/* fill the 3x3 bloc of VOF centered at cell iCell */
void getCellNgbVOF_unifrom(int iCell, Real cc[][3])
{
   int ngbCell, southCell, northCell;
   int iLv, nbLv;

   // center cell
   cc[1][1] = vof[iCell];
   ngbCell = cellNb[0][iCell];
   // west cell
   cc[0][1] = vof[ngbCell];
   // south-west cell
   southCell = cellNb[2][ngbCell];
   cc[0][0] = vof[southCell];
   // north-west cell
   northCell = cellNb[3][ngbCell];
   cc[0][2] = vof[northCell];
   // est cell
   ngbCell = cellNb[1][iCell];
   cc[2][1] = vof[ngbCell];
   // south-est cell
   southCell = cellNb[2][ngbCell];
   cc[2][0] = vof[southCell];
   // north-est cell
   northCell = cellNb[3][ngbCell];
   cc[2][2] = vof[northCell];
   // south cell
   ngbCell = cellNb[2][iCell];
   cc[1][0] = vof[ngbCell];
   // north cell
   ngbCell = cellNb[3][iCell];
   cc[1][2] = vof[ngbCell];
}

void getCellNgbVOF(int iCell, Real cc[][3])
{
   int ngbCell, southCell, northCell;
   int iLv, nbLv;

   // center cell
   cc[1][1] = vof[iCell];
   ngbCell = cellNb[0][iCell];

   // west cell
   cc[0][1] = getVOF(ngbCell, iCell, 0);
   // south-west cell
   southCell = cellNb[2][ngbCell];
   cc[0][0] = getVOF(southCell, ngbCell, 2);
   // north-west cell
   northCell = cellNb[3][ngbCell];
   cc[0][2] = getVOF(northCell, ngbCell, 3);

   // est cell
   ngbCell = cellNb[1][iCell];
   cc[2][1] = getVOF(ngbCell, iCell, 1);
   // south-est cell
   southCell = cellNb[2][ngbCell];
   cc[2][0] = getVOF(southCell, ngbCell, 2);
   // north-est cell
   northCell = cellNb[3][ngbCell];
   cc[2][2] = getVOF(northCell, ngbCell, 3);

   // south cell
   ngbCell = cellNb[2][iCell];
   cc[1][0] = getVOF(ngbCell, iCell, 2);
   // north cell
   ngbCell = cellNb[3][iCell];
   cc[1][2] = getVOF(ngbCell, iCell, 3);
}

/* fill the 3x3 bloc of VOF centered at cell iCell */
void _getCellNgbVOF(int iCell, Real cc[][3])
{
   int i, j, k;
   int iOct;
   int nbs[2], nb_pointer;
   iOct = iCell % 4;

   // offet depending on the values of iCell % 4 to make iCell centre of 3x3
   int px[] = {1, 0, 1, 0};
   int pz[] = {1, 1, 0, 0};

   i = px[iCell % 4];
   j = pz[iCell % 4];
   cc[i][j] = vof[4 * iOct + 0];
   cc[i + 1][j] = vof[4 * iOct + 1];
   cc[i][j + 1] = vof[4 * iOct + 2];
   cc[i + 1][j + 1] = vof[4 * iOct + 3];

   nb_pointer = 0;
   for (int dir = 0; dir < 4; dir++)
   {
      // loking for outer neighbour
      if (morton_lookup[dir][iCell % 4] > 4)
      {
         nbs[nb_pointer++] = dir;
      }
   }

   // horizontal Nb
   nbs[0] = cellNb[nbs[0]][iCell];
}

Real getVOF(int nbCell, int iCell, int dir)
{
   int cOct, k;
   if (octLv[nbCell / cellNumberInOct] > octLv[iCell / cellNumberInOct])
   {
      cOct = cellChOct[nbCell];
      if (cOct != 0)
      {
         k = morton_lookup[dir][iCell % 4];
         if (k < 4)
         {
            printf("*******************************\n");
            printf("error inside neighbour selected\n\n");
            exit(1);
         }
         return vof[4 * cOct + (k - 4)];
      }
      int prCell = octPrCell[nbCell];
      Real cc[3][3];
      Real mm1, mm2, mx, mz, mxp, mzp, V1, V2, alpha;
      int i, j;

      i = 1;
      j = 1;
      getCellNgbVOF(nbCell, cc);

      mm1 = cc[i - 1][j - 1] + 2. * cc[i - 1][j] + cc[i - 1][j + 1];
      mm2 = cc[i + 1][j - 1] + 2. * cc[i + 1][j] + cc[i + 1][j + 1];
      mxp = mm1 - mm2;
      mm1 = cc[i - 1][j - 1] + 2. * cc[i][j - 1] + cc[i + 1][j - 1];
      mm2 = cc[i - 1][j + 1] + 2. * cc[i][j + 1] + cc[i + 1][j + 1];
      mzp = mm1 - mm2;

      mx = fabs(mxp) + 1.e-50;
      mz = fabs(mzp) + 1.e-50;
      mm2 = DMAX(mx, mz);
      mx = mx / mm2;
      mz = mz / mm2;

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

      Real px[] = {0, 0.5, 0, 0.5};
      Real pz[] = {0, 0, 0.5, 0.5};

      if (mxp < 0)
      {
         // px = {0.5, 0, 0.5, 0};
         px[0] = px[2] = 0.5;
         px[1] = px[3] = 0;
      }

      if (mzp < 0)
      {
         // pz = {0.5, 0.5, 0, 0};
         pz[0] = pz[1] = 0.5;
         pz[2] = pz[3] = 0;
      }
      int pos = morton_lookup[dir][iCell % nbNumberOfOct] - 4;
      // rescaling and repositioning using alpha (n.r' = 2*(alpha - n.p))
      mm2 = 2 * (alpha - mx * px[pos] - mz * pz[pos]);
      return VOL2(mx, mz, mm2, 1);
   }
   return vof[nbCell];
}

void getChildVOF(int iOct, Real list[4], int prCell)
{ // gets vof of child if it exists if not it calculated from parent cell
  // setting iOct to 0 will force the child vof calculation

   // printf("childVOF *****************\n");
   // printf("iOct = %d\n", iOct);

   if (iOct)
   {
      list[0] = vof[4 * iOct + 0];
      list[1] = vof[4 * iOct + 1];
      list[2] = vof[4 * iOct + 2];
      list[3] = vof[4 * iOct + 3];
      return;
   }

   if (vof[prCell] == 0)
   {
      list[0] = 0;
      list[1] = 0;
      list[2] = 0;
      list[3] = 0;
      return;
   }

   if (vof[prCell] == 1)
   {
      list[0] = 1;
      list[1] = 1;
      list[2] = 1;
      list[3] = 1;
      return;
   }

   Real cc[3][3], mm1, mm2, mx, mz, mxp, mzp, V1, V2, alpha;
   int i, j;

   getCellNgbVOF(prCell, cc);

   i = 1;
   j = 1;
   mm1 = cc[i - 1][j - 1] + 2. * cc[i - 1][j] + cc[i - 1][j + 1];
   mm2 = cc[i + 1][j - 1] + 2. * cc[i + 1][j] + cc[i + 1][j + 1];
   mxp = mm1 - mm2;
   mm1 = cc[i - 1][j - 1] + 2. * cc[i][j - 1] + cc[i + 1][j - 1];
   mm2 = cc[i - 1][j + 1] + 2. * cc[i][j + 1] + cc[i + 1][j + 1];
   mzp = mm1 - mm2;

   mx = fabs(mxp) + 1.e-50;
   mz = fabs(mzp) + 1.e-50;
   mm2 = DMAX(mx, mz);
   mx = mx / mm2;
   mz = mz / mm2;

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

   Real px[] = {0, 0.5, 0, 0.5};
   Real pz[] = {0, 0, 0.5, 0.5};

   if (mxp < 0)
   {
      // px = {0.5, 0, 0.5, 0};
      px[0] = px[2] = 0.5;
      px[1] = px[3] = 0;
   }

   if (mzp < 0)
   {
      // pz = {0.5, 0.5, 0, 0};
      pz[0] = pz[1] = 0.5;
      pz[2] = pz[3] = 0;
   }

   // printf("alpha = %f\n", alpha);
   // printf("(mx, mz) = (%f %f)\n\n", mxp, mzp);

   for (int k = 0; k < 4; k++)
   {
      // rescaling and repositioning using alpha (n.r' = 2*(alpha - n.p))
      mm2 = 2 * (alpha - mx * px[k] - mz * pz[k]);
      list[k] = VOL2(mx, mz, mm2, 1);
   }
}

void getCellNgbVOF_6x6(int iOct, Real cc[][6])
{ // generates 6x6 grid of data from neighbours of octs
   int ngbOct, ngbOctCell, southOct, southOctCell, northOct, northOctCell, prCell;
   Real list[4];

   prCell = octPrCell[iOct];

   // center cells
   cc[2][2] = vof[4 * iOct + 0];
   cc[3][2] = vof[4 * iOct + 1];
   cc[2][3] = vof[4 * iOct + 2];
   cc[3][3] = vof[4 * iOct + 3];

   // west cell
   ngbOctCell = cellNb[0][prCell];
   ngbOct = cellChOct[ngbOctCell];
   // ngbOct = octNb[0][iOct];
   getChildVOF(ngbOct, list, ngbOctCell);
   cc[0][2] = list[0];
   cc[1][2] = list[1];
   cc[0][3] = list[2];
   cc[1][3] = list[3];

   // south-west cell
   southOctCell = cellNb[2][ngbOctCell];
   southOct = cellChOct[southOctCell];
   getChildVOF(southOct, list, southOctCell);
   cc[0][0] = list[0];
   cc[1][0] = list[1];
   cc[0][1] = list[2];
   cc[1][1] = list[3];

   // north-west cell
   northOctCell = cellNb[3][ngbOctCell];
   northOct = cellChOct[northOctCell];
   getChildVOF(northOct, list, northOctCell);
   cc[0][4] = list[0];
   cc[1][4] = list[1];
   cc[0][5] = list[2];
   cc[1][5] = list[3];

   // east cell
   ngbOctCell = cellNb[1][prCell];
   ngbOct = cellChOct[ngbOctCell];
   getChildVOF(ngbOct, list, ngbOctCell);
   cc[4][2] = list[0];
   cc[5][2] = list[1];
   cc[4][3] = list[2];
   cc[5][3] = list[3];

   // south-east cell
   southOctCell = cellNb[2][ngbOctCell];
   southOct = cellChOct[southOctCell];
   getChildVOF(southOct, list, southOctCell);
   cc[4][0] = list[0];
   cc[5][0] = list[1];
   cc[4][1] = list[2];
   cc[5][1] = list[3];

   // north-east cell
   northOctCell = cellNb[3][ngbOctCell];
   northOct = cellChOct[northOctCell];
   getChildVOF(northOct, list, northOctCell);
   cc[4][4] = list[0];
   cc[5][4] = list[1];
   cc[4][5] = list[2];
   cc[5][5] = list[3];

   // south cell
   ngbOctCell = cellNb[2][prCell];
   ngbOct = cellChOct[ngbOctCell];
   getChildVOF(ngbOct, list, ngbOctCell);
   cc[2][0] = list[0];
   cc[3][0] = list[1];
   cc[2][1] = list[2];
   cc[3][1] = list[3];

   // north cell
   ngbOctCell = cellNb[3][prCell];
   ngbOct = cellChOct[ngbOctCell];
   getChildVOF(ngbOct, list, ngbOctCell);
   cc[2][4] = list[0];
   cc[3][4] = list[1];
   cc[2][5] = list[2];
   cc[3][5] = list[3];
}

void printCellNgbVOF(int iCell)
{
   int i, j;
   Real cc[3][3];

   getCellNgbVOF(iCell, cc);
   printf("VOF around cell %d\n", iCell);

   for (j = 2; j >= 0; j--)
   {
      printf("%g %g %g\n", cc[0][j], cc[1][j], cc[2][j]);
   }
}