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
void getCellNgbVOF(int iCell, Real cc[][3])
{
   int ngbCell, southCell, northCell;

   // TODO fix this

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

void getChildVOF(int iOct, Real list[4], int prCell)
{ // gets vof of child if it exists if not it calculated from parent cell
   // TODO fix this not all cell has vof of (0, 1)

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
   mm2 = DMAX(mxp, mzp);
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
/*
void calArea()
{
   Real px[] = {0, 0, 0.5, 0.5, 0};
   Real py[] = {0, 0.5, 0.5, 0, 0};

   Real pointsx[5];
   Real pointsy[5];
   int points_counter;

   Real x_offset = 0;
   Real y_offset = 0;

   bool x_dir = 1;

   Real area, x0, y0, x1, y1;

   for (i = 0; i < 4; i++)
   {
      for (int pos = 0; pos < 4; pos++)
      {
         x0 = px[i] + px[pos];
         y0 = py[i] + px[pos];

         if (mx * x0 + mz * y0 < alpha)
         {
            //   points.append([x0, y0])
            pointsx[points_counter] = x0;
            pointsy[points_counter] = y0;
            points_counter++;
         }

         x1 = px[i + 1] + px[pos];
         y1 = py[i + 1] + px[pos];

         if (x_dir)
         {
            if (y1 < y0)
            {
               mm2 = y0;
               y0 = y1;
               y1 = mm2;
            }
            // checking if line crossed along y
            mm2 = (alpha - mx * x0) / mz;
            if ((y0 < mm2) && (mm2 < y1))
            {
               // points.append([x0, yp])
               pointsx[points_counter] = x0;
               pointsy[points_counter] = y0;
               points_counter++;
            }
            x_dir = !x_dir; // flip direction
         }
         else
         {
            if (x1 < x0)
            {
               mm2 = x0;
               x0 = x1;
               x1 = mm2;
            }
            // checking if line crossed along x
            mm2 = (alpha - mz * y0) / mx;
            if ((x0 < mm2) && (mm2 < x1))
            {
               // points.append([xp, y0])
               pointsx[points_counter] = x0;
               pointsy[points_counter] = y0;
               points_counter++;
            }
            x_dir = !x_dir; // flip direction
         }
      }

      list[i] = fabs(polyArea(pointsx, pointsy, points_counter));
   }

   // fixing ordering to be Z-ordering
   mm2 = list[1];
   list[1] = list[3];
   list[3] = mm2;

   mm2 = list[3];
   list[3] = list[2];
   list[2] = mm2;
}
*/