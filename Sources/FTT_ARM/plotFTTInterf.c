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