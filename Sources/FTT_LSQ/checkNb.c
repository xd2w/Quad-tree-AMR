#include <stdio.h>
#include <stdlib.h>
#include "ftt.h"
#include "nrutil.h"

/* 
Check the neighbourhood among cells.
0 -> west neighbour, 1 -> east neighbour
2 -> south neighbour, 3 -> neighbour
Negative number neighbour means no neighbour.
*/

void checkNb(void)
{
  int iCell, iOct, iLv;

  for(iLv=1; iLv <= maxLevel; iLv++)
  {
    for(iCell=0; iCell<numberOfCells; iCell++)
    {
      iOct = iCell/cellNumberInOct; 
      if(cellType[iCell] == 0 && octLv[iOct]) checkCellNb(iCell);
    }  
  }  
 
}
/* check an interior cell neigbourhood */
void checkCellNb(int iCell)
{
  int iLv, iOct, prCell, ngbOct, ngbLv, ngbCell;

  iOct = iCell/4;
  iLv = octLv[iOct];
  prCell = octPrCell[iOct];
// west neighbour
   ngbCell = cellNb[0][iCell]; 
   if(ngbCell >= numberOfCells || ngbCell < -1)
   {
     printf("cell %d type %d west neighbour is %d\n", iCell, cellType[iCell],
             ngbCell);
     printf("ngbCell == maxNumberOfCells || ngbCell < -1\n");
     exit(1);
   }
   if(ngbCell>=0)
   {
     ngbOct = ngbCell/4;
     ngbLv = octLv[ngbOct];
     if(ngbLv != iLv && (ngbLv+1) != iLv)
     {
       printf("cell %d type %d west neighbour is %d\n", iCell, cellType[iCell],
               ngbCell);
       printf("!error in cell %d neighbourhood: ngbLv = %d iLv = %d\n",
               iCell, ngbLv, iLv);
       exit(1);
     }
     if(ngbLv == iLv)
     {
       if(cellNb[1][ngbCell] != iCell)
       {
         printf("error in cell neighbourhood: case 1\n"); exit(1);
       } 
     }
     else
     {
       if(cellNb[1][ngbCell] != prCell)
       {
         printf("error in cell neighbourhood: case 1\n"); exit(1);
       } 
     }  
   }
// east neighbour
   ngbCell = cellNb[1][iCell];
   if(ngbCell >= numberOfCells || ngbCell < -1)
   {
     printf("cell %d type %d east neighbour is %d\n", iCell, cellType[iCell],
             ngbCell);
     printf("ngbCell == maxNumberOfCells || ngbCell < -1\n");
     exit(1);
   }
   if(ngbCell>=0)
   {
     ngbOct = ngbCell/4;
     ngbLv = octLv[ngbOct];
     if(ngbLv != iLv && (ngbLv+1) != iLv)
     {
       printf("cell %d type %d east neighbour is %d\n", iCell, cellType[iCell],
               ngbCell);
       printf("!error in cell %d neighbourhood: ratio >2, ngbLv = %d iLv = %d\n",
               iCell, ngbLv, iLv);
       exit(1);
     }
     if(ngbLv == iLv)
     {
       if(cellNb[0][ngbCell] != iCell)
       {
         printf("error in cell east neighbourhood: case 1\n"); exit(1);
       }
     }
     else
     {
       if(cellNb[0][ngbCell] != prCell)
       {
         printf("error in cell east neighbourhood: case 1\n"); exit(1);
       }
     }

  }     

// south neighbour
   ngbCell = cellNb[2][iCell];
   if(ngbCell >= numberOfCells || ngbCell < -1)
   {
     printf("cell %d type %d south neighbour is %d\n", iCell, cellType[iCell],
             ngbCell);
     printf("ngbCell %d == maxNumberOfCells %d || ngbCell < -1\n", ngbCell, numberOfCells);
     printFTTCell(iCell);
     printFTTCell(ngbCell);
     drawNgbCells(iCell);
     exit(1);
   }
   if(ngbCell>=0)
   {
     ngbOct = ngbCell/4;
     ngbLv = octLv[ngbOct];
     if(ngbLv != iLv && (ngbLv+1) != iLv)
     {
       printf("cell %d type %d south neighbour is %d\n", iCell, cellType[iCell],
               ngbCell);
       printf("!error in cell %d neighbourhood: ngbLv = %d iLv = %d\n",
               iCell, ngbLv, iLv);
       exit(1);
     }
     if(ngbLv == iLv)
     {
       if(cellNb[3][ngbCell] != iCell)
       {
         printf("error in cell south neighbourhood: case 1\n"); exit(1);
       }
     }
     else
     {
       if(cellNb[3][ngbCell] != prCell)
       {
         printf("error in cell south neighbourhood: case 1\n"); exit(1);
       }
     }
  }
// north neighbour
   ngbCell = cellNb[3][iCell];
   if(ngbCell >= numberOfCells || ngbCell < -1)
   {
     printf("cell %d type %d north neighbour is %d\n", iCell, cellType[iCell],
             ngbCell);
     printf("ngbCell == maxNumberOfCells || ngbCell < -1\n");
     exit(1);
   }
   if(ngbCell>=0)
   {
     ngbOct = ngbCell/4;
     ngbLv = octLv[ngbOct];
     if(ngbLv != iLv && (ngbLv+1) != iLv)
     {
       printf("cell %d type %d north neighbour is %d\n", iCell, cellType[iCell],
               ngbCell);
       printf("!error in cell %d neighbourhood: ngbLv = %d iLv = %d\n",
               iCell, ngbLv, iLv);
       exit(1);
     }
     if(ngbLv == iLv)
     {
       if(cellNb[2][ngbCell] != iCell)
       {
         printf("error in cell north neighbourhood: case 1\n"); exit(1);
       }
     }
     else
     {
       if(cellNb[2][ngbCell] != prCell)
       {
         printf("error in cell north neighbourhood: case 1\n"); exit(1);
       }
     }

  }

}
 
