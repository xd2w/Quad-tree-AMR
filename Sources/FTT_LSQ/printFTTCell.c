#include <stdio.h>
#include <stdlib.h>
#include "ftt.h"
#include "nrutil.h"

void printFTTOct(int iOct)
{
  int prCell, octCell, flag;
  octCell = 4*iOct;
  flag = cellFlag[octCell];
  prCell = octPrCell[iOct];
  printf("FTT information Oct %d flag %d \n", iOct, flag);
  printf("Parent cell %d Parent Oct %d \n", 
                                prCell, prCell/4);
  printf("consists of cells: %d %d %d %d\n", octCell  , octCell+1, 
                                             octCell+2, octCell+3);
  printf("\n");
}
void printFTTCell(int iCell)
{
  int iOct, cLv, chCell, flag, type;
  Real xLeft, yLeft, xRight, yRight;
 
  iOct = iCell/cellNumberInOct; 
  cLv = octLv[iOct];
  flag = cellFlag[iCell];
  type = cellType[iCell];
  printf("FTT information on Cell %d\n", iCell);
  printf("Oct %d cLv %d flag %d type %d\n", iOct, cLv, flag, type);
  xLeft = xCell[iCell]; yLeft = yCell[iCell];
  xRight = xLeft+ dxCell[cLv];
  yRight = yLeft+ dyCell[cLv];
  printf("rectangle box\n");
  printf("%g %g\n", xLeft, yLeft);
  printf("%g %g\n", xRight, yLeft);
  printf("%g %g\n", xRight, yRight);
  printf("%g %g\n", xLeft, yRight);
  printf("%g %g\n", xLeft, yLeft);

  if(iOct)
  {
    printf("Cell %d's parent cell is %d\n", iCell, octPrCell[iOct]); 
  }
  if(cellChOct[iCell])
  {
    printf("Cell %d's child Oct is Oct %d\n", iCell, cellChOct[iCell]);
     chCell = 4*cellChOct[iCell];
    printf("Cell %d's children Cells are %d %d %d %d\n", 
            iCell, chCell, chCell+1, chCell+2, chCell+3);
  }
  else
  {
    printf("Cell %d has no child\n", iCell);
  }
  printf("\n");

}

