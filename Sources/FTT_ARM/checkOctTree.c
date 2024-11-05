#include <stdio.h>
#include <stdlib.h>
#include "ftt.h"
#include "nrutil.h"

void checkOctTree(void)
{
  int id, iCell, iLv, iOct, iFlag, octCell, chOct, prCell;



  for(iCell=0; iCell<numberOfCells; iCell++)
  {
     iOct = cellChOct[iCell];
     if(iOct >= numberOfOcts)
     {
       printf("error in Oct Tree at cell %d: cellChOct = %d >= %d\n", 
               iCell, iOct, numberOfOcts);
       exit(1);
     }
  } 
  for(iOct=1; iOct<numberOfOcts; iOct++)
  {
    iCell = octPrCell[iOct];
    if(iCell >= numberOfCells)
    {
      printf("error in Oct Tree at oct %d: octPrCell = %d >= %d\n",
             iOct, iCell, numberOfCells);
       exit(1);
    }
  }


  for(iCell=0; iCell<numberOfCells; iCell++)
  {
     iOct = cellChOct[iCell];
     if(iOct) // no leave cell
     {
       prCell = octPrCell[iOct];
       if(prCell != iCell)
       {
         printf("error in parent-child correspondence at cell %d\n", iCell);
         printf("child Oct is %d whose cell is %d\n", iOct, prCell);
         exit(1);
       }
     }
  }

  for(iOct=1; iOct<numberOfOcts; iOct++)
  {
    prCell = octPrCell[iOct];
    chOct = cellChOct[prCell];
    if(iOct != chOct)
    {
      printf("error in parent-child correspondence at oct %d\n", iOct);
      printf("parent cell is %d whose child oct is %d\n", prCell, chOct);
      exit(1);
    }
  }
}
