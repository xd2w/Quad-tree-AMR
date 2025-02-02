#include <stdio.h>
#include <stdlib.h>
#include "ftt.h"

extern void splitFlagCellsAtLevel(int level);
extern void oct_PrCellFlagAtLvel(int level);

void reMesh(int itNb)
{
  int level;
 
  setCellInt1DZero(cellFlag);
  setOctInt1DZero(octFlag);
  flagInterfCells();
//  for(level=maxLevel; level>maxLevel-1; level--)
  for(level=maxLevel; level>minLevel; level--)
  {
    setOctInt1DZeroAtLevel(octFlag, level);
    cell_OctFlagAtLevel(level);
    setCellInt1DZeroAtLevel(cellFlag, level);
    propagateOctFlagAtLevel(level);
   // propagateOctFlagAtLevel(level);
    splitFlagCellsAtLevel(level-1);
    binCollectionAtLevel(level);
    oct_PrCellFlagAtLvel(level);
    establishNb();
    checkNb();
  }
  restrField(vof);
/*
  drawNgbCells(510);
  printCellNgbVOF(510);
  level=maxLevel-1;
  setOctInt1DZeroAtLevel(octFlag, level);
  cell_OctFlagAtLevel(level);
  plotFlagOctsAtLevel(0, level);
    propagateOctFlagAtLevel(level);
  plotFlagOctsAtLevel(19, level);
    plotFlagCellsAtLevel(19, level-1);
    propagateOctFlagAtLevel(level);
    splitFlagCellsAtLevel(level-1);
  plotFlagOctsAtLevel(20, level);
    plotFlagCellsAtLevel(22, level-1);
    plotNgFlagOctsAtLevel(20, level);
    binCollectionAtLevel(level);
    plotFlagCellsAtLevel(23, level-1);
*/
}
