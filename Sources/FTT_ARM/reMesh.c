#include <stdio.h>
#include <stdlib.h>
#include "ftt.h"

extern void splitFlagCellsAtLevel(int level);
extern void oct_PrCellFlagAtLvel(int level);

// flags the interface cells and refines the flaged cells
void reMesh(int itNb)
{
  int level;

  // previous flag is cleared
  // setCellInt1DZero(cellFlag);
  setOctInt1DZero(octFlag);

  // interface cells are flaged
  printf("calling flag interface");
  flagInterfCells();

  //  for(level=maxLevel; level>maxLevel-1; level--)
  for (level = maxLevel; level > minLevel; level--)
  {
    setOctInt1DZeroAtLevel(octFlag, level);   // making sure no flag from previous level is carried over?
    cell_OctFlagAtLevel(level);               // flagging Octs of the flaged cell
    setCellInt1DZeroAtLevel(cellFlag, level); // wipe flag of all cells
    propagateOctFlagAtLevel(level);           // propagating the flag to neighbouring cells
    propagateOctFlagAtLevel(level);           // im gussing this was here to make refineing less sudden
    propagateOctFlagAtLevel(level);
    // propagateOctFlagAtLevel(level);
    splitFlagCellsAtLevel(level - 1); // split flagged cells
    binCollectionAtLevel(level);
    oct_PrCellFlagAtLvel(level);
    establishNb();
    checkNb();
  }
  // restrField(vof);
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
