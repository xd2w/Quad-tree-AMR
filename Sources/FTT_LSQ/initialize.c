#include <stdio.h>
#include <stdlib.h>
#include "ftt.h"
#include "pfplib.h"

extern void readData(void);

void initialize(void)
{
  int restart =0;
  FILE *fp;

  ifetch("restart", &restart);

  readData();
  initMemory();
  if(restart == 0)
  {
    initFTT();
    refineFTT();
    //plotFTT(0);
    //plotSFC(0);
    //exit(1);
    octPropagation();
    establishNb();
    bcOctTree();
    establishNb();
    bcOctTree();
    establishNb();
    bcOctTree();
    establishNb();
    bcOctTree();
    establishNb();
    bcOctTree();
    checkNb();
  }
  else
  {
    fp = fopen("fttMesh","r");
    readFTT(fp);
    seekCell(0, .5736, 0.3854);
    establishNb();
    checkNb();
    drawPrCells(1861);
    drawNgbCells(1861);
    plotFTT(990);
    fclose(fp);
    exit(1);
  }
    checkOctTree();

  return;
}
