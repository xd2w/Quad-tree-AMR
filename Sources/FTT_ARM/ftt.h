/* Fully Threaded Tree Algorithms, Khokhlov, JCP, vol.143, pp519-543, 1998 */

#include "variable.h"
#ifndef FTT_H
#define FTT_H

/* ftt tree relations */
#if (ocTree) /* 3D */
#define cellNumberInOct 8
#define nbNumberOfOct 6
#else /* 2D */
#define cellNumberInOct 4
#define nbNumberOfOct 4
#endif

static int hilbert_map[4][4] =
    {
        0, 2, 3, 1, // H
        0, 1, 3, 2, // A
        3, 2, 0, 1, // B
        3, 1, 0, 2  // C
};

static int hilbert_production[4][4] =
    {
        1, 2, 0, 0, // H -> ABHH
        0, 1, 3, 1, // A -> HACA
        2, 0, 2, 3, // B -> BHBC
        3, 3, 1, 2  // C -> CCAB
};

extern int maxLevel, minLevel, maxNumberOfOcts, maxNumberOfCells;
extern int levelNumber, numberOfOcts, numberOfCells; // actual number

/* cell size at different level */
extern Real1D dxCell, dyCell, dzCell;
/* physical domain */
extern Real Lx, Ly, Lz;

/* oct variables */
extern Int1D octLv, octFlag;
extern Int1D octPrCell, octMark;
extern Int2D octNb;

/* cell variables */
extern Int1D cellChOct;
extern Int1D cellFlag, cellMark;
extern Int1D cellType;
extern Int2D cellNb;
extern Int1D cellHilb;

/* corner coordinates of Cells */
extern Real1D xCell;
extern Real1D yCell;
#if (ocTree) /* 3D */
extern Real1D zCell;
#endif

/* physical quantities */
// velocity at cell corners
// extern Real1D u;
// extern Real1D v;
// extern Real1D w;
// velocity flux at cell faces
// extern Real1D U;
// extern Real1D V;
// extern Real1D W;
// pressure at cell centers
// extern Real1D p;
// extern Real1D dive;
// extern Real1D vof;
// extern Real1D work1, work2, work3;
extern Real1D vx;
extern Real1D vy;
extern Real1D xCircle;
extern Real1D yCircle;
extern int numberOfCirclePoints;
extern int maxNumberOfCirclePoints;
extern Int1D hilbLeaves;

extern void initFTT(void);
extern void initMemory(void);
extern void computeCoord(void);
extern void computeCoordAtLevel(int Lv);
extern void initVOF(int it);
extern Real computeVOF(int iCell, int it);
extern Real polyArea(Real1D xList, Real1D yList, int count);
extern void initialize(void);
extern void splitCell(int iCell);
extern void splitFlagCells(void);
extern void balanceFTT(int balanceLevel, int levelDrop);
extern void balanceFTTxP(int iCell, int levelDrop);
extern void balanceFTTxN(int iCell, int levelDrop);
extern void balanceFTTyP(int iCell, int levelDrop);
extern void balanceFTTyN(int iCell, int levelDrop);
extern void refineFTT(void);
extern void plotSFC(int ndata);
extern void plotFTTCell(int iCell, FILE *fp);
extern void plotHilbertSFC(int ndata);
extern void plotFTTCellHilbert(int iCell, FILE *fsfc);
extern void plotFTTCellSFC(int iCell, FILE *fsfc);
extern void plotFlagCellsAtLevel(int ndata, int level);
extern void drawFTTOct(int iOct);
extern void plotOctMesh(int ndata);
extern void plotFlagOct(int iOct, FILE *fp);
extern void plotFlagOctsAtLevel(int ndata, int level);
extern void plotFlagOcts(int ndata);
extern void plotNgFlagOcts(int ndata);
extern void plotNgFlagOctsAtLevel(int ndata, int level);
extern void printFTTCell(int iCell);
extern void printFTTOct(int iOct);
extern void plotFTT(int ndata);
extern void saveFTT(FILE *fp);
extern void readFTT(FILE *fp);
extern void plotFlagFTT(int ndata);
extern void plotFlagCell(int iCell, FILE *fp);
extern void fttStatistics(void);
extern void cell_OctFlag(void);
extern void cell_OctFlagAtLevel(int level);
extern void oct_PrCellFlag(void);
extern void oct_PrCellFlagAtLevel(int level);
extern void restrFlagAtLevel(int lev);
extern void restrFlag(void);
extern void restrField(Real1D field);
extern void restrFieldAtLevel(Real1D field, int lev);
extern void markVOF(void);
extern void flagOct(void);
extern void garbageCollection(void);
extern void binCollection(void);
extern void binCollectionAtLevel(int level);
extern void drawNgbCells(int iCell);
extern void drawPrCells(int iCell);
extern void drawChCells(int iCell, FILE *fp);
extern int seekCell(int iCell, Real x, Real y);
extern void checkOctTree(void);
extern void establishNb(void);
extern void establishNbAtBase(void);
extern void establishNbAtLevel(int nghLv);
extern void checkNb(void);
extern void checkCellNb(int iCell);
extern void octPropagation(void);
extern void bcOctTree(void);
extern void bcOctTreeAtCell(int iCell, int direction);
extern void bcOctTreeAtWestCell(int iCell);
extern void bcOctTreeAtEstCell(int iCell);
extern void bcOctTreeAtSouthCell(int iCell);
extern void bcOctTreeAtNorthCell(int iCell);
extern void plotFTTInterf(int ndata);
extern void getCellNgbVOF(int iCell, Real cc[][3]);
extern void plotCellInterf(int iCell, FILE *fp);
extern void printCellNgbVOF(int iCell);
extern void octTreeXSwp(int iCell);
extern void octTreeYSwp(int iCell);
extern Real VOL2(Real mx, Real mz, Real alpha, Real b);
extern void plic(void);
extern void computeXVOF(void);
extern void computeYVOF(void);
extern void copyCellInt1D(Int1D from, Int1D to);
extern void copyOctInt1D(Int1D from, Int1D to);
extern void flagInterfCells(void);
extern void propagateFlag(int dir);
extern void propagateOctFlag(void);
extern void propagateOctFlagAtLevel(int level);
extern void setCellInt1DZero(Int1D val);
extern void setCellInt1DZeroAtLevel(Int1D val, int level);
extern void setOctInt1DZeroAtLevel(Int1D val, int level);
extern void setOctInt1DZero(Int1D val);
extern void reMesh(int itNb);

extern void initPotential(int itNb);
extern void initCircle();
extern void computePotential(int iCell, int itNb);
extern Real computeVX(Real x, Real y);
extern Real computeVY(Real x, Real y);
extern void propagateCircles(float dt);
extern void stream();
extern void plotCircle(int ndata);

#endif
