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

// [direction][iCell % 4]
// if >3 its outside the Oct (subtract 4 to get index on outside neighbour)
static int morton_lookup[4][4] =
    {
        5, 0, 7, 2,
        1, 4, 3, 6,
        6, 7, 0, 1,
        2, 3, 4, 5};

static double est1[20] = {
    +8.809173 * 1e-3,
    -7.302646 * 1e-2,
    +1.476918 * 1e-1,
    +8.127531 * 1e-1,
    +1.662098 * 1e-1,
    -3.158622 * 1e-1,
    -7.899257 * 1e-2,
    -2.954707 * 1e-1,
    -8.268886 * 1e-1,
    +1.961225 * 1e-2,
    -1.107777 * 1e-1,
    -4.160499 * 1e-4,
    -1.507973 * 1e-2,
    -1.792648 * 1e-4,
    +8.268877 * 1e-1,
    +6.325443 * 1e-1,
    -6.501391 * 1e-2,
    +1.579836 * 1e-1,
    +2.944518 * 1e-6,
    +4.704595 * 1e-6,
};

static double est2[20] = {
    +8.744777 * 1e-3,
    -4.773088 * 1e-2,
    +1.911265 * 1e-1,
    +1.051815 * 100,
    +9.073978 * 1e-2,
    -3.656622 * 1e-1,
    -1.242422 * 1e-1,
    -3.824833 * 1e-1,
    -2.037791 * 1e-1,
    -9.380094 * 1e-1,
    -6.043060 * 1e-2,
    -1.001143 * 1e-3,
    -5.859334 * 1e-2,
    -3.742121 * 1e-4,
    +2.037529 * 1e-1,
    +7.331659 * 1e-1,
    +4.143641 * 1e-1,
    +2.482401 * 1e-1,
    +3.065368 * 1e-4,
    +1.256257 * 1e-4,
};

static double est3[20] = {
    +6.411956 * 1e-3,
    -1.217681 * 1e-1,
    +1.213353 * 1e-1,
    +1.091137 * 100,
    +3.268512 * 1e-1,
    -1.658780 * 1e-1,
    +9.856956 * 1e-2,
    -2.427756 * 1e-1,
    -1.522460 * 1e-1,
    -8.981913 * 1e-1,
    -2.178954 * 1e-1,
    -2.737967 * 1e-4,
    -2.349028 * 100,
    -2.355644 * 1e-5,
    +1.522527 * 1e-1,
    +3.322431 * 1e-1,
    +2.703116 * 1e-1,
    -1.982763 * 1e-1,
    +1.363379 * 1e-3,
    -3.048371 * 1e-7,
};

extern int maxLevel, minLevel, minIntfLevel, maxNumberOfOcts, maxNumberOfCells;
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
extern Real1D mxCell;
extern Real1D mzCell;
extern Real1D alphaCell;

/* corner coordinates of Cells */
extern Real1D xCell;
extern Real1D yCell;
#if (ocTree) /* 3D */
extern Real1D zCell;
#endif

/* physical quantities */
// velocity at cell corners
extern Real2D u;
extern Real2D v;
#if (ocTree) /* 3D */
extern Real2D w;
#endif
// velocity flux at cell faces
extern Real1D U;
extern Real1D V;
extern Real1D W;
// pressure at cell centers
extern Real1D p;
extern Real1D dive;
extern Real1D vof;
extern Real1D temp_vof;
extern Real1D work1, work2, work3;

// newly added variables for ellipse case

// coefficient of ellipse Axx + Bxy + Cyy = radius^2 [A, B, C]
extern Real init_VOF_coefs[3];

extern Real CFL;
extern Real global_dt;
extern Real t_total;
extern int kappaMode; // 0 - Barick, 1 - HF, 2 - Meier
extern Real runtime;

// for thoretical lines
extern int nThPoints;
extern Real1D xThPoints;
extern Real1D yThPoints;

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
extern void plotCellGrad(int iCell, FILE *fpGd);
extern void plotCellGradSmoothed(int iCell, FILE *fpGd);
extern void plotCellGrad_4x4(int iCell, FILE *fpGd);
extern void plotCellGradSmoothed_4x4(int iCell, FILE *fpGd);
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
extern void propagateFlagAtLevel(int dir, int level);
extern void propagateOctFlag(void);
extern void propagateOctFlagAtLevel(int level);
extern void setCellInt1DZero(Int1D val);
extern void setCellInt1DZeroAtLevel(Int1D val, int level);
extern void setOctInt1DZeroAtLevel(Int1D val, int level);
extern void setOctInt1DZero(Int1D val);
extern void reMesh(int itNb);
extern void balanceCellsAround(int iCell);

// extern Real computeVOF_ellipse(int iCell, int itNb);
// extern void refineFTT_ellipse(void);
extern void getCellNgbVOF_6x6(int iOct, Real cc[][6]);
extern void getChildVOF(int iOct, Real list[4], int prCell);
extern Real getVOF(int nbCell, int iCell, int dir);

extern void curvature_6x6(Real cc[][6], double kappas[4], Real dx, Real dy);
extern Real curvature_5x5(Real cc[][6], int ip, int jp, Real dx, Real dy);
extern void plotCurvatureAtLevel(int ndata, int level);
extern void plotCurvatureAtLeafCells(int ndata);
extern void show6x6VofGrid(int iCell);
extern void plotFTTInterfAtLevel(int ndata, int level);
extern void plotFTTAtLevel(int ndata, int level);
extern void flagInterfLeaves(void);
extern void calcWorksX(int iCell, Real vofVal, Real alpha, Real mx, Real mz, int invx, int invz, Real s1, Real s2);
// extern void calcWorksX(int iCell, Real vofVal, Real alpha, Real mx, Real mz, int invx, int invz, Real s1, Real s2, int scale1, int scale2, Real s11, Real s12, Real s21, Real s22);
extern void calcWorksY(int iCell, Real vofVal, Real alpha, Real mx, Real mz, int invx, int invz, Real s1, Real s2);
extern void calcWorksXFull(int iCell, Real s1, Real s2);
// extern void calcWorksXFull(int iCell, Real s1, Real s2, int scale1, int scale2, Real s11, Real s12, Real s21, Real s22);
extern void calcWorksYFull(int iCell, Real s1, Real s2);
extern void plotVOF(int ndata);
extern void plotVOFAtLevel(int ndata, int level);
extern void setTimeStep(void);
extern void smooth2ndOrd(Real from[][6], Real to[][4]);
extern void smoothUnifrom(Real from[][6], Real to[][4]);
extern void copycc(Real cc[][6], Real ccs[][4]);
extern Real equation_analytical_curvature(Real x, Real y, Real xc, Real yc);

extern Real kappaBarickALELike(int iCell, Real cc[][6]);
extern Real kappaBarickALELike_wider(int iCell, Real cc[][6]);
extern Real kappaMeier(int iCell);
extern Real kappaHF(int iCell, Real cc[][6]);

extern Real computeVX(Real x, Real y);
extern Real computeVY(Real x, Real y);
extern void computeVelocityAtLeaves(void);

extern void setPLICPramForAll(void);
extern void setPLICPramForOne(int iCell, Real cc[][3]);
extern Real getSubVOF(int iCell, int iLocal);

extern void splitCell_smart(int iCell);

extern void plotVOFContour(int ndata);
extern void plotCellGradAtIntf(int ndata);

extern void printcc3(Real cc[3][3]);
extern void printcc6(Real cc[6][6]);

extern void getCellNgbTempVOF(int iCell, Real cc[][3]);
extern void copyCellReal1D(Real1D from, Real1D to);

extern void getCellNgbVOF_unifrom(int iCell, Real cc[][3]);
extern void smooth1D(void);

extern void initTheoreticalInterf(void);
extern void plotTheoreticalInterf(int ndata);
extern void advTheoreticalInterf(void);

extern void constructHilbertIndex(void);

#endif
