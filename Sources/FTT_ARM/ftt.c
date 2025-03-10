#include "variable.h"

/* NOTE BOOK (V), page 52
  fully threaded tree data structure:
  each Oct has four cells, which share a same parent cell.
  cellChOct --> pointing to the child Oct of a cell
  octPrCell --> pointing to the parent cell of an Oct
  octPrCell[0] is not defined
  When cellChOct = 0, cell does not have children.
  When cellChOct !=0, cell has 4 children cells which
  located between cellChOct*4 and cellChOct*4+3
*/

/* ftt tree relations */
#if (ocTree) /* 3D */
#define cellNumberInOct 8
#define nbNumberOfOct 6
#else /* 2D */
#define cellNumberInOct 4
#define nbNumberOfOct 4
#endif
int maxLevel, minLevel, maxNumberOfOcts, maxNumberOfCells;
int levelNumber, numberOfOcts, numberOfCells; // actual number

/* cell size at different level */
Real1D dxCell, dyCell, dzCell;
/* physical domain */
Real Lx, Ly, Lz;

/* oct variables */
Int1D octLv, octFlag;
Int1D octPrCell, octMark;
Int2D octNb;

/* cell variables */
Int1D cellChOct;
Int1D cellFlag, cellMark;
Int1D cellType;
Int2D cellNb;
Int1D cellHilb;

/* corner coordinates  of Cells*/
Real1D xCell;
Real1D yCell;
#if (ocTree) /* 3D */
Real1D zCell;
#endif

/* physical quantities */
// velocity at cell corners
Real1D u;
Real1D v;
Real1D w;
// velocity flux at cell faces
Real1D U;
Real1D V;
Real1D W;
// pressure at cell centers
Real1D p;
Real1D dive;
Real1D vof;
Real1D work1, work2, work3;

// newly added variables for ellipse case
Real init_VOF_coefs[3];
Real1D temp_vof;
Real CFL;
Real global_dt;
