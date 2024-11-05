#include "variable.h"
#ifndef BOX_H
#define BOX_H
struct point
{
  Real x, y;
#if (ocTree) /* 3D */
  Real z;
#endif
};

struct box
{
  struct point pt1, pt2;
};

extern void printPoint(struct point p);
extern void printBox(struct box b);
extern int ptWithinBox(struct point p, struct box b);
extern int ptListIntersectBox(struct point *pList, int nPoint, struct box b);

#endif
