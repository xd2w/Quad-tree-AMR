#include <stdio.h>
#include <stdlib.h>
#include "box.h"

int ptWithinBox(struct point p, struct box b)
{
  return b.pt1.x <= p.x && p.x <= b.pt2.x &&
         b.pt1.y <= p.y && p.y <= b.pt2.y
#if (ocTree) /* 3D */
      && b.pt1.z <= p.z && p.z <= b.pt2.z
#endif
;
}

int ptListIntersectBox(struct point *pList, int nPoint, struct box b)
{
  int ipt, true;
  true = 0;
  for(ipt=0; ipt<nPoint; ipt++)
  {
    if(ptWithinBox(pList[ipt], b)) true ++;
  }
  return true;
}
void printPoint(struct point p)
{
#if (ocTree) /* 3D */
  printf("point: %g %g %g\n", p.x, p.y, p.z);
#else
  printf("point: %g %g\n", p.x, p.y);
#endif
}
void printBox(struct box b)
{
  printf("Box\n");
  printPoint(b.pt1);
  printPoint(b.pt2);
}

