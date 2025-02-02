#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "pfplib.h"
#include "ftt.h"
#include "nrutil.h"
#include "box.h"
#include "pfplib.h"


Real computeVOF(int iCell, int itNb)
{
  int iOct, cLv, i, count;
  Real xc, yc, radius;
  Real dx, dy, xcoor[5], ycoor[5];
  Real xList[8], yList[8], color[5], lambda, fraction;
  xc = 0.75;
  yc = 0.5;
  radius = 0.1;
  dfetch("xc", &xc);
  dfetch("yc", &yc);
  dfetch("radius", &radius);
  /* position of the circle at itNb interations */
  xc +=.0075*itNb;
  yc +=.005*itNb;

  printf("\ncircle center at (%g,%g), radius = %g itNb %d\n", xc, yc, radius, itNb);

  iOct = iCell/cellNumberInOct;
  cLv = octLv[iOct];

  dx = dxCell[cLv]; 
  dy = dyCell[cLv]; 
  printf("dx %g dy %g\n", dx, dy);

  xcoor[0] = xcoor[4] = xCell[iCell];
  ycoor[0] = ycoor[4] = yCell[iCell];
  xcoor[1] = xCell[iCell] + dx;
  ycoor[1] = yCell[iCell];
  xcoor[2] = xCell[iCell] + dx;
  ycoor[2] = yCell[iCell] + dy;
  xcoor[3] = xCell[iCell];
  ycoor[3] = yCell[iCell] + dy;
  for(i=0; i<5; i++)
  {
    printf("box %d %g %g \n", i, xcoor[i], ycoor[i]);
    color[i] = (xcoor[i]-xc)*(xcoor[i]-xc) 
             + (ycoor[i]-yc)*(ycoor[i]-yc) - radius*radius; 
  } 
  if(color[0]<0.0)
   if(color[1]<0.0)
    if(color[2]<0.0)
     if(color[3]<0.0) return 1.0;
  if(color[0]>0.0)
   if(color[1]>0.0)
    if(color[2]>0.0)
     if(color[3]>0.0) return 0.0;

  count = 0;
  printf("!!!! count %d\n", count); 
  for(i=0; i<5; i++)
  {
    printf("box %d %g %g color %g\n", i, xcoor[i], ycoor[i], color[i]);
  } 
  for(i=0; i<4; i++)
  {
    if(color[i] <=0.0 )
    {
      xList[count] = xcoor[i];
      yList[count] = ycoor[i];
      count++;
      printf("count %d %g %g color %g\n", count, xcoor[i], ycoor[i], color[i]);
    }
    if(color[i]*color[i+1] < 0.0)
    {
      // interpolation 
      printf("color2 %d %g %g %g\n", i, color[i], color[i+1], color[i]*color[i+1]);
      lambda = fabs(color[i])/(fabs(color[i]) + fabs(color[i+1])+1e-50);
      xList[count] = (1.0-lambda)*xcoor[i]+ lambda*xcoor[i+1];
      yList[count] = (1.0-lambda)*ycoor[i]+ lambda*ycoor[i+1];
      count++;
    }
    
  }
  printf("COUNT %d\n", count); 
  //if(iCell==55 && itNb==0) exit(1);
  if(count) 
  {
    xList[count] = xList[0];
    yList[count] = yList[0];
  }
  if(count)
  {
    cellFlag[iCell] = 1;
    fraction = polyArea(xList, yList, count);
if(iCell==55 && itNb==0)
{
    printf("count %d\n", count);
    for(i=0; i<=count; i++)
    {
      printf("%g %g\n", xList[i],yList[i]);
    }
    printf("\n");
  //  exit(1); 
}
    printf("box\n");
    for(i=0; i<=4; i++)
    {
      printf("%g %g\n", xcoor[i],ycoor[i]);
    }
  }
  else
  {
    printf("box outside circle\n"); exit(1);
  }
//  exit(1);
  fraction /= dx*dy;
//  printf("fraction %g\n", fraction);
   if(fraction >= 1.0)
   {
     printf("error in polyArea %g iCell %d itNb %d cLv %d\n", fraction,
      iCell, itNb, cLv); exit(1);
   }
  return fraction;
}

Real polyArea(Real1D xList, Real1D yList, int count)
{
   int i;
   Real x0, x1, x2;
   Real y0, y1, y2, area;
   area = 0.0;
   x0 = xList[0];
   y0 = yList[0];
 //  printf("count %d \n", count);
   for(i=1; i<count-1; i++)
   {
//     printf("i %d\n", i);
     x1 = xList[i]; y1 = yList[i];
     x2 = xList[i+1]; y2 = yList[i+1];
     area += 0.5*((x1-x0)*(y2-y0)+(y1-y0)*(x0-x2));
   }
   return area;
}


Real computeVOF_ellipse(int iCell, int itNb, Real a, Real b, Real c)
{
  int iOct, cLv, i, count;
  Real xc, yc, radius;
  Real dx, dy, xcoor[5], ycoor[5];
  Real xList[8], yList[8], color[5], lambda, fraction;
  xc = 0.75;
  yc = 0.5;
  radius = 0.1;
  dfetch("xc", &xc);
  dfetch("yc", &yc);
  dfetch("radius", &radius);
  /* position of the circle at itNb interations */
  xc +=.0075*itNb;
  yc +=.005*itNb;

  printf("\ncircle center at (%g,%g), radius = %g itNb %d\n", xc, yc, radius, itNb);

  iOct = iCell/cellNumberInOct;
  cLv = octLv[iOct];

  dx = dxCell[cLv]; 
  dy = dyCell[cLv]; 
  printf("dx %g dy %g\n", dx, dy);

  xcoor[0] = xcoor[4] = xCell[iCell];
  ycoor[0] = ycoor[4] = yCell[iCell];
  xcoor[1] = xCell[iCell] + dx;
  ycoor[1] = yCell[iCell];
  xcoor[2] = xCell[iCell] + dx;
  ycoor[2] = yCell[iCell] + dy;
  xcoor[3] = xCell[iCell];
  ycoor[3] = yCell[iCell] + dy;
  for(i=0; i<5; i++)
  {
    printf("box %d %g %g \n", i, xcoor[i], ycoor[i]);
    color[i] = a*(xcoor[i]-xc)*(xcoor[i]-xc) 
             + b*(xcoor[i]-xc)*(ycoor[i]-yc)
             + c*(ycoor[i]-yc)*(ycoor[i]-yc) - radius*radius; 
  } 
  if(color[0]<0.0)
   if(color[1]<0.0)
    if(color[2]<0.0)
     if(color[3]<0.0) return 1.0;
  if(color[0]>0.0)
   if(color[1]>0.0)
    if(color[2]>0.0)
     if(color[3]>0.0) return 0.0;

  count = 0;
  printf("!!!! count %d\n", count); 
  for(i=0; i<5; i++)
  {
    printf("box %d %g %g color %g\n", i, xcoor[i], ycoor[i], color[i]);
  } 
  for(i=0; i<4; i++)
  {
    if(color[i] <=0.0 )
    {
      xList[count] = xcoor[i];
      yList[count] = ycoor[i];
      count++;
      printf("count %d %g %g color %g\n", count, xcoor[i], ycoor[i], color[i]);
    }
    if(color[i]*color[i+1] < 0.0)
    {
      // interpolation 
      printf("color2 %d %g %g %g\n", i, color[i], color[i+1], color[i]*color[i+1]);
      lambda = fabs(color[i])/(fabs(color[i]) + fabs(color[i+1])+1e-50);
      xList[count] = (1.0-lambda)*xcoor[i]+ lambda*xcoor[i+1];
      yList[count] = (1.0-lambda)*ycoor[i]+ lambda*ycoor[i+1];
      count++;
    }
    
  }
  printf("COUNT %d\n", count); 
  //if(iCell==55 && itNb==0) exit(1);
  if(count) 
  {
    xList[count] = xList[0];
    yList[count] = yList[0];
  }
  if(count)
  {
    cellFlag[iCell] = 1;
    fraction = polyArea(xList, yList, count);
if(iCell==55 && itNb==0)
{
    printf("count %d\n", count);
    for(i=0; i<=count; i++)
    {
      printf("%g %g\n", xList[i],yList[i]);
    }
    printf("\n");
  //  exit(1); 
}
    printf("box\n");
    for(i=0; i<=4; i++)
    {
      printf("%g %g\n", xcoor[i],ycoor[i]);
    }
  }
  else
  {
    printf("box outside circle\n"); exit(1);
  }
//  exit(1);
  fraction /= dx*dy;
//  printf("fraction %g\n", fraction);
   if(fraction >= 1.0)
   {
     printf("error in polyArea %g iCell %d itNb %d cLv %d\n", fraction,
      iCell, itNb, cLv); exit(1);
   }
  return fraction;
}