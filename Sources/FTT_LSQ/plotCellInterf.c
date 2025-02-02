#define MAX(x, y)       ((x) > (y) ? (x) : (y))
#define MIN(x, y)       ((x) < (y) ? (x) : (y))
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "ftt.h"
#include "nrutil.h"

void plotCellInterf(int iCell, FILE *fp) 
{
  int i,j,invx,invz,imax,jmax;
  Real mx,mz,mm1,mm2,V1,V2,alpha;
  int iOct, cLv;
  Real cc[3][3], xcoor, ycoor, dx, dy;
  float x1,z1,x2,z2;
  float umax, vmax;
  float *interf;

  getCellNgbVOF(iCell, cc);
  xcoor = xCell[iCell];
  ycoor = yCell[iCell];
  iOct = iCell/cellNumberInOct;
  cLv = octLv[iOct];
  dx = dxCell[cLv];
  dy = dyCell[cLv];

  i=1; j=1;
            /* PLIC: normal to the interface */
            mm1 = cc[i-1][j-1] + 2.0*cc[i-1][j] + cc[i-1][j+1];
            mm2 = cc[i+1][j-1] + 2.0*cc[i+1][j] + cc[i+1][j+1];
            mx = mm1 - mm2;
            mm1 = cc[i-1][j-1] + 2.0*cc[i][j-1] + cc[i+1][j-1];
            mm2 = cc[i-1][j+1] + 2.0*cc[i][j+1] + cc[i+1][j+1];
            mz = mm1 - mm2;
 // symmetry
            invx = 0;
            if (mx < 0)
              {
                mx = -mx;
                invx = 1;
              }
            invz = 0;
            if (mz < 0)
              {
                mz = -mz;
                invz = 1;
              }
            mx += 1e-50;
            mz += 1e-50;
            mm2 = MAX(mx, mz);
            mx /= mm2;
            mz /= mm2;

            mm1 = MIN(mx, mz);
            V1 = 0.5 * mm1;
            V2 = 1 - V1;
          // determine alpha from volume fraction 
            if (cc[i][j] < V1)
            {
                alpha = sqrt(2.0*cc[i][j]*mm1);
                x1 = alpha/mx; z1 = 0.0;
                x2 = 0.0;      z2 = alpha/mz;
                if (x1 > 1.0) x1 = 1.0;
                if (z2 > 1.0) z2 = 1.0;
            } 
            else if(cc[i][j] < V2)
            { 
               alpha = cc[i][j] + 0.5 * mm1;
               if(mx > mz)
               {
                   x1 = alpha/mx; x2 = (alpha - mz)/mx;
                   z1 = 0.0; z2 = 1.0;
               }
               else
               {
                   x1 = 0.0; x2 = 1.0;
                   z1=  alpha/mz;  z2 = (alpha - mx)/mz;
               }
            }
            else
            { 
                alpha = sqrt(2.0*(1.0-cc[i][j])*mm1);
                x1 = 1.0 - alpha/mx; z1 = 1.0;
                x2 = 1.0;      z2 = 1.0 - alpha/mz;
                if (x1 < 0.0) x1 = 0.0;
                if (z2 < 0.0) z2 = 0.0;
            }
        // use symmetry
            if(invx != 0)
            {
                x1 = 1.0 - x1; x2 = 1.0 - x2;
            }
            if(invz != 0)
            {
                z1 = 1.0 - z1; z2 = 1.0 -z2;
            }
            x1 = xcoor+dx*x1;
            z1 = ycoor+dy*z1;
            x2 = xcoor+dx*x2;
            z2 = ycoor+dy*z2;
           fprintf(fp,"%f %f \n",x1, z1);
           fprintf(fp,"%f %f \n",x2, z2);
           fprintf(fp,"\n");

    return ;
}
