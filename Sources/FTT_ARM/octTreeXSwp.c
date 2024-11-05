#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "ftt.h"
#include "nrutil.h"

#define MAX(x, y)       ((x) > (y) ? (x) : (y))
#define MIN(x, y)       ((x) < (y) ? (x) : (y))

/* ------------------------------------------------------------------- */
/* u is the flux on the cell faces and cc volume fraction */
/* work1, work2, work3 temporal variables providing working allocations */
void octTreeXSwp(int iCell)
{
  int i,j,inv;
  Real mx,mz,alpha,s1,s2,mm1,mm2,V1,V2;
  Real cc[3][3];

  getCellNgbVOF(iCell, cc);
  i = 1; j = 1; 

/*
   connect to your code 
      s1 = u[i][j]*tau/(deltax*deltay);
      s2 = u[i+1][j]*tau/(deltax*deltay);
*/
   s1 = 0.1; s2 = 0.1;

      if(cc[i][j] == 0.0)
      {
	work1[iCell]=0.0;
	work2[iCell]=0.0;
	work3[iCell]=0.0;
      }
      else if(cc[i][j] == 1.0){
	work1[iCell] = DMAX(-s1,0.0);
	work2[iCell] = 1.0 - DMAX(s1,0.0) + DMIN(s2,0.0);
	work3[iCell] = DMAX(s2,0.0);
      }
      else
      {

	/* normal to the interface */
	mm1 = cc[i-1][j-1] + 2.*cc[i-1][j] + cc[i-1][j+1];
	mm2 = cc[i+1][j-1] + 2.*cc[i+1][j] + cc[i+1][j+1];
	mx = mm1 - mm2;
	mm1 = cc[i-1][j-1] + 2.*cc[i][j-1] + cc[i+1][j-1];
	mm2 = cc[i-1][j+1] + 2.*cc[i][j+1] + cc[i+1][j+1];
	mz = mm1 - mm2;
	
	inv = 0;
	if(mx < 0.)
        {
	  mm1 = -s1;
	  s1 = -s2;
	  s2 = mm1;
	  mx = -mx;
	  inv = 1;
	}
	mx = mx +1.e-50; 
	mz = fabs(mz) + 1.e-50;
	mm2 = DMAX(mx,mz);
	mx = mx/mm2;
	mz = mz/mm2;

	/* get alpha to determine the equation of the interface */
	mm1 = DMIN(mx,mz);
        /* compute the two critical volume fraction */
	V1 = 0.5*mm1;
	V2 = 1.0 - V1;
	if (cc[i][j] <= V1)
	  alpha = sqrt(2.0*cc[i][j]*mm1);
	else if (cc[i][j] <= V2) 
	  alpha = cc[i][j] + 0.5*mm1;
	else
	  alpha = mm1 + 1.0 - sqrt(2.0*(1.0-cc[i][j])*mm1);

	/* the new equation of the interface after advection */
	mx = mx/(1.0-s1+s2);
	alpha = alpha + mx*s1;
 
	/* calculate works */
	mm1 = DMAX(-s1,0.0);
	mm2 = alpha + mx*mm1;
	work1[iCell] = VOL2(mx,mz,mm2,mm1);

	mm1 = 1.0 - DMAX(s1,0.0) + DMIN(s2,0.0);
	mm2 = alpha - mx*DMAX(s1,0.0);
	work2[iCell] = VOL2(mx,mz,mm2,mm1);

	mm1 = DMAX(s2,0.0);
	mm2 = alpha - mx;
	work3[iCell] = VOL2(mx,mz,mm2,mm1);

	/* symmetry conditions */
	if(inv == 1) 
        {
	  mm1 = work1[iCell];
	  work1[iCell] = work3[iCell];
	  work3[iCell] = mm1;
	}
      }

    return;
}

