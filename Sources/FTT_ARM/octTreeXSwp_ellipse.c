#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "ftt.h"
#include "nrutil.h"

#include "pfplib.h"

#define MAX(x, y) ((x) > (y) ? (x) : (y))
#define MIN(x, y) ((x) < (y) ? (x) : (y))
#define SWAP(a, b) \
	{              \
		temp = a;  \
		a = b;     \
		b = temp;  \
	}

/* ------------------------------------------------------------------- */
/* u is the flux on the cell faces and cc volume fraction */
/* work1, work2, work3 temporal variables providing working allocations */
void octTreeXSwp(int iCell)
{
	int i, j, invx, invz, iLv;
	Real mx, mz, alpha, s1, s2, mm1, mm2, V1, V2;
	Real cc[3][3];

	getCellNgbVOF(iCell, cc);
	i = 1;
	j = 1;

	/*
		  s1 = (u[i][j]/dx)*dt(cfl);
		  s2 = (u[i+1][j]/dx)*dt(cfl);
	*/
	// Real ux = 0.8; // change to actual u later
	// dfetch("ux", &ux);
	Real ux1, ux2, x, y, dx, dy;

	iLv = octLv[iCell / 4];

	x = xCell[iCell];
	y = yCell[iCell];
	dx = dxCell[iLv];
	dy = dyCell[iLv];

	ux1 = (computeVX(x, y) + computeVX(x, y + dy)) * 0.5;
	ux2 = (computeVX(x + dx, y) + computeVX(x + dx, y + dy)) * 0.5;

	// # TODO: do this with u[iCell] to be flexible
	// also use that to do vof tranfer of larger to small with higher res v

	s1 = (ux1 / dxCell[iLv]) * global_dt;
	s2 = (ux2 / dxCell[iLv]) * global_dt;

	// printf("\ns1, s2 = %f\n", s1);

	// TODO adjust s1 and s2 depending on cfl
	// same s1 s2 for different sized cells makes no sense

	if (cc[i][j] == 0.0)
	{
		return;
	}
	else if (cc[i][j] == 1.0)
	{
		calcWorksXFull(iCell, s1, s2);
	}
	else
	{

		/* normal to the interface */
		mm1 = cc[i - 1][j - 1] + 2. * cc[i - 1][j] + cc[i - 1][j + 1];
		mm2 = cc[i + 1][j - 1] + 2. * cc[i + 1][j] + cc[i + 1][j + 1];
		mx = mm1 - mm2;
		mm1 = cc[i - 1][j - 1] + 2. * cc[i][j - 1] + cc[i + 1][j - 1];
		mm2 = cc[i - 1][j + 1] + 2. * cc[i][j + 1] + cc[i + 1][j + 1];
		mz = mm1 - mm2;

		invx = 0;
		if (mx < 0.)
		{
			mm1 = -s1;
			s1 = -s2;
			s2 = mm1;
			mx = -mx;
			invx = 1;
		}

		invz = 0;
		if (mz < 0.0)
		{
			invz = 1;
		}

		mx = mx + 1.e-50;
		mz = fabs(mz) + 1.e-50;
		mm2 = DMAX(mx, mz);
		mx = mx / mm2;
		mz = mz / mm2;

		/* get alpha to determine the equation of the interface */
		mm1 = DMIN(mx, mz);
		/* compute the two critical volume fraction */
		V1 = 0.5 * mm1;
		V2 = 1.0 - V1;
		if (cc[i][j] <= V1)
			alpha = sqrt(2.0 * cc[i][j] * mm1);
		else if (cc[i][j] <= V2)
			alpha = cc[i][j] + 0.5 * mm1;
		else
			alpha = mm1 + 1.0 - sqrt(2.0 * (1.0 - cc[i][j]) * mm1);

		/* the new equation of the interface after advection */
		mx = mx / (1.0 - s1 + s2);
		alpha = alpha + mx * s1;

		calcWorksX(iCell, vof[iCell], alpha, mx, mz, invx, invz, s1, s2);
	}

	return;
}

void calcWorksX(int iCell, Real vofVal, Real alpha, Real mx, Real mz, int invx, int invz, Real s1, Real s2)
{
	int leftNb, rightNb, temp;
	Real V1, V3, mm1, mm2;
	int ltop, lbot, rtop, rbot;

	// 2 3 < ltop  rtop > 2 3
	// 0 1 < lbot  rbot > 0 1

	// the child position if the destination
	ltop = 3;
	lbot = 1;
	rtop = 2;
	rbot = 0;

	leftNb = cellNb[0][iCell];
	rightNb = cellNb[1][iCell];
	if (leftNb == 0 || rightNb == 0)
	{
		printf("***************************************\n");
		printf("Error in calcWorksX: invalid neighbours\n\n");
		exit(1);
	}

	if (invx)
	{ // inverting left and right
		SWAP(ltop, rtop)
		SWAP(lbot, rbot)
		SWAP(leftNb, rightNb)
	}

	if (invz)
	{ // inverting top and bottom
		SWAP(ltop, lbot)
		SWAP(rtop, rbot)
	}

	if (vofVal >= 1)
	{
		return;
	}

	// complicated case
	if (0 < vofVal && vofVal < 1)
	{
		mm1 = 1.0 - DMAX(s1, 0.0) + DMIN(s2, 0.0);
		mm2 = alpha - mx * DMAX(s1, 0.0);
		temp_vof[iCell] += VOL2(mx, mz, mm2, mm1);

		if (octLv[leftNb / 4] == octLv[iCell / 4] && cellChOct[leftNb] != 0)
		{ // if neighbour is smaller
			// printf("\nworks\n");
			// top
			mm1 = DMAX(-s1, 0.0);
			mm2 = 2 * (alpha + mx * mm1 - mz * 0.5);
			temp_vof[4 * cellChOct[leftNb] + ltop] += VOL2(mx, mz, mm2, 2 * mm1);
			// bottom
			mm1 = DMAX(-s1, 0.0);
			mm2 = 2 * (alpha + mx * mm1);
			temp_vof[4 * cellChOct[leftNb] + lbot] += VOL2(mx, mz, mm2, 2 * mm1);
		}
		else
		{
			mm1 = DMAX(-s1, 0.0);
			mm2 = alpha + mx * mm1;
			if (octLv[leftNb / 4] < octLv[iCell / 4])
			{
				temp_vof[leftNb] += 0.25 * VOL2(mx, mz, mm2, mm1);
			}
			else
			{
				temp_vof[leftNb] += VOL2(mx, mz, mm2, mm1);
			}
		}

		if (octLv[rightNb / 4] == octLv[iCell / 4] && cellChOct[rightNb] != 0)
		{ // if neighbour is smaller
			// printf("\nworks\n");
			mm1 = DMAX(s2, 0.0);
			mm2 = 2 * (alpha - mx - mz * 0.5);
			temp_vof[4 * cellChOct[rightNb] + rtop] += VOL2(mx, mz, mm2, 2 * mm1);

			mm1 = DMAX(s2, 0.0);
			mm2 = 2 * (alpha - mx);
			temp_vof[4 * cellChOct[rightNb] + rbot] += VOL2(mx, mz, mm2, 2 * mm1);
		}
		else
		{
			mm1 = DMAX(s2, 0.0);
			mm2 = (alpha - mx);
			if (octLv[rightNb / 4] < octLv[iCell / 4])
			{
				temp_vof[rightNb] += 0.25 * VOL2(mx, mz, mm2, mm1);
			}
			else
			{
				temp_vof[rightNb] += VOL2(mx, mz, mm2, mm1);
			}
		}

		return;
	}
}

void calcWorksXFull(int iCell, Real s1, Real s2)
{
	int leftNb, rightNb, temp;
	Real V1, V3, mm1, mm2;
	int ltop, lbot, rtop, rbot;

	// 2 3 < ltop  rtop > 2 3
	// 0 1 < lbot  rbot > 0 1

	// the child position if the destination
	ltop = 3;
	lbot = 1;
	rtop = 2;
	rbot = 0;

	V1 = DMAX(-s1, 0.0);
	temp_vof[iCell] += 1.0 - DMAX(s1, 0.0) + DMIN(s2, 0.0);
	V3 = DMAX(s2, 0.0);

	leftNb = cellNb[0][iCell];
	rightNb = cellNb[1][iCell];
	if (leftNb == 0 || rightNb == 0)
	{
		printf("***************************************\n");
		printf("Error in calcWorksX: invalid neighbours\n\n");
		exit(1);
	}

	// left
	if ((octLv[leftNb / 4] == octLv[iCell / 4]) && (cellChOct[leftNb] != 0))
	{
		// printf("\nworks\n");
		// 2 instead of 4 because V3 = 0.5 should fill up both
		temp_vof[4 * cellChOct[leftNb] + ltop] += 2 * V1; // 0 1 <
		temp_vof[4 * cellChOct[leftNb] + lbot] += 2 * V1; // 2 3 <
	}
	else
	{
		if (octLv[leftNb / 4] < octLv[iCell / 4])
		{
			temp_vof[leftNb] += 0.25 * V1;
		}
		else
		{
			temp_vof[leftNb] += V1;
		}
	}

	// right
	if ((octLv[rightNb / 4] == octLv[iCell / 4]) && (cellChOct[rightNb] != 0))
	{
		// printf("\nworks\n");
		// 2 instead of 4 because V3 = 0.5 should fill up both
		temp_vof[4 * cellChOct[rightNb] + rtop] += 2 * V3; // > 0 1
		temp_vof[4 * cellChOct[rightNb] + rbot] += 2 * V3; // > 2 3
	}
	else
	{
		if (octLv[rightNb / 4] < octLv[iCell / 4])
		{
			temp_vof[rightNb] += 0.25 * V3;
		}
		else
		{
			temp_vof[rightNb] += V3;
		}
	}
}

void octTreeXSwp_6x6(int iCell)
{
	int i, j, inv;
	Real mx, mz, alpha, s1, s2, mm1, mm2, V1, V2;
	Real cc[6][6];

	getCellNgbVOF_6x6(iCell, cc);
	i = 1;
	j = 1;

	/*
	   connect to your code
		  s1 = u[i][j]*tau/(deltax*deltay);
		  s2 = u[i+1][j]*tau/(deltax*deltay);
	*/
	s1 = 0.1;
	s2 = 0.1;

	if (cc[i][j] == 0.0)
	{
		work1[iCell] = 0.0;
		work2[iCell] = 0.0;
		work3[iCell] = 0.0;
	}
	else if (cc[i][j] == 1.0)
	{
		work1[iCell] = DMAX(-s1, 0.0);
		work2[iCell] = 1.0 - DMAX(s1, 0.0) + DMIN(s2, 0.0);
		work3[iCell] = DMAX(s2, 0.0);
	}
	else
	{

		/* normal to the interface */
		mm1 = cc[i - 1][j - 1] + 2. * cc[i - 1][j] + cc[i - 1][j + 1];
		mm2 = cc[i + 1][j - 1] + 2. * cc[i + 1][j] + cc[i + 1][j + 1];
		mx = mm1 - mm2;
		mm1 = cc[i - 1][j - 1] + 2. * cc[i][j - 1] + cc[i + 1][j - 1];
		mm2 = cc[i - 1][j + 1] + 2. * cc[i][j + 1] + cc[i + 1][j + 1];
		mz = mm1 - mm2;

		inv = 0;
		if (mx < 0.)
		{
			mm1 = -s1;
			s1 = -s2;
			s2 = mm1;
			mx = -mx;
			inv = 1;
		}
		mx = mx + 1.e-50;
		mz = fabs(mz) + 1.e-50;
		mm2 = DMAX(mx, mz);
		mx = mx / mm2;
		mz = mz / mm2;

		/* get alpha to determine the equation of the interface */
		mm1 = DMIN(mx, mz);
		/* compute the two critical volume fraction */
		V1 = 0.5 * mm1;
		V2 = 1.0 - V1;
		if (cc[i][j] <= V1)
			alpha = sqrt(2.0 * cc[i][j] * mm1);
		else if (cc[i][j] <= V2)
			alpha = cc[i][j] + 0.5 * mm1;
		else
			alpha = mm1 + 1.0 - sqrt(2.0 * (1.0 - cc[i][j]) * mm1);

		/* the new equation of the interface after advection */
		mx = mx / (1.0 - s1 + s2);
		alpha = alpha + mx * s1;

		/* calculate works */
		mm1 = DMAX(-s1, 0.0);
		mm2 = alpha + mx * mm1;
		work1[iCell] = VOL2(mx, mz, mm2, mm1);

		mm1 = 1.0 - DMAX(s1, 0.0) + DMIN(s2, 0.0);
		mm2 = alpha - mx * DMAX(s1, 0.0);
		work2[iCell] = VOL2(mx, mz, mm2, mm1);

		mm1 = DMAX(s2, 0.0);
		mm2 = alpha - mx;
		work3[iCell] = VOL2(mx, mz, mm2, mm1);

		/* symmetry conditions */
		if (inv == 1)
		{
			mm1 = work1[iCell];
			work1[iCell] = work3[iCell];
			work3[iCell] = mm1;
		}
	}

	return;
}
