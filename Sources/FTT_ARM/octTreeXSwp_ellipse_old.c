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

// s1(i=1) and s2(i=2)
// Real getNbSi(int iCell, int nbCell, int i, int invx)
// {
// 	Real ux1, ux2, x, y, dx, dy, s1, s2;
// 	int iLv, nbLv;
//
// 	iLv = octLv[iCell / 4];
// 	nbLv = octLv[nbCell / 4];
//
// 	x = xCell[iCell];
// 	y = yCell[iCell];
// 	dx = dxCell[iLv];
// 	dy = dyCell[iLv];
//
// 	// # TODO: do this with u[iCell] to be flexible
// 	if (i == 1)
// 	{
// 		ux1 = (computeVX(x, y) + computeVX(x, y + dy)) * 0.5;
// 		s1 = (ux1 / dxCell[nbLv]) * global_dt;
// 		return ux1;
// 	}
//
// 	if (i == 2)
// 	{
// 		ux2 = (computeVX(x + dx, y) + computeVX(x + dx, y + dy)) * 0.5;
// 		s2 = (ux2 / dxCell[nbLv]) * global_dt;
// 		return ux2;
// 	}
//
// 	printf("\n\nError in octTreeXSwp: getS_i wrong i\n\n");
// 	exit(1);
// 	return 0;
// }

void calcs1s2Child(int iCell, int iLocal, Real cs1, Real cs2, Real calpha, Real cmx, int invx)
{

	// calpha = alpha - mx * s1;
	// cmx = mx * (1.0 - s1 + s2);

	int iLv;
	Real ux1, ux2, x, y, dx, dy;

	float px[] = {0, 1, 0, 1};
	float py[] = {0, 0, 1, 1};

	x = xCell[iCell];
	y = yCell[iCell];

	iLv = octLv[iCell / 4];
	dx = 0.5 * dxCell[iLv];
	dy = 0.5 * dyCell[iLv];

	x += dx * px[iLocal];
	y += dy * py[iLocal];

	// TODO: proof of concept change later
	ux1 = (computeVX(x, y) + computeVX(x, y + dy)) * 0.5;
	ux2 = (computeVX(x + dx, y) + computeVX(x + dx, y + dy)) * 0.5;

	cs1 = (ux1 / dx) * global_dt;
	cs2 = (ux2 / dx) * global_dt;

	if (invx)
	{
		ux1 = -cs1;
		cs1 = -cs2;
		cs2 = ux1;
	}

	cmx = cmx / (1.0 - cs1 + cs2);
	calpha = calpha + cmx * cs1;
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

	Real ux1, ux2, x, y, dx, dy;

	iLv = octLv[iCell / 4];

	x = xCell[iCell];
	y = yCell[iCell];
	dx = dxCell[iLv];
	dy = dyCell[iLv];

	ux1 = (computeVX(x, y) + computeVX(x, y + dy)) * 0.5;
	ux2 = (computeVX(x + dx, y) + computeVX(x + dx, y + dy)) * 0.5;

	s1 = (ux1 / dxCell[iLv]) * global_dt;
	s2 = (ux2 / dxCell[iLv]) * global_dt;

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
	Real V1, V3, mm1, mm2, ux1, ux2, x, y, dx, dy, cs1, cs2, calpha, cmx;
	int ltop, lbot, rtop, rbot, iLv, nbCell;

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

			// get neighbour's s2 unless x inverted
			temp = invx ? 1 : 2;
			calpha = alpha;
			cmx = mx;
			cs1 = s1;
			cs2 = s2;

			// getting original alpha and mx
			// calpha = alpha - mx * s1;
			// cmx = mx * (1.0 - s1 + s2);
			// calcs1s2Child(iCell, rtop, cs1, cs2, calpha, cmx, invx);

			// top
			mm1 = DMAX(-cs1, 0.0);
			mm2 = 2 * (calpha + cmx * mm1 - mz * 0.5);
			temp_vof[4 * cellChOct[leftNb] + ltop] += VOL2(cmx, mz, mm2, 2 * mm1);

			// calpha = alpha - mx * s1;
			// cmx = mx * (1.0 - s1 + s2);
			// calcs1s2Child(iCell, rbot, cs1, cs2, calpha, cmx, invx);

			// bottom
			mm1 = DMAX(-cs1, 0.0);
			mm2 = 2 * (calpha + cmx * mm1);
			temp_vof[4 * cellChOct[leftNb] + lbot] += VOL2(cmx, mz, mm2, 2 * mm1);
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

			// get neighbour's s1 unless x inverted
			temp = invx ? 2 : 1;
			calpha = alpha;
			cmx = mx;
			cs1 = s1;
			cs2 = s2;

			// top
			// calpha = alpha - mx * s1;
			// cmx = mx * (1.0 - s1 + s2);
			// calcs1s2Child(iCell, ltop, cs1, cs2, calpha, cmx, invx);

			mm1 = DMAX(cs2, 0.0);
			mm2 = 2 * (calpha - cmx - mz * 0.5);
			temp_vof[4 * cellChOct[rightNb] + rtop] += VOL2(cmx, mz, mm2, 2 * mm1);

			// bottom
			// calpha = alpha - mx * s1;
			// cmx = mx * (1.0 - s1 + s2);
			// calcs1s2Child(iCell, lbot, cs1, cs2, calpha, cmx, invx);

			mm1 = DMAX(cs2, 0.0);
			mm2 = 2 * (calpha - cmx);
			temp_vof[4 * cellChOct[rightNb] + rbot] += VOL2(cmx, mz, mm2, 2 * mm1);
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

	Real calpha, cmx, cs1, cs2;

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
