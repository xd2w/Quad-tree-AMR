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
	int i, j, invx, invz, scale1, scale2;
	int iLv, leftNb, rightNb;
	Real mx, mz, alpha, s1, s2, mm1, mm2, V1, V2;
	Real s11, s12, s21, s22, temp;
	Real cc[3][3];

	// Nighbours
	leftNb = cellNb[0][iCell];
	rightNb = cellNb[1][iCell];
	if (leftNb == 0 || rightNb == 0)
	{
		printf("***************************************\n");
		printf("Error in octTreeXSwap: invalid neighbours\n\n");
		exit(1);
	}

	// 1 for normal, 2 for twice as big, -1 for smaller
	scale1 = 1;
	scale2 = 1;

	Real ux1, ux2, x, y, dx, dy;
	iLv = octLv[iCell / 4];

	x = xCell[iCell];
	y = yCell[iCell];
	dx = dxCell[iLv];
	dy = dyCell[iLv];

	// #TODO : change later
	ux1 = (computeVX(x, y) + computeVX(x, y + dy)) * 0.5;
	ux2 = (computeVX(x + dx, y) + computeVX(x + dx, y + dy)) * 0.5;

	// default s1 s2
	s1 = (ux1 / dxCell[iLv]) * global_dt;
	s2 = (ux2 / dxCell[iLv]) * global_dt;

	// if neighbour is smaller
	if (octLv[leftNb / 4] == octLv[iCell / 4] && cellChOct[leftNb] != 0)
	{ // left
		scale1 = -1;
		ux1 = computeVX(x, y + 0.5 * dy);
		s1 = (ux1 / dxCell[iLv]) * global_dt; // more acurate s1

		ux1 = (computeVX(x, y) + computeVX(x, y + 0.5 * dy)) * 0.5;
		s11 = (ux1 / dxCell[iLv]) * global_dt; // uncorrected in scale
		ux1 = (computeVX(x, y + 0.5 * dy) + computeVX(x, y + dy)) * 0.5;
		s12 = (ux1 / dxCell[iLv]) * global_dt; // uncorrected in scale
	}
	if (octLv[rightNb / 4] == octLv[iCell / 4] && cellChOct[rightNb] != 0)
	{ // right
		scale1 = -1;
		ux2 = computeVX(x + dx, y + 0.5 * dy);
		s2 = (ux2 / dxCell[iLv]) * global_dt; // more acurate s2

		ux2 = (computeVX(x + dx, y) + computeVX(x + dx, y + 0.5 * dy)) * 0.5;
		s21 = (ux2 / dxCell[iLv]) * global_dt; // uncorrected in scale
		ux2 = (computeVX(x + dx, y + 0.5 * dy) + computeVX(x + dx, y + dy)) * 0.5;
		s22 = (ux2 / dxCell[iLv]) * global_dt; // uncorrected in scale
	}

	// if twice as big
	if (octLv[leftNb / 4] < octLv[iCell / 4])
	{ // left
		scale1 = 2;
	}
	if (octLv[rightNb / 4] < octLv[iCell / 4])
	{ // rigth
		scale2 = 2;
	}

	// colour func
	getCellNgbVOF(iCell, cc);

	if (cc[i][j] == 0.0)
	{
		return;
	}
	else if (cc[i][j] == 1.0)
	{
		calcWorksXFull(iCell, s1, s2, scale1, scale2, s11, s12, s21, s22);
	}
	else
	{
		i = 1;
		j = 1;
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

			mm1 = -s11;
			s11 = -s21;
			s21 = mm1;
			mm1 = -s12;
			s12 = -s22;
			s22 = mm1;

			mx = -mx;
			invx = 1;
		}

		invz = 0;
		if (mz < 0.0)
		{

			invz = 1;
			// SWAP(s11, s12)
			// SWAP(s21, s22)
			mm1 = s11;
			s11 = s12;
			s12 = mm1;
			mm1 = s21;
			s21 = s22;
			s22 = mm1;
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

		calcWorksX(iCell, vof[iCell], alpha, mx, mz, invx, invz, s1, s2, scale1, scale2, s11, s12, s21, s22);
	}
}

void calcWorksX(int iCell, Real vofVal, Real alpha, Real mx, Real mz, int invx, int invz, Real s1, Real s2, int scale1, int scale2, Real s11, Real s12, Real s21, Real s22)
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

	// if (vofVal < 0)
	// {
	// 	return;
	// }

	mm1 = 1.0 - DMAX(s1, 0.0) + DMIN(s2, 0.0);
	mm2 = alpha - mx * DMAX(s1, 0.0);
	temp_vof[iCell] += VOL2(mx, mz, mm2, mm1);

	if (scale1 < 0)
	{ // if neighbour is smaller

		// top
		mm1 = DMAX(-s11, 0.0);
		mm2 = 2 * (alpha + mx * mm1 - mz * 0.5);
		temp_vof[4 * cellChOct[leftNb] + ltop] += VOL2(mx, mz, mm2, 2 * mm1);

		// bottom
		mm1 = DMAX(-s12, 0.0);
		mm2 = 2 * (alpha + mx * mm1);
		temp_vof[4 * cellChOct[leftNb] + lbot] += VOL2(mx, mz, mm2, 2 * mm1);
	}
	else
	{
		mm1 = DMAX(-s1, 0.0);
		mm2 = alpha + mx * mm1;
		temp_vof[leftNb] += VOL2(mx, mz, mm2, mm1) / (scale1 * scale1);
	}

	if (scale2 < 0)
	{ // if neighbour is smaller

		// top
		mm1 = DMAX(s21, 0.0);
		mm2 = 2 * (alpha - mx - mz * 0.5);
		temp_vof[4 * cellChOct[rightNb] + rtop] += VOL2(mx, mz, mm2, 2 * mm1);

		// bottom
		mm1 = DMAX(s22, 0.0);
		mm2 = 2 * (alpha - mx);
		temp_vof[4 * cellChOct[rightNb] + rbot] += VOL2(mx, mz, mm2, 2 * mm1);
	}
	else
	{
		mm1 = DMAX(s2, 0.0);
		mm2 = (alpha - mx);
		temp_vof[rightNb] += VOL2(mx, mz, mm2, mm1) / (scale2 * scale2);
	}

	return;
}

void calcWorksXFull(int iCell, Real s1, Real s2, int scale1, int scale2, Real s11, Real s12, Real s21, Real s22)
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
	if (scale1 < 0)
	{
		// 2 instead of 4 because V3 = 0.5 should fill up both
		temp_vof[4 * cellChOct[leftNb] + ltop] += 2 * DMAX(-s11, 0.0); // 0 1 <
		temp_vof[4 * cellChOct[leftNb] + lbot] += 2 * DMAX(-s12, 0.0); // 2 3 <
	}
	else
	{
		temp_vof[leftNb] += V1 / (scale1 * scale1);
	}

	// right
	if (scale2 < 0)
	{
		// printf("\nworks\n");
		// 2 instead of 4 because V3 = 0.5 should fill up both
		temp_vof[4 * cellChOct[rightNb] + rtop] += 2 * DMAX(s21, 0.0); // > 0 1
		temp_vof[4 * cellChOct[rightNb] + rbot] += 2 * DMAX(s22, 0.0); // > 2 3
	}
	else
	{
		temp_vof[rightNb] += V3 / (scale2 * scale2);
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
