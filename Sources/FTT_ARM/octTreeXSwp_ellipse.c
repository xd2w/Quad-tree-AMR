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

	int leftNb = cellNb[0][iCell];
	int rightNb = cellNb[1][iCell];
	if (leftNb == 0 || rightNb == 0)
	{
		printf("***************************************\n");
		printf("Error in calcWorksX: invalid neighbours\n\n");
		exit(1);
	}

	Real ux11, ux12, ux21, ux22;
	Real s11, s12, s21, s22;

	ux11 = (computeVX(x, y) + computeVX(x, y + 0.5 * dy)) * 0.5;
	ux12 = (computeVX(x, y + 0.5 * dy) + computeVX(x, y + dy)) * 0.5;

	ux21 = (computeVX(x + dx, y) + computeVX(x + dx, y + 0.5 * dy)) * 0.5;
	ux22 = (computeVX(x + dx, y + 0.5 * dy) + computeVX(x + dx, y + dy)) * 0.5;

	s11 = (ux11 / dxCell[iLv]) * global_dt;
	s12 = (ux12 / dxCell[iLv]) * global_dt;
	s21 = (ux21 / dxCell[iLv]) * global_dt;
	s22 = (ux22 / dxCell[iLv]) * global_dt;

	if (octLv[leftNb / 4] == octLv[iCell / 4] && cellChOct[leftNb] != 0)
	{
		s1 = 0.5 * (s11 + s12);
	}

	if (octLv[rightNb / 4] == octLv[iCell / 4] && cellChOct[rightNb] != 0)
	{
		s2 = 0.5 * (s21 + s22);
	}

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
	Real V1, V2, V3, mm1, mm2, ux1, ux2, x, y, dx, dy, cs1, cs2, calpha, cmx, vol;
	int ltop, lbot, rtop, rbot, iLv, nbCell;

	// -----------
	// s12	s22
	// s11	s21

	Real s11, s12, s21, s22;
	Real ux11, ux12, ux21, ux22;
	// Real x, y, dx, dy;

	iLv = octLv[iCell / 4];

	x = xCell[iCell];
	y = yCell[iCell];
	dx = dxCell[iLv];
	dy = dyCell[iLv];

	ux11 = (computeVX(x, y) + computeVX(x, y + 0.5 * dy)) * 0.5;
	ux12 = (computeVX(x, y + 0.5 * dy) + computeVX(x, y + dy)) * 0.5;

	ux21 = (computeVX(x + dx, y) + computeVX(x + dx, y + 0.5 * dy)) * 0.5;
	ux22 = (computeVX(x + dx, y + 0.5 * dy) + computeVX(x + dx, y + dy)) * 0.5;

	s11 = (ux11 / dxCell[iLv]) * global_dt;
	s12 = (ux12 / dxCell[iLv]) * global_dt;
	s21 = (ux21 / dxCell[iLv]) * global_dt;
	s22 = (ux22 / dxCell[iLv]) * global_dt;

	if (invx)
	{
		mm1 = -s11;
		s11 = -s21;
		s21 = mm1;
		mm1 = -s12;
		s12 = -s22;
		s22 = mm1;
	}

	if (invz)
	{
		mm1 = s11;
		s11 = s12;
		s12 = mm1;
		mm1 = s21;
		s21 = s22;
		s22 = mm1;
	}

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
		V2 = vofVal;

		if (octLv[leftNb / 4] == octLv[iCell / 4] && cellChOct[leftNb] != 0)
		{ // if neighbour is smaller

			// top
			mm1 = DMAX(-s12, 0.0);
			mm2 = 2 * (alpha + mx * mm1 - mz * 0.5);
			vol = VOL2(mx, mz, mm2, 2 * mm1);
			temp_vof[4 * cellChOct[leftNb] + ltop] += vol;
			V2 -= 0.25 * vol;

			// bottom
			mm1 = DMAX(-s11, 0.0);
			mm2 = 2 * (alpha + mx * mm1);
			vol = VOL2(mx, mz, mm2, 2 * mm1);
			temp_vof[4 * cellChOct[leftNb] + lbot] += vol;
			V2 -= 0.25 * vol;
		}
		else
		{
			mm1 = DMAX(-s1, 0.0);
			mm2 = alpha + mx * mm1;
			if (octLv[leftNb / 4] < octLv[iCell / 4])
			{
				// splitCell_smart(leftNb);
				// calcWorksX(iCell, vof[iCell], alpha, mx, mz, invx, invz, s1, s2);
				// return;
				vol = VOL2(mx, mz, mm2, mm1);
				temp_vof[leftNb] += 0.25 * vol;
				V2 -= vol;
			}
			else
			{
				vol = VOL2(mx, mz, mm2, mm1);
				temp_vof[leftNb] += vol;
				V2 -= vol;
			}
		}

		if (octLv[rightNb / 4] == octLv[iCell / 4] && cellChOct[rightNb] != 0)
		{ // if neighbour is smaller

			// top
			mm1 = DMAX(s22, 0.0);
			mm2 = 2 * (alpha - mx - mz * 0.5);
			vol = VOL2(mx, mz, mm2, 2 * mm1);
			temp_vof[4 * cellChOct[rightNb] + rtop] += vol;
			V2 -= 0.25 * vol;

			// bottom
			mm1 = DMAX(s21, 0.0);
			mm2 = 2 * (alpha - mx);
			vol = VOL2(mx, mz, mm2, 2 * mm1);
			temp_vof[4 * cellChOct[rightNb] + rbot] += vol;
			V2 -= 0.25 * vol;
		}
		else
		{
			mm1 = DMAX(s2, 0.0);
			mm2 = (alpha - mx);
			if (octLv[rightNb / 4] < octLv[iCell / 4])
			{
				// splitCell_smart(rightNb);
				// calcWorksX(iCell, vof[iCell], alpha, mx, mz, invx, invz, s1, s2);
				// return;
				vol = VOL2(mx, mz, mm2, mm1);
				temp_vof[rightNb] += 0.25 * vol;
				V2 -= vol;
			}
			else
			{
				vol = VOL2(mx, mz, mm2, mm1);
				temp_vof[rightNb] += vol;
				V2 -= vol;
			}
		}

		// V2 = DMAX(V2, 0.0);
		// mm1 = 1.0 - DMAX(s1, 0.0) + DMIN(s2, 0.0);
		// mm2 = alpha - mx * DMAX(s1, 0.0);
		// if (fabs(VOL2(mx, mz, mm2, mm1) - V2) > 1e-10)
		// {
		// 	printf("\nlarge variation in streaming VOF\n");
		// 	// exit(1);
		// }

		// temp_vof[iCell] += V2;

		return;
	}
}

void calcWorksXFull(int iCell, Real s1, Real s2)
{
	int leftNb, rightNb, temp;
	Real V1, V2, V3, mm1, mm2;
	int ltop, lbot, rtop, rbot;

	// -----------

	Real s11, s12, s21, s22, cs1, cs2;
	Real ux11, ux12, ux21, ux22;
	Real x, y, dx, dy;

	int iLv = octLv[iCell / 4];

	x = xCell[iCell];
	y = yCell[iCell];
	dx = dxCell[iLv];
	dy = dyCell[iLv];

	ux11 = (computeVX(x, y) + computeVX(x, y + 0.5 * dy)) * 0.5;
	ux12 = (computeVX(x, y + 0.5 * dy) + computeVX(x, y + dy)) * 0.5;

	ux21 = (computeVX(x + dx, y) + computeVX(x + dx, y + 0.5 * dy)) * 0.5;
	ux22 = (computeVX(x + dx, y + 0.5 * dy) + computeVX(x + dx, y + dy)) * 0.5;

	s11 = (ux11 / dxCell[iLv]) * global_dt;
	s12 = (ux12 / dxCell[iLv]) * global_dt;
	s21 = (ux21 / dxCell[iLv]) * global_dt;
	s22 = (ux22 / dxCell[iLv]) * global_dt;

	// 2 3 < ltop  rtop > 2 3
	// 0 1 < lbot  rbot > 0 1

	// the child position if the destination
	ltop = 3;
	lbot = 1;
	rtop = 2;
	rbot = 0;

	V1 = DMAX(-s1, 0.0);
	// temp_vof[iCell] += 1.0 - DMAX(s1, 0.0) + DMIN(s2, 0.0);
	// V2 = 1;
	cs1 = s1;
	cs2 = s2;
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
		temp_vof[4 * cellChOct[leftNb] + ltop] += 2 * DMAX(-s12, 0.0); // 0 1 <
		temp_vof[4 * cellChOct[leftNb] + lbot] += 2 * DMAX(-s11, 0.0); // 2 3 <
		// temp_vof[iCell] += 0.5 * ((fmax(s1, 0.0)) - (fmax(s11, 0.0)));
		// temp_vof[iCell] += 0.5 * ((fmax(s1, 0.0)) - (fmax(s12, 0.0)));
		cs1 = 0.5 * (s11 + s12);
	}
	else
	{
		if (octLv[leftNb / 4] < octLv[iCell / 4])
		{
			temp_vof[leftNb] += 0.25 * V1;
			V2 -= V1;
		}
		else
		{
			temp_vof[leftNb] += V1;
			V2 -= V1;
		}
	}

	// right
	if ((octLv[rightNb / 4] == octLv[iCell / 4]) && (cellChOct[rightNb] != 0))
	{
		// printf("\nworks\n");
		// 2 instead of 4 because V3 = 0.5 should fill up both
		temp_vof[4 * cellChOct[rightNb] + rtop] += 2 * DMAX(s22, 0.0); // > 0 1
		temp_vof[4 * cellChOct[rightNb] + rbot] += 2 * DMAX(s21, 0.0); // > 2 3
		// temp_vof[iCell] -= 0.5 * (fmin(s2, 0.0) - fmin(s21, 0.0));
		// temp_vof[iCell] -= 0.5 * (fmin(s2, 0.0) - fmin(s22, 0.0));
		// V2 -= 0.5 * DMAX(s21, 0.0);
		// V2 -= 0.5 * DMAX(s22, 0.0);
		cs2 = 0.5 * (s21 + s22);
	}
	else
	{
		if (octLv[rightNb / 4] < octLv[iCell / 4])
		{
			temp_vof[rightNb] += 0.25 * V3;
			V2 -= V3;
		}
		else
		{
			temp_vof[rightNb] += V3;
			V2 -= V3;
		}
	}
	V2 = DMAX(V2, 0.0);
	// if (fabs(1.0 - DMAX(s1, 0.0) + DMIN(s2, 0.0) - V2) > 1e-10)
	// {
	// 	printf("\nlarge variation in streaming VOF FULL\n");
	// 	// exit(1);
	// }
	// temp_vof[iCell] += V2;
	temp_vof[iCell] += 1.0 - DMAX(cs1, 0.0) + DMIN(cs2, 0.0);
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
