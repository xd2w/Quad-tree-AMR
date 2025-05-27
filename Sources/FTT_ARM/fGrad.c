#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "ftt.h"
#include "nrutil.h"

#include "pfplib.h"

void plotCellGrad(int iCell, FILE *fpGd)
{
	int i, j, invx, invz, iLv;
	Real mx, mz, alpha, s1, s2, mm1, mm2, V1, V2;
	Real cc[3][3], cc6[6][6], cc4[4][4];

	getCellNgbVOF(iCell, cc);
	i = 1;
	j = 1;

	Real x, y, dx, dy;

	iLv = octLv[iCell / 4];

	x = xCell[iCell];
	y = yCell[iCell];
	dx = dxCell[iLv];
	dy = dyCell[iLv];

	/* normal to the interface */
	mm1 = cc[i - 1][j - 1] + 2. * cc[i - 1][j] + cc[i - 1][j + 1];
	mm2 = cc[i + 1][j - 1] + 2. * cc[i + 1][j] + cc[i + 1][j + 1];
	mx = 0.5 * (mm1 - mm2) / dx;
	mm1 = cc[i - 1][j - 1] + 2. * cc[i][j - 1] + cc[i + 1][j - 1];
	mm2 = cc[i - 1][j + 1] + 2. * cc[i][j + 1] + cc[i + 1][j + 1];
	mz = 0.5 * (mm1 - mm2) / dy;

	getCellNgbVOF_6x6(iCell / 4, cc6);
	smoothUnifrom(cc6, cc4);

	int px[] = {1, 2, 1, 2};
	int pz[] = {1, 1, 2, 2};

	i = px[iCell % 4];
	j = pz[iCell % 4];

	if (0 < cc4[i][j] && cc4[i][j] < 1)
	{
		if (kappaMode == 1)
		{
			mm2 = kappaHF(iCell, cc);
		}
		else if (kappaMode == 2)
		{
			mm2 = kappaMeier(iCell);
		}
		else
		{
			mm2 = kappaBarickALELike(iCell, cc);
		}
		// mm2 = kappaBarickALELike(iCell, cc6);
		// mm2 = kappaHF(iCell, cc6);
		// mm2 = kappaBarickALELike_wider(iCell, cc6);
	}
	else
	{
		mm2 = 0;
	}

	fprintf(fpGd, "%g %g ", x + .5 * dx, y + .5 * dy);
	fprintf(fpGd, "%g %g %g\n", mx, mz, mm2);
	fprintf(fpGd, "\n");

	return;
}

void plotCellGradSmoothed(int iCell, FILE *fpGd)
{
	int i, j, invx, invz, iLv;
	Real mx, mz, alpha, s1, s2, mm1, mm2, V1, V2;
	Real ccs[6][6], cc[4][4];

	getCellNgbVOF_6x6(iCell / 4, ccs);
	smooth2ndOrd(ccs, cc);
	int px[] = {1, 1, 2, 2};
	int pz[] = {1, 2, 1, 2};
	i = px[iCell % 4];
	j = pz[iCell % 4];

	Real x, y, dx, dy;

	iLv = octLv[iCell / 4];

	x = xCell[iCell];
	y = yCell[iCell];
	dx = dxCell[iLv];
	dy = dyCell[iLv];

	/* normal to the interface */
	mm1 = cc[i - 1][j - 1] + 2. * cc[i - 1][j] + cc[i - 1][j + 1];
	mm2 = cc[i + 1][j - 1] + 2. * cc[i + 1][j] + cc[i + 1][j + 1];
	mx = 0.5 * (mm1 - mm2) / dx;
	mm1 = cc[i - 1][j - 1] + 2. * cc[i][j - 1] + cc[i + 1][j - 1];
	mm2 = cc[i - 1][j + 1] + 2. * cc[i][j + 1] + cc[i + 1][j + 1];
	mz = 0.5 * (mm1 - mm2) / dy;

	mm2 = kappaBarickALELike(iCell, ccs);

	fprintf(fpGd, "%g %g ", x + .5 * dx, y + .5 * dy);
	fprintf(fpGd, "%g %g %g\n", mx, mz, mm2);
	fprintf(fpGd, "\n");

	return;
}

void plotCellGrad_4x4(int iCell, FILE *fpGd)
{
	int i, j, invx, invz, iLv;
	Real mx, mz, alpha, s1, s2, mm1, mm2, V1, V2;
	Real cc[3][3], cc6[6][6];

	getCellNgbVOF(iCell, cc);
	i = 1;
	j = 1;

	Real x, y, dx, dy;

	iLv = octLv[iCell / 4];

	x = xCell[iCell];
	y = yCell[iCell];
	dx = dxCell[iLv];
	dy = dyCell[iLv];

	/* normal to the interface */
	mx = (cc[i + 1][j] + cc[i + 1][j + 1] - cc[i][j] - cc[i][j + 1]) / (2 * dx + 1e-50);
	mz = (cc[i][j + 1] + cc[i + 1][j + 1] - cc[i][j] - cc[i + 1][j]) / (2 * dy + 1e-50);

	getCellNgbVOF_6x6(iCell / 4, cc6);
	mm2 = kappaBarickALELike(iCell, cc6);

	fprintf(fpGd, "%g %g ", x + .5 * dx, y + .5 * dy);
	fprintf(fpGd, "%g %g %g\n", mx, mz, mm2);
	fprintf(fpGd, "\n");

	return;
}

void plotCellGradSmoothed_4x4(int iCell, FILE *fpGd)
{
	int i, j, invx, invz, iLv;
	Real mx, mz, alpha, s1, s2, mm1, mm2, V1, V2;
	Real ccs[6][6], cc[4][4];

	getCellNgbVOF_6x6(iCell / 4, ccs);
	smooth2ndOrd(ccs, cc);
	// smoothUnifrom(ccs, cc);
	int px[] = {1, 1, 2, 2};
	int pz[] = {1, 2, 1, 2};
	i = px[iCell % 4];
	j = pz[iCell % 4];

	Real x, y, dx, dy;

	iLv = octLv[iCell / 4];

	x = xCell[iCell];
	y = yCell[iCell];
	dx = dxCell[iLv];
	dy = dyCell[iLv];

	/* normal to the interface */
	mx = (ccs[i + 1][j] + ccs[i + 1][j + 1] - ccs[i][j] - ccs[i][j + 1]) / (2 * dx + 1e-50);
	mz = (ccs[i][j + 1] + ccs[i + 1][j + 1] - ccs[i][j] - ccs[i + 1][j]) / (2 * dy + 1e-50);

	mm2 = kappaBarickALELike(iCell, ccs);

	fprintf(fpGd, "%g %g ", x + .5 * dx, y + .5 * dy);
	fprintf(fpGd, "%g %g %g\n", mx, mz, mm2);
	fprintf(fpGd, "\n");

	return;
}

void plotVOFContourCell(FILE *fp, int iCell, int intrinsicMaxLev)
{
	Real dx, dy, mindx, mindy, xp, yp, x, y;

	xp = xCell[iCell];
	yp = yCell[iCell];
	dx = dxCell[octLv[iCell / 4]];
	dy = dyCell[octLv[iCell / 4]];
	mindx = dxCell[intrinsicMaxLev];
	mindy = dyCell[intrinsicMaxLev];

	y = yp + 0.5 * mindy;
	while (y < yp + dy)
	{
		x = xp + 0.5 * mindx;
		while (x < xp + dx)
		{
			fprintf(fp, "%f %f %f\n", x, y, vof[iCell]);
			x += mindx;
		}
		y += mindy;
	}
}

void plotVOFContour(int ndata)
{
	int intrinsicMaxLev = 0;
	int iOct, iCell;
	char fname[] = "DATA/vofcont.000";
	fname[13] = '0' + ndata / 100;
	fname[14] = '0' + (ndata / 10) % 10;
	fname[15] = '0' + ndata % 10;
	FILE *fp = fopen(fname, "w");

	for (iOct = 0; iOct < numberOfOcts; iOct++)
	{
		if (intrinsicMaxLev < octLv[iOct])
		{
			intrinsicMaxLev = octLv[iOct];
		}
	}

	for (iCell = 0; iCell < numberOfCells; iCell++)
	{
		if (cellChOct[iCell] == 0)
		{
			plotVOFContourCell(fp, iCell, intrinsicMaxLev);
		}
	}
	fclose(fp);
}

void plotCellGradAtIntf(int ndata)
{
	int iCell;
	char fname[] = "DATA/fgrd.000";
	fname[10] = '0' + ndata / 100;
	fname[11] = '0' + (ndata / 10) % 10;
	fname[12] = '0' + ndata % 10;
	FILE *fp = fopen(fname, "w");

	for (iCell = 0; iCell < numberOfCells; iCell++)
	{
		if (cellChOct[iCell] == 0)
		{
			// if (0 < vof[iCell] && vof[iCell] < 1)
			// {
			plotCellGrad(iCell, fp);
			// plotCellGrad_4x4(iCell, fp);
			// plotCellGradSmoothed(iCell, fp);
			// plotCellGradSmoothed_4x4(iCell, fp);
			// }
		}
	}
}
