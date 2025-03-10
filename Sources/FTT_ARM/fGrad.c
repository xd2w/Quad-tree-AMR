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
	Real cc[3][3];

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

	fprintf(fpGd, "%g %g ", x + .5 * dx, y + .5 * dy);
	fprintf(fpGd, "%g %g\n", mx, mz);
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

	fprintf(fpGd, "%g %g ", x + .5 * dx, y + .5 * dy);
	fprintf(fpGd, "%g %g\n", mx, mz);
	fprintf(fpGd, "\n");

	return;
}

void plotCellGrad_4x4(int iCell, FILE *fpGd)
{
	int i, j, invx, invz, iLv;
	Real mx, mz, alpha, s1, s2, mm1, mm2, V1, V2;
	Real cc[3][3];

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

	fprintf(fpGd, "%g %g ", x + .5 * dx, y + .5 * dy);
	fprintf(fpGd, "%g %g\n", mx, mz);
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

	fprintf(fpGd, "%g %g ", x + .5 * dx, y + .5 * dy);
	fprintf(fpGd, "%g %g\n", mx, mz);
	fprintf(fpGd, "\n");

	return;
}
