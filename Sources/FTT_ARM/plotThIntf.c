#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "ftt.h"
#include "nrutil.h"
#include "pfplib.h"

void plotTheoreticalInterf(int ndata)
{
  int i;
  char fname[] = "DATA/thintf.000";
  FILE *fp;

  i = ndata;
  int i1 = i % 10;
  i /= 10;
  int i2 = i % 10;
  i /= 10;
  int i3 = i % 10;
  fname[12] = '0' + i3;
  fname[13] = '0' + i2;
  fname[14] = '0' + i1;

  fp = fopen(fname, "w");
  if (fp == NULL)
  {
    printf("Error opening file %s\n", fname);
    exit(1);
  }

  // fprintf(fp, "%d\n", nThPoints + 1);

  for (int i = 0; i < nThPoints; i++)
  {
    fprintf(fp, "%g %g\n", xThPoints[i], yThPoints[i]);
  }
  fprintf(fp, "%g %g\n", xThPoints[0], yThPoints[0]);
  fprintf(fp, "\n");
}

void advTheoreticalInterf(void)
{
  Real x, y;

  for (int i = 0; i < nThPoints; i++)
  {
    x = xThPoints[i];
    y = yThPoints[i];

    xThPoints[i] = x + computeVX(x, y) * global_dt;
    yThPoints[i] = y + computeVY(x, y) * global_dt;
  }
}

void initTheoreticalInterf(void)
{
  int i, j, nCircle;
  Real x, y, dx, dy;
  Real xc, yc, radius;
  Real delta, a, b, A, B, C;
  Real theta, phi;

  nThPoints = 1000;

  xc = 0.5;
  yc = 0.75;
  radius = 0.2;
  // dfetch("xc", &xc);
  // dfetch("yc", &yc);
  // dfetch("radius", &radius);

  A = init_VOF_coefs[0];
  B = init_VOF_coefs[1];
  C = init_VOF_coefs[2];

  phi = 0.5 * atan2(B, A - C);
  if (phi > 1e50)
  {
    phi = 0;
  }

  a = radius / (sqrt(A * cos(phi) * cos(phi) + B * cos(phi) * sin(phi) + C * sin(phi) * sin(phi)) + 1e-50);
  b = radius / (sqrt(A * sin(phi) * sin(phi) - B * cos(phi) * sin(phi) + C * cos(phi) * cos(phi)) + 1e-50);

  delta = 2 * M_PI / nThPoints;

  for (i = 0; i < nThPoints; i++)
  {
    theta = (i + 0.5) * delta;
    xThPoints[i] = xc + a * cos(theta) * cos(phi) - b * sin(theta) * sin(phi);
    yThPoints[i] = yc + a * cos(theta) * sin(phi) + b * sin(theta) * cos(phi);

    // printf("%f %f\n", xThPoints[i], yThPoints[i]);
  }
  // exit(1);
}
