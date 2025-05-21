#include "variable.h"

Real Lx;
Real Ly;
Real dx;
Real dy;
int nx;
int ny;
Real1D x;
Real1D y;
Real dt;
Real t_total;

Real2D vof;
Real2D temp_vof;
Real2D ux;
Real2D uy;

// possion solver
Real2D phi;
Real2D ux_star;
Real2D uy_star;

// LBM

// static double w[9] = {4.0 / 9.0, 1.0 / 9.0, 1.0 / 9.0, 1.0 / 9.0, 1.0 / 9.0, 1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0};
// static double e_x[9] = {0, 0, 1, 0, -1, -1, 1, 1, -1};
// static double e_y[9] = {0, 1, 0, -1, 0, 1, 1, -1, -1};

Real2D f0;
Real2D f1;
Real2D f2;
Real2D f3;
Real2D f4;
Real2D f5;
Real2D f6;
Real2D f7;
Real2D f8;

Real2D f0_2;
Real2D f1_2;
Real2D f2_2;
Real2D f3_2;
Real2D f4_2;
Real2D f5_2;
Real2D f6_2;
Real2D f7_2;
Real2D f8_2;

Real2D rho;
Real tau;