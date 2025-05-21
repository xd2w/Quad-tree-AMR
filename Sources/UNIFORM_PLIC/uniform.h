#include "variable.h"

#ifndef UNIFORM_H
#define UNIFORM_H

extern Real Lx;
extern Real Ly;
extern Real dx;
extern Real dy;
extern int nx;
extern int ny;
extern Real dt;
extern Real t_total;

extern Real1D x;
extern Real1D y;
extern Real2D vof;
extern Real2D temp_vof;
extern Real2D ux;
extern Real2D uy;

// possion solver
extern Real2D phi;
extern Real2D ux_star;
extern Real2D uy_star;

// LBM

// extern Real w[9];
// extern Real e_x[9];
// extern Real e_y[9];
// extern Real2D f[9];
// extern Real2D f_2[9];

extern Real2D f0;
extern Real2D f1;
extern Real2D f2;
extern Real2D f3;
extern Real2D f4;
extern Real2D f5;
extern Real2D f6;
extern Real2D f7;
extern Real2D f8;

extern Real2D f0_2;
extern Real2D f1_2;
extern Real2D f2_2;
extern Real2D f3_2;
extern Real2D f4_2;
extern Real2D f5_2;
extern Real2D f6_2;
extern Real2D f7_2;
extern Real2D f8_2;

extern Real2D rho;
extern Real tau;

extern void init_uniform_mesh(void);
extern void plic(void);
extern void plot(int ndata);
extern void plot_vof(int ndata);
extern void solveV(void);
extern void init_solver_vars(void);
extern void set_BC(void);
extern void set_IC(void);

extern Real VOL2(Real mx, Real mz, Real alpha, Real b);
extern void copy_vof(Real2D from, Real2D to);
extern void clear_vof(Real2D from);
extern void correct_vof(Real2D from);
extern void plicx(Real2D from, Real2D to);
extern void plicy(Real2D from, Real2D to);
extern void plot_one(FILE *fp, int i, int j);

#endif