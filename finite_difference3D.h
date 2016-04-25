
#ifndef FINITE_DIFFERENCE_H_INCLUDED
#define FINITE_DIFFERENCE_H_INCLUDED

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "tools.h"
#include "constants.h"

double Godunov(double signage, double a, double b, double c, double d, double e, double f);
int FD_DX(double ***phi, int i, int j, int k, const int N, double H, double *Dxp, double *Dxn, int mode);
int FD_DY(double ***phi, int i, int j, int k, const int N, double H, double *Dyp, double *Dyn, int mode);
int FD_DZ(double ***phi, int i, int j, int k, const int N, double H, double *Dzp, double *Dzn, int mode);
double XDER_C(double ***level, int i, int j, int k, int N, double H);
double YDER_C(double ***level, int i, int j, int k, int N, double H);
double ZDER_C(double ***level, int i, int j, int k, int N, double H);
double XDER2(double ***level, int i, int j, int k, int N, double H);
double YDER2(double ***level, int i, int j, int k, int N, double H);
double ZDER2(double ***level, int i, int j, int k, int N, double H);
double secondpartialxy(double ***level, int i, int j, int k, int N, double H);
double secondpartialxz(double ***level, int i, int j, int k, int N, double H);
double secondpartialyz(double ***level, int i, int j, int k, int N, double H);



#endif
