
#ifndef WENO_H_INCLUDED
#define WENO_H_INCLUDED


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "tools.h"



#define DEFAULT_WENO_ORDER 3

//	int N = 512;
//	double H = 8.0/N;

//	extern int N;
//	extern double H;


//	void WENO_DX(double **uvalue, const int i, const int j, int N, double H, double *ux_p, double *ux_n, const int order);
//	void WENO_DY(double **uvalue, const int i, const int j, int N, double H, double *ux_p, double *ux_n, const int order);


void WENO_DX(double ***uvalue, const int i, const int j, const int k, const int N, const double H, double *ux_p, double *ux_n, const int order);
void WENO_DY(double ***uvalue, const int i, const int j, const int k, const int N, const double H, double *ux_p, double *ux_n, const int order);
void WENO_DZ(double ***uvalue, const int i, const int j, const int k, const int N, const double H, double *ux_p, double *ux_n, const int order);


void weno(double *dfp, double *dfn, double fi[], const double dx, const int order);



#endif
