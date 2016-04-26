
#ifndef FUNCTIONS_H_DEFINED
#define FUNCTIoNS_H_DEFINED


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "PDE.h"
#include "frame.h"
#include "constants.h"
#include "tools.h"

			//	BOTH OF THE ABOVE DEPEND ON THE WAVE NUMBER	//
#define PSI_ETA 0.012	//	determine how smooth the point source is, the smaller the more concentrated 	//
#define DIST_MODE 0 	//	0 for circle(s), 1 for level set function	//
#define NONHOM_MODE 0	//	1 for nonhomogeneous mode, 0 for homogeneous mode.	//
#define DELTA_MODE "cosine"
#define J_M 2		//	the mode for J(z)	//
#define ONE_SPHERE 1
#define MOMENT_ORDER 1

#define CONSTA1 -759.2781934172483 
#define CONSTB1 446.2604260472818 
 
#define CONSTA2 14317.969703708994 
#define CONSTB2 -16509.044867497203 
#define CONSTC2 4480.224717878304 

extern int N;
extern int SIZE;
extern double H;
extern double WAVE;
extern double G_DELTA;

//	double cal_interpgrad(int mode, double xj, double yj, double zj);
double cal_partial(int mode, double xi, double yi, double zi, double xj, double yj, double zj);
//	double cal_regpartial(double x1, double y1, double x2, double y2, double tau);

//	double cal_J(int mode, double x, double y, double z);
double cal_delta(double epsilon, double x);
double cal_nonhom(double zstarx, double zstary, double zstarz, double A[]);
double cal_f(int m, int l, double x);

//	double gradientdx(double x, double y, double z);
//	double gradientdy(double x, double y, double z);
//	double gradientdz(double x, double y, double z);
//	double secondpartialxx(double x, double y, double z);
//	double secondpartialyy(double x, double y, double z);
//	double secondpartialzz(double x, double y, double z);
//	double secondpartialxy(double x, double y, double z);
//	double secondpartialxz(double x, double y, double z);
//	double secondpartialyz(double x, double y, double z);
//	double laplace(double x, double y, double z);

//	double GaussianCurv(double x, double y, double z);
//	double MeanCurv(double x, double y, double z);


double cal_phi(double xi, double yi, double zi, double xj, double yj, double zj);
double phiimg(double xi, double yi, double zi, double xj, double yj, double zj);
double psi(double xj, double yj, double zj);
double psiimg(double xj, double yj, double zj);


//	double starx(double zx, double zy, double zz);
//	double stary(double zx, double zy, double zz);
//	double starz(double zx, double zy, double zz);

double distance(double x, double y, double z);

double PolyWeight(double t, double tau);
double SineWeight(double t, double tau);

#endif
