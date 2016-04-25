

#ifndef INTERPOLATE3D_H_INCLUDED
#define INTERPOLATE3D_H_INCLUDED


#include <stdio.h>
#include "frame.h"
#include "finite_difference3D.h"
#include "distance_application3D.h"


int cal_interpgrad(int N, double H, double xj, double yj, double zj, double ***gradx, double ***grady, double ***gradz, \
		double *interpolatex, double *interpolatey, double *interpolatez);

double cal_curv( double ***dist, double x, double y, double z, int N, double H);
double cal_GaussianCurv(double ***dist, double x, double y, double z, int N, double H);

double cal_dist(double ***dist, double x, double y, double z, int N, double H);
int cal_BothCurvs(double x, double y, double z, double H, double ***mcurv, double ***gcurv, double *calmcurv, double *calgcurv);

#endif
