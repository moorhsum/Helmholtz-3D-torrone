
#ifndef DISTANCE_APPLICATION3D_H_INCLUDED
#define DISTANCE_APPLICATION3D_H_INCLUDED

#include "frame.h"
#include "finite_difference3D.h"
#include "weno.h"

double GaussianCurv(double ***dist, int i, int j, int k, int N, double H);
double MeanCurv(double ***dist, int i, int j, int k, int N, double H);
int BothCurvs(double ***dist, int i, int j, int k, int N, double H, double *meancurv, double *gausscurv);
int ASSIGN_GRADIENT(int N, double H, double ***dist, double ***graddx, double ***graddy, double ***graddz, int weno_mode);

#endif
