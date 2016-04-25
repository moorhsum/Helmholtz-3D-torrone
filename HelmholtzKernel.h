
#ifndef HELMHOLTZ_KERNEL_H_DEFINED
#define HELMHOLTZ_KERNEL_H_DEFINED

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
#include "constants.h"
#include "frame.h"
#include "tools.h"
#include "functions.h"
#include "interpolate3D.h"
#include "PDE.h"

extern int N, SIZE;
extern double H, HCUBED;
extern double WAVE;

int new_HelmholtzKernel(int size, double delta, double epsilon, double *zxlist, double *zylist, double *zzlist, \
	double *zxstarlist, double *zystarlist, double *zzstarlist, double *dist, double *gradx, double *grady, \
	double *gradz, double *curv, double *Jacobian, \
	double *PolyW, double *SineW, double *PolyResult, double *SineResult, \
	double **PolyKernelH, double **PolyKernelHimg, double **SineKernelH, double **SineKernelHimg);

int new_HelmholtzKernelCombo(int size, double delta, double epsilon, double ***level, double ***leveldx, \
	double ***leveldy, double ***leveldz, double ***levelmcurv, double ***levelgcurv, \
	double *zxlist, double *zylist, double *zzlist, \
	double *zxstarlist, double *zystarlist, double *zzstarlist, double *dist, double *gradx, double *grady, \
	double *gradz, double *Jacobian, double *regfactor, double *regfactorimg, \
	double *PolyW, double *SineW, double *PolyResult, double *SineResult, \
	double **PolyKernelH, double **PolyKernelHimg, double **SineKernelH, double **SineKernelHimg, double wavenum, \
	double eta);

int DoubleLayerKernel(int size, double epsilon, double ***level, double ***leveldx, double ***leveldy, double ***leveldz, \
	double ***levelmcurv, double ***levelgcurv, double *zxlist, double *zylist, double *zzlist, \
	double *zxstarlist, double *zystarlist, double *zzstarlist, double *dist, double *gradx, double *grady, \
	double *gradz, double *delta, double *Jacobian, double *regfactor, double *regfactorimg, \
	double *result, double **Kernel, double **Kernelimg, double wavenum);


#endif
