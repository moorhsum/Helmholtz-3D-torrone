
#ifndef PRINTSLICE_H_DEFINED
#define PRINTSLICE_H_DEFINED

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "PDE.h"
#include "functions.h"
#include "tools.h"
#include "frame.h"
#include "interpolate3D.h"
#include "constants.h"

#define FORCING 1	//	force to print the interior as irrelevant		//

extern int N, SIZE;
extern double H, HCUBED;
extern int N2;
extern double H2;
extern double WAVE;
extern double XSLICE, YSLICE, ZSLICE;

int print_apruz(int size, double ***level, double ***leveldx, double ***leveldy, double ***leveldz, \
	double ***levelmcurv, double ***levelgcurv);
int print_apruy(int size, double ***level, double ***leveldx, double ***leveldy, double ***leveldz, \
	double ***levelmcurv, double ***levelgcurv);
int print_aprux(int size, double ***level, double ***leveldx, double ***leveldy, double ***leveldz, \
	double ***levelmcurv, double ***levelgcurv);

#endif
