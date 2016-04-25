
#ifndef PDE_H_DEFINED
#define PDE_H_DEFINED


#define LAPLACE 0
#define DIRICHLET 0
#define INTERIOR_MODE 0

#define INCIDENTALWAVE 1//	if true, has incidental wave, and partialu and partialuimg would be different	//
			//	The exact solution will (usually) not be available, we do sanity check at far field	//
#define PLANEWAVEX -1.0
#define PLANEWAVEY 0.0
#define PLANEWAVEZ 0.0

#define FACTOR 2.0	//	the factor for epsilon = factor * h	//
#define REG_MODE 1	//	to use regularized partial or not	//
#define POLY_TEST 0	//	test polynomial weight
#define SINE_TEST 1	//	test sine weight
#define TRAD_ONESIDED 0	//	onesided tubular region for traditional way	//
#define COMBO_ONESIDED 0//	onesided tubular region for combo way	//
#define COMBO_INTERIOR 0

#define DIST_FILEMODE 1	//	Reading distance function from file	//
#define BOUNDARY_UTEST 0//	Test solutions on boundary	//
#define TEST_SINGLE_LAYER 1//	Test solving with single layer formulation	//
#define TEST_NEW_ESTIMATOR 0//	test combo approach

#define TAU_FACTOR 2.0	//	Regularization
#define COMBO 1

#endif
