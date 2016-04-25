#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/resource.h>
#include <time.h>
#include <string.h>
#include <omp.h>

#include "frame.h"
#include "constants.h"
#include "BICGSTAB.h"
#include "HelmholtzKernel.h"
#include "PDE.h"
#include "tools.h"
#include "functions.h"
#include "printslice.h"

#define VERIFY 0
#define ABS_FDTEST 0	//	1 if the finite difference is set to be absolute distance, else a constant multiple of h

#define DRAW 0		//	print apru (values on the whole frame) or not		//
#define TEST_BESSEL 0	//	test Bessel function for accuracy			//

//	#define TAU_FACTOR 0.1	//	DETERMINES TAU = h * TAU_FACTOR		//
#define GFACTOR 1.0	//	G_delta = GFACTOR * h^GPOWER		//
#define GPOWER 1.0	//	G_delta = GFACTOR * h^GPOWER            //
#define FASTSIZEBOUND 100000	//	for matrix size smaller than this, store the whole thing for speed	//
#define TEST_CONDITION 0	//	print out the whole matrix to test on condition number	//

//	#define TEST_NEW_ESTIMATOR 1	//	test new estimator for solving PDE	//

	#define ABS_MODE 0		//	0 for k|nabla d|_1 H, 1 for absolute, 2 for kH	//
	#define ABS_EPSILON 0.295
	#define ABS_EPSILON0 0.3
	#define ABS_DELTA0 0.005
//	#define COMBO 1
#define DEFAULT_ETA -1.0		//	u = ( dG/dn - i eta G ) * beta , G = (-i/4) H0(k|x|)		//
#define DEFAULT_OMP_CACHE_CHUNKSIZE 16
#define OMP_SIZE_THRESHOLD 1000
#define DEFAULT_OMP_THREADS 12	//	default threads used for parallelism 

//	Fourrier coefficients	//
#define COEFFICIENTA { {0.0, 0.0,0.0,0.0}, {-7.0, 3.0, 0.0, 0.0}, {-5.0,3.0, 5.0,0.0}, {6.0,-9.0,7.0,1.0} }
#define COEFFICIENTB { {0.0, 0.0,0.0,0.0}, {0.0, 8.0, 0.0, 0.0}, {0.0, -5.0, -4.0, 0.0}, {0.0, 4.0, -4.0, 8.0}}

int TOTALPOINTS = 0;
int TOTALNEWPOINTS = 0;
int N = 20;
int SIZE = 21;
double H, HCUBED;
double G_DELTA;
int N2 = 301;
int SIZE2 = 302;
double H2;
double WAVE = 15.0;
double XSLICE = 0.0;
double YSLICE = 0.0;
double ZSLICE = 0.0;

double NEARFACTOR = 2.0, FARFACTOR = 3.0;
double EPS_ORDER = 0.5, EPSFACTOR = 6.5, DELTAFACTOR = 0.1, DELTA_ORDER = 1.0;
//	double EPSFACTOR2 = 0.5, DELTAFACTOR2 = 0.0625;
	double EPSFACTOR2 = 1.0, DELTAFACTOR2 = 0.0625;
double ABS_NEARDIST = 0.05, ABS_FARDIST = 0.1;


double ***DIST, ***DIST_DX, ***DIST_DY, ***DIST_DZ, ***MCURV, ***GCURV;

double *Density, *Densityimg, *newPolyDensity, *newPolyDensityimg, *newSineDensity, *newSineDensityimg;
double **newPolyKernelH, **newPolyKernelHimg, **newSineKernelH, **newSineKernelHimg;
double *newzx, *newzy, *newzz, *newzstarx, *newzstary, *newzstarz;
double *newdist, *newgradx, *newgrady, *newgradz, *newcurv, *newJacobian, *newRegfactor, *newRegfactorimg;
double *PolyW, *SineW, *PolyResult, *SineResult;
double *newBoundaryValue, *newBoundaryValueimg;

double get_C(int size, int pi, int pj);
double get_Ccomp(int size, int pi, int pj, double epsilon, double tempC[]);
double omp_get_Ccomp(int size, int pi, int pj, double epsilon, double tempC[]);
int print_Creg(int size);
int BiCGSTAB(int size);
int BiCGcomplex(int size, double epsilon);
int parallel_BiCGcomplex(int size, double epsilon, int chunk);


int print_f(int size, char *filename, char *filenameimg, double delta, double epsilon, int mode);
double cal_u(double x, double y, double z);
double cal_uimg(double x, double y, double z);
double cal_partialu(double x, double y, double z);
double cal_partialuimg(double x, double y, double z);
double cal_vdl(double x, double y, double z);
double cal_vdlimg(double x, double y, double z);

int print_kernel(int size, double samplex, double sampley, double samplez);
int NeumannTest(int size, double delta0, double epsilon0, double samplex, double sampley, double samplez);
int Combo_NeumannTest(int size, double samplex, double sampley, double samplez, double eta);

int main(int argc, char *argv[])
{
	int weno_mode = 5;
	int i, j, k;
	int adjust;
	double testx;
	int RESOLVE = 0, RESOLVE_NEW = 0;
	double delta0, epsilon0, eta;
	double thisdist, thisgradx, thisgrady;

	//	test erf	//
	//	printf("Test erf function. Enter x\n");
	//	scanf("%lf", &testx);
	//	printf("erf(%lf) = %lf.\n", testx, erf(testx));

	FILE *fpdensity, *fpdensityimg, *fpnewPolydensity, *fpnewPolydensityimg, *fpnewSinedensity, *fpnewSinedensityimg;
	FILE *fpPolyResult, *fpSineResult, *fpnewpointdata;
	FILE *fptestPoly, *fptestSine, *fptestx, *fptestb, *fptestPolyimg, *fptestSineimg, *fptestximg, *fptestbimg;

	int counter = 0, newcounter = 0;
	double zx, zy, zz;
//	double h = (INT_U - INT_L)/N;
	double epsilon = FACTOR * H;
	double samplex, sampley, samplez;
	double inJ0, outJ0, inY0, outY0;
	double temp, funcf, funcfder, coefp, coefq;
	time_t begin, end, seconds;

	char denfilename[50], denimgfilename[50], newPolydenfilename[50], newPolydenimgfilename[50];
	char newSinedenfilename[50], newSinedenimgfilename[50], PolyResultfilename[50], SineResultfilename[50];
	char newpointdatafilename[50];
	char filenameb[50], filenamebimg[50], newfilenameb[50], newfilenamebimg[50];

	char distfilename[50], verdistfilename[50], vergradx[50], vergrady[50], vergradz[50], vermcurv[50], vergcurv[50];
	FILE *fpdist, *fpdist2, *fpgradx, *fpgrady, *fpgradz, *fpmcurv, *fpgcurv;

	const rlim_t kStackSize = 512 * 1024 * 1024;   // min stack size = 128 MB 
	struct rlimit rl; 
	int result; 
	int start_nthreads, nthreads, using_nthreads;
	int chunk = DEFAULT_OMP_CACHE_CHUNKSIZE;


	nthreads = DEFAULT_OMP_THREADS;
	H = (double)(INT_U-INT_L)/(double)(N);
	HCUBED = H*H*H;
	G_DELTA = GFACTOR * pow(H, GPOWER);
	H2 = (double)(INT_U-INT_L)/(double)(N2);

	eta = DEFAULT_ETA;

//	printf("Here.\n");
	
//	omp_set_dynamic(0);
//	printf("Now has %d threads.\n", omp_get_num_threads());
	#pragma omp parallel
	{
		#pragma omp master
		{
			start_nthreads = omp_get_num_threads();
			printf("Originally has %d threads.\n", start_nthreads);
		}
//		#pragma omp master
//		{
//			omp_set_num_threads(nthreads);
//		}
//		#pragma omp master
//		{
//			using_nthreads = omp_get_num_threads();
//			printf("Set to use %d threads.\n", using_nthreads);
//		}
	}

	if (argc > 1)
	{
		N = atof(argv[1]);
		SIZE = N + 1;
		H = (INT_U - INT_L)/N;
		HCUBED = H*H*H;
		G_DELTA = GFACTOR * pow(H, GPOWER);
		if ( LAPLACE == 0 )
		{
			if (argc > 2)
			{
				WAVE = atof(argv[2]);
			}
			else
			{		
				printf("Enter Wave number k for the helmoholtz problem.\n");
				scanf("%lf", &WAVE);
			}
		}
	}

	if ( (LAPLACE == 0) && (N*PI/(2.0*WAVE) < 6.0) )
	{
		printf("Capturing one wave with roughly %lf points, which may not be enough and result in large error.\n", N*PI/(2.0*WAVE));
		printf("Ideally each wave should have more than six points. In this case, we need at least N = %lf.\n", 12.0*WAVE/PI);
		printf("Adjusting N?(1 for yes, 0 for no) : ");
		scanf("%d", &adjust);

		while (adjust!=0)
		{
			if (adjust==1)
			{
				printf("Please enter new N.\n");
				scanf("%d", &N);
				printf("Using N = %d.\n", (int)(N));
				if (N*PI/WAVE < 12.0)
				{
					printf("The new N = %d still may not be enough.\n", (int)(N));
					printf("Ideally each wave should have more than six points. \n");
					printf("In this case, we need at least N = %lf.\n", 12.0*WAVE/PI);
					printf("Adjusting N? (1 for yes/0 for no) : ");
					scanf("%d", &adjust);
					continue;
				}
				break;
			}
			else
			{
				printf("Invalid input. 1/0 only. Current N = %d\n", (int) (N));
				printf("Ideally each wave should have more than six points. \n");
				printf("In this case, we need at least N = %lf.\n", 12.0*WAVE/PI);
				printf("Adjusting N?(1 for yes/0 for no) : ");
				scanf("%d", &adjust);
			}
		}
	}
	if (TEST_BESSEL==1)
	{
		inJ0 = 0.0;
		printf("Input a number to test BesselJ0 functions. negative to exit.\n");
		while (inJ0 >= -0.0000001)
		{
			scanf("%lf", &inJ0);
			outJ0 = BesselJ0(inJ0);
			printf("Bessel function J0(%lf) is %.10lf.\n", inJ0, outJ0);
		}
	
		inY0 = 0.0;
		printf("Input a number to test BesselY0 functions. negative to exit.\n");
		while (inY0 >= -0.0000001)
		{
			scanf("%lf", &inY0);
			outY0 = BesselY0(inY0);
			printf("Bessel function Y0(%lf) is %.10lf.\n", inY0, outY0);
		}
	}
	result = getrlimit(RLIMIT_STACK, &rl);
	if (result == 0) 
	{ 
        	if (rl.rlim_cur < kStackSize) 
	        { 
	        	rl.rlim_cur = kStackSize; 
		        result = setrlimit(RLIMIT_STACK, &rl); 
		        if (result != 0) 
	  	      	{
		                fprintf(stderr, "setrlimit returned result = %d\n", result); 
		        }
		}
	}
//	printf("The stack size now is %d Mbytes.\n", int(rl.rlim_cur)/1048576 );

	if ( ( (DIRICHLET == 0) && (LAPLACE == 1) ) || ( (LAPLACE == 0) && (INTERIOR_MODE == 1) ) )
		printf("Check mode. Dirichlet = %d, Laplce = %d, Interior mode = %d.\n", DIRICHLET, LAPLACE, INTERIOR_MODE);

	if (argc < 6)
	{
		printf("Please enter the sample point to solve.\n");
		scanf("%lf%lf%lf", &samplex, &sampley, &samplez);
	}
	else
	{
		samplex = atof(argv[3]);
		sampley = atof(argv[4]);
		samplez = atof(argv[5]);

		if (argc == 8)
		{
			if (ABS_FDTEST == 1)
			{
				ABS_NEARDIST = atof(argv[6]);
				ABS_FARDIST = atof(argv[7]);
			}
			else
			{
				NEARFACTOR = atof(argv[6]);
				FARFACTOR = atof(argv[7]);
			}
		}
	}
	if (ABS_MODE == 1)
	{
//		epsilon = ABS_EPSILON;
		epsilon0 = ABS_EPSILON0;
		delta0 = ABS_DELTA0;
	}
	else
	{
//		epsilon = FACTOR * h;
		if (fabs(EPS_ORDER - 1.) < 1e-6)
			epsilon0 = EPSFACTOR * H;
		else
			epsilon0 = EPSFACTOR2 * pow(H, EPS_ORDER);
		if (fabs(DELTA_ORDER - 1.) < 1e-6 )
			delta0 = DELTAFACTOR * H;
		else
			delta0 = DELTAFACTOR2 * pow(H, DELTA_ORDER);
//		epsilon = epsilon0 - delta0;
	}
	

	printf("[delta0, epsilon0] = [%lf, %lf]\n", delta0, epsilon0);
	if ((TRAD_ONESIDED == 0) && (COMBO_ONESIDED == 1) )
	{
		coefp = 3.;
		if (COMBO_INTERIOR == 1)
			coefq = -0.5*( pow(1.-delta0,3.0) - pow(1.-epsilon0, 3.0) );
		else
			coefq = -0.5*( pow(1.+epsilon0,3.0) - pow(1.+delta0, 3.0) );
		epsilon = (epsilon0 - delta0)/2.0;
		temp = epsilon;
		for (i = 0; i < 10; i++)
		{
			funcf = pow(temp,3.0) + 3.*temp + coefq;
			funcfder = 3.*temp*temp + 3.;
			temp = temp - funcf/funcfder;
//			printf("epsilon = %lf, funcf = %lf, funcfder = %lf, residual = %.4E.\n", \
//			temp, funcf, funcfder, temp*temp*temp + 3.*temp + coefq);
		}
		epsilon = temp;
//		epsilon = sqrt((epsilon0*epsilon - delta0*delta0)/2.);
	}
	else
	{
		coefp = 3.;
		coefq = -1. * ( pow(1.+epsilon0,3.0) - pow(1.+delta0, 3.0) );
		epsilon = epsilon0 - delta0;
		temp = epsilon;
		for (i = 0; i < 10; i++)
		{
			funcf = pow(temp,3.0) + 3.*temp + coefq;
			funcfder = 3.*temp*temp + 3.;
			temp = temp - funcf/funcfder;
//			printf("epsilon = %lf, residual = %.4E.\n", temp, temp*temp*temp + 3.*temp + coefq);
		}
//		epsilon = temp;
		epsilon = sqrt(epsilon0*epsilon - delta0*delta0);
	}

////////////////////////////////////////////////////////////

	epsilon = FACTOR * H;

////////////////////////////////////////////////////////////

	printf("Using N = %d, H = %lf, [delta0, epsilon0] = [%lf, %lf], epsilon = %lf.\n", \
		N, H, delta0, epsilon0, epsilon);

	if ( LAPLACE == 0 )
	{
		if (INTERIOR_MODE==1)
			printf("Solving interior Helmholtz equation where k = %lf.\n", WAVE);
		else
			printf("Solving exterior Helmholtz equation where k = %lf.\n", WAVE);
		printf("This means using roughly %d points to capture a wave.\n", (int)(N*PI/(2.0*WAVE)));
	}
	else
	{
		if (INTERIOR_MODE==1)
			printf("Solving interior Laplace equation.\n");
		else
			printf("Solving exterior Laplace equation.\n");
	}
//	if (TEST_NEW_ESTIMATOR == 1)
//	{
//		printf("Delta = %lfH^%.2lf, epsilon = %lfH^%.2lf.\n", DELTAFACTOR, DELTA_ORDER, EPSFACTOR, EPS_ORDER);
//		printf("h = %lf, [delta0, epsilon0] = [%lf, %lf]\n", h, delta0, epsilon0);
//	}

	DIST = (double ***) malloc(SIZE * sizeof(double **));
	DIST_DX = (double ***) malloc(SIZE * sizeof(double **));
	DIST_DY = (double ***) malloc(SIZE * sizeof(double **));
	DIST_DZ = (double ***) malloc(SIZE * sizeof(double **));
	MCURV = (double ***) malloc(SIZE * sizeof(double **));
	GCURV = (double ***) malloc(SIZE * sizeof(double **));

//	#pragma omp parallel for private(i,j) schedule(static, 2)
	for ( i = 0; i < SIZE; i++)
	{
		DIST[i] = (double **) malloc(SIZE * sizeof(double *));
		DIST_DX[i] = (double **) malloc(SIZE * sizeof(double *));
		DIST_DY[i] = (double **) malloc(SIZE * sizeof(double *));
		DIST_DZ[i] = (double **) malloc(SIZE * sizeof(double *));
		MCURV[i] = (double **) malloc(SIZE * sizeof(double *));
		GCURV[i] = (double **) malloc(SIZE * sizeof(double *));
		for ( j = 0; j < SIZE; j++)
		{
			DIST[i][j] = (double *) malloc(SIZE * sizeof(double ));
			DIST_DX[i][j] = (double *) malloc(SIZE * sizeof(double ));
			DIST_DY[i][j] = (double *) malloc(SIZE * sizeof(double ));
			DIST_DZ[i][j] = (double *) malloc(SIZE * sizeof(double ));
			MCURV[i][j] = (double *) malloc(SIZE * sizeof(double ));
			GCURV[i][j] = (double *) malloc(SIZE * sizeof(double ));
		}
	}

	if (DIST_FILEMODE == 1)
	{
		sprintf(distfilename, "distance%d.txt", N);
		if ((fpdist = fopen(distfilename, "r")) == NULL)
		{
			printf("Open file %s error.\n", distfilename);
			exit(0);
		}
		for (i = 0; i < SIZE; i++)
		{
		for (j = 0; j < SIZE; j++)
		{
		for (k = 0; k < SIZE; k++)
		{
			fscanf(fpdist, "%lf", &DIST[i][j][k]);
		}
		}
		}
		fclose(fpdist);
	}
	else
	{
		#pragma omp parallel for private(i,j,k,zx,zy,zz) schedule(static, 2)
		for (i = 0; i < SIZE; i++)
		{
			zy = INT_L + (double)(i) * H;
		for (j = 0; j < SIZE; j++)
		{
			zx = INT_L + (double)(j) * H;
		for (k = 0; k < SIZE; k++)
		{
			zz = INT_L + (double)(k) * H;
			DIST[i][j][k] = distance(zx, zy, zz);
		}
		}
		}
	}

	ASSIGN_GRADIENT(N, H, DIST, DIST_DX, DIST_DY, DIST_DZ, weno_mode);


	for (i = 0; i < SIZE; i++)
	{
	for (j = 0; j < SIZE; j++)
	{
	for (k = 0; k < SIZE; k++)
	{
		BothCurvs(DIST, i, j, k, N, H, &MCURV[i][j][k], &GCURV[i][j][k]);
	}
	}
	}

	if (VERIFY == 1)
	{
		sprintf(verdistfilename, "vdistance%d.txt", N);
		sprintf(vergradx, "vgradx%d.txt", N);
		sprintf(vergrady, "vgrady%d.txt", N);
		sprintf(vergradz, "vgradz%d.txt", N);
		sprintf(vermcurv, "vmcurv%d.txt", N);
		sprintf(vergcurv, "vgcurv%d.txt", N);
	
		if ((fpdist2 = fopen(verdistfilename, "w")) == NULL)
		{
			printf("Open %s error.\n", verdistfilename);
			exit(0);
		}
		fpgradx = fopen(vergradx, "w");
		fpgrady = fopen(vergrady, "w");
		fpgradz = fopen(vergradz, "w");
		fpmcurv = fopen(vermcurv, "w");
		fpgcurv = fopen(vergcurv, "w");
		for (i = 0; i < SIZE; i++)
		{
		for (j = 0; j < SIZE; j++)
		{
		for (k = 0; k < SIZE; k++)
		{
			fprintf(fpdist2, "%.12lf ", DIST[i][j][k]);
			fprintf(fpgradx, "%.12lf ", DIST_DX[i][j][k]);
			fprintf(fpgrady, "%.12lf ", DIST_DY[i][j][k]);
			fprintf(fpgradz, "%.12lf ", DIST_DZ[i][j][k]);
			fprintf(fpmcurv, "%.12lf ", MCURV[i][j][k]);
			fprintf(fpgcurv, "%.12lf ", GCURV[i][j][k]);
		}
		}
			fprintf(fpdist2, "\n");
			fprintf(fpgradx, "\n");
			fprintf(fpgrady, "\n");
			fprintf(fpgradz, "\n");
			fprintf(fpmcurv, "\n");
			fprintf(fpgcurv, "\n");
		}
		fclose(fpdist2);
		fclose(fpgradx);
		fclose(fpgrady);
		fclose(fpgradz);
		fclose(fpmcurv);
		fclose(fpgcurv);
	}


//	#pragma omp parallel for private(i,j,k) schedule(static, 2)
//	for (i = 0; i < SIZE; i++)
//	{
//		double dxx, dyy, dzz, dxy, dxz, dyz;
//	for (j = 0; j < SIZE; j++)
//	{
//	for (k = 0; k < SIZE; k++)
//	{
//		dxx = XDER2(DIST, i, j, k, N, H);
///		dyy = YDER2(DIST, i, j, k, N, H);
//		dzz = ZDER2(DIST, i, j, k, N, H);
//		dxy = secondpartialxy(DIST, i, j, k, N, H);
//		dyz = secondpartialyz(DIST, i, j, k, N, H);
//		dxz = secondpartialxz(DIST, i, j, k, N, H);
//		MCURV[i][j][k] = -0.5 * (dxx + dyy + dzz);
//		GCURV[i][j][k] = dxx*dyy + dyy*dzz + dxx*dzz - (dxy*dxy + dyz*dyz + dxz*dxz);
//		if ( (DIST[i][j][k] != DIST[i][j][k]) || (DIST_DX[i][j][k] != DIST_DX[i][j][k]) || \
//			(DIST_DY[i][j][k] != DIST_DY[i][j][k]) || (DIST_DZ[i][j][k] != DIST_DZ[i][j][k]) || \
//			(MCURV[i][j][k] != MCURV[i][j][k]) || (GCURV[i][j][k]!=GCURV[i][j][k]) )
//		{
//			printf("At (%d, %d, %d), dist %lf, dx %lf, dy %lf, dz %lf, mcurv %lf, gcurv %lf.\n", \
//			i, j, k, DIST[i][j][k], DIST_DX[i][j][k], DIST_DY[i][j][k], DIST_DZ[i][j][k], \
//			MCURV[i][j][k], GCURV[i][j][k]);
//			exit(0);
//		}
//	}
//	}
//	}


	counter = 0;
	newcounter = 0;

//	if (WAVE < 8.0)
//		eta = -2./sqrt( PI*PI + 4.* pow(log(WAVE/2.) + EULER, 2.0));
//	else
//		eta = -1. * WAVE;

	printf("Eta is equal to %lf.\n", eta);



	for ( i = 0; i < SIZE; i++)
	{
//		zz = INT_L + H * (double)(i);
	for ( j = 0; j < SIZE; j++)
	{
//		zy = INT_L + H * (double)(j);
	for ( k = 0; k < SIZE; k++)
	{
//		zx = INT_L + H * (double)(k);
//		thisdist = distance(zx, zy, zz);
		thisdist = DIST[i][j][k];		

		if (TRAD_ONESIDED == 0)
		{
			if (fabs(thisdist) <= epsilon)
			{
			//		if (counter < 30)
			//			printf("zx = %.10lf, zy = %.10lf, zz = %.10lf.\n", zx, zy, zz);
				counter++;
			}
		}
		else
		{
			if ( (thisdist <= epsilon) && (thisdist >= 0.) )
				counter++;
		}
		if (COMBO_ONESIDED == 0)
		{
			if ( (fabs(thisdist) >= delta0) && (fabs(thisdist) <= epsilon0) )
			{
				newcounter++;
			}
		}
		else
		{
			if (COMBO_INTERIOR == 1)
			{
				if ( (thisdist >= delta0) && (thisdist<=epsilon0) )
					newcounter++;
			}
			else
			{
				if ( (thisdist <= -1.*delta0) && (thisdist>=-1.*epsilon0) )
					newcounter++;
			}
		}
	}
	}
	}
	printf("The matrix size is %d by %d, new matrix size is %d by %d.\n", counter, counter, newcounter, newcounter);
	TOTALPOINTS = counter;
	TOTALNEWPOINTS = newcounter;

////////////////////////////////////////////////////////////////////////////////////////////
//	return 0;
////////////////////////////////////////////////////////////////////////////////////////////

	begin = time(NULL);

	sprintf(filenameb, "%sbN%dk%d.txt", (ABS_MODE==1)?"Const_":"Rel_", (int)(N), (int)(WAVE));
	sprintf(filenamebimg, "%sbimgN%dk%d.txt", (ABS_MODE==1)?"Const_":"Rel_", (int)(N), (int)(WAVE));
	print_f(counter, filenameb, filenamebimg, -1., epsilon, 0);

//	printf("Here.\n");

	newBoundaryValue = (double *) malloc(newcounter * sizeof(double));
	newBoundaryValueimg = (double *) malloc(newcounter * sizeof(double));
	if (TEST_NEW_ESTIMATOR == 1)
	{
		sprintf(newfilenameb, "%snewbN%dk%d.txt", (ABS_MODE==1)?"Const_":"Rel_", (int)(N), (int)(WAVE));
		sprintf(newfilenamebimg, "%snewbN%dk%dimg.txt", (ABS_MODE==1)?"Const_":"Rel_", (int)(N), (int)(WAVE));
		print_f(newcounter, newfilenameb, newfilenamebimg, delta0, epsilon0, 1);

	}

	end = time(NULL);
	seconds = end - begin;
//	printf("Printing Boundary condition took %ld minutes %ld seconds.\n", seconds/60, seconds%60);


	Density = (double *) malloc(counter * sizeof(double));
	Densityimg = (double *) malloc(counter * sizeof(double));

	RESOLVE = 0;

	if (LAPLACE == 1)
	{
		sprintf(denfilename, "%sDensityN%d.txt", (ABS_MODE==1)?"Const_":"Rel_", (int)(N));
		sprintf(newPolydenfilename, "newPolyDensityN%ddel%.1lfeps%.1lf.txt", (int)(N), DELTAFACTOR, EPSFACTOR);
		sprintf(newSinedenfilename, "newSineDensityN%ddel%.1lfeps%.1lf.txt", (int)(N), DELTAFACTOR, EPSFACTOR);
		sprintf(PolyResultfilename, "PolyResultN%ddel%.1lfeps%.1lf.txt", (int)(N), DELTAFACTOR, EPSFACTOR);
		sprintf(SineResultfilename, "SineResultN%ddel%.1lfeps%.1lf.txt", (int)(N), DELTAFACTOR, EPSFACTOR);
		sprintf(newpointdatafilename, "newPointdataN%ddel%.1lfeps%.1lf.txt", (int)(N), DELTAFACTOR, EPSFACTOR);
	}
	else
	{
		sprintf(denfilename, "%sDensityN%dk%d.txt", (ABS_MODE==1)?"Const_":"Rel_", (int)(N), (int)(WAVE));
		sprintf(newPolydenfilename, "%snewPolyDensityN%dk%ddel%.1lfeps%.1lf.txt", \
			(ABS_MODE==1)?"Const_":"Rel_", (int)(N), (int)(WAVE), DELTAFACTOR, EPSFACTOR);
		sprintf(newSinedenfilename, "%snewSineDensityN%dk%ddel%.1lfeps%.1lf.txt", \
			(ABS_MODE==1)?"Const_":"Rel_", (int)(N), (int)(WAVE), DELTAFACTOR, EPSFACTOR);
		sprintf(PolyResultfilename, "%sPolyResultN%ddel%.1lfeps%.1lf.txt", \
			(ABS_MODE==1)?"Const_":"Rel_", (int)(N), DELTAFACTOR, EPSFACTOR);
		sprintf(SineResultfilename, "%sSineResultN%ddel%.1lfeps%.1lf.txt", \
			(ABS_MODE==1)?"Const_":"Rel_", (int)(N), DELTAFACTOR, EPSFACTOR);
		sprintf(newpointdatafilename, "%snewPointdataN%ddel%.1lfeps%.1lf.txt", \
			(ABS_MODE==1)?"Const_":"Rel_", (int)(N), DELTAFACTOR, EPSFACTOR);
	}
	

	if (TEST_CONDITION == 1)
	{
//		BiCGcomplex(counter);

		if (POLY_TEST == 1)
		{
			newPolyDensity = (double *) malloc(newcounter * sizeof(double));
			newPolyDensityimg = (double *) malloc(newcounter * sizeof(double));
			PolyW = (double *) malloc(newcounter * sizeof(double));
			PolyResult = (double *) malloc(newcounter * sizeof(double));

			newPolyKernelH = (double **) malloc(newcounter * sizeof(double *));
			newPolyKernelHimg = (double **) malloc(newcounter * sizeof(double *));

			for ( i = 0; i < newcounter; i++)
			{
				newPolyKernelH[i] = (double *) malloc(newcounter * sizeof(double));
				newPolyKernelHimg[i] = (double *) malloc(newcounter * sizeof(double));
			}
		}
		if (SINE_TEST == 1)
		{
			newSineDensity = (double *) malloc(newcounter * sizeof(double));
			newSineDensityimg = (double *) malloc(newcounter * sizeof(double));
			SineW = (double *) malloc(newcounter * sizeof(double));
			SineResult = (double *) malloc(newcounter * sizeof(double));
			newSineKernelH = (double **) malloc(newcounter * sizeof(double *));
			newSineKernelHimg = (double **) malloc(newcounter * sizeof(double *));
	
			for ( i = 0; i < newcounter; i++)
			{
				newSineKernelH[i] = (double *) malloc(newcounter * sizeof(double));
				newSineKernelHimg[i] = (double *) malloc(newcounter * sizeof(double));
			}
		}

		newzx = (double *) malloc(newcounter * sizeof(double));
		newzy = (double *) malloc(newcounter * sizeof(double));
		newzz = (double *) malloc(newcounter * sizeof(double));
		newzstarx = (double *) malloc(newcounter * sizeof(double));
		newzstary = (double *) malloc(newcounter * sizeof(double));
		newzstarz = (double *) malloc(newcounter * sizeof(double));
		newdist = (double * ) malloc(newcounter * sizeof(double));
		newgradx = (double *) malloc(newcounter * sizeof(double));
		newgrady = (double *) malloc(newcounter * sizeof(double));
		newgradz = (double *) malloc(newcounter * sizeof(double));
		newcurv = (double *) malloc(newcounter * sizeof(double));
		newJacobian = (double *) malloc(newcounter * sizeof(double));
		newRegfactor = (double *) malloc(newcounter * sizeof(double));
		newRegfactorimg = (double *) malloc(newcounter * sizeof(double));


		fptestb = fopen("b.txt", "w");
		fptestbimg = fopen("bimg.txt", "w");
		for ( i = 0; i < newcounter; i++)
		{
			fprintf(fptestb, "%.12lf\n", newBoundaryValue[i]);
			fprintf(fptestbimg, "%.12lf\n", newBoundaryValueimg[i]);
		}
		fclose(fptestb);
		fclose(fptestbimg);

		if (COMBO == 0)
		{
			new_HelmholtzKernel(newcounter, delta0, epsilon0, newzx, newzy, newzz, newzstarx, newzstary, \
			newzstarz, newdist, newgradx, newgrady, newgradz, newcurv, newJacobian, PolyW, SineW, \
			PolyResult, SineResult, newPolyKernelH, newPolyKernelHimg, newSineKernelH, newSineKernelHimg);
		}
		else
		{
			new_HelmholtzKernelCombo(newcounter, delta0, epsilon0, DIST, DIST_DX, DIST_DY, DIST_DZ, MCURV, \
			GCURV, newzx, newzy, newzz, newzstarx, newzstary, \
			newzstarz, newdist, newgradx, newgrady, newgradz, newJacobian, newRegfactor, \
			newRegfactorimg, PolyW, SineW, \
			PolyResult, SineResult, newPolyKernelH, newPolyKernelHimg, newSineKernelH, newSineKernelHimg, \
			WAVE, eta);
		}
		if (POLY_TEST == 1)
		{
			fptestPoly = fopen("Poly.txt", "w");
			fptestPolyimg = fopen("Polyimg.txt", "w");		
			for ( i = 0; i < newcounter; i++)
			{
				for ( j = 0; j < newcounter; j++)
				{
					fprintf(fptestPoly, "%lf\t", newPolyKernelH[i][j]);
					fprintf(fptestPolyimg, "%lf\t", newPolyKernelHimg[i][j]);
				}
				fprintf(fptestPoly, "\n");
				fprintf(fptestPolyimg, "\n");
			}
			fclose(fptestPoly);
			fclose(fptestPolyimg);
			for ( i = 1; i < newcounter; i++)
			{
				free(newPolyKernelH[i]);
				free(newPolyKernelHimg[i]);
			}
			free(newPolyDensity);
			free(newPolyDensityimg);
			free(PolyResult);
			free(PolyW);
			free(newPolyKernelH);
			free(newPolyKernelHimg);
		}
		if (SINE_TEST == 1)
		{
			fptestSine = fopen("Sine.txt", "w");
			fptestSineimg = fopen("Sineimg.txt", "w");
			for ( i = 0; i < newcounter; i++)
			{
				for ( j = 0; j < newcounter; j++)
				{
					fprintf(fptestSine, "%lf\t", newSineKernelH[i][j]);
					fprintf(fptestSineimg, "%lf\t", newSineKernelHimg[i][j]);
				}
				fprintf(fptestSine, "\n");
				fprintf(fptestSineimg, "\n");
			}
			fclose(fptestSine);
			fclose(fptestSineimg);
			for ( i = 1; i < newcounter; i++)
			{
				free(newSineKernelH[i]);
				free(newSineKernelHimg[i]);
			}
			free(newSineDensity);
			free(newSineDensityimg);
			free(SineResult);
			free(SineW);
			free(newSineKernelH);
			free(newSineKernelHimg);
		}
		free(newzx);
		free(newzy);
		free(newzz);
		free(newzstarx);
		free(newzstary);
		free(newzstarz);
		free(newgradx);
		free(newgrady);
		free(newgradz);
		free(newdist);
		free(newJacobian);
		free(newcurv);
		free(newRegfactor);
		free(newRegfactorimg);

		BiCGcomplex(counter, epsilon);
	}

//////////////////////////////////////////////////////////////////////////////////////////


	if (TEST_SINGLE_LAYER == 1)
	{
	if ((fpdensity = fopen( denfilename, "r"))==NULL)
	{
		printf("Open file %s error.\n", denfilename);
		if (LAPLACE == 1)
			printf("File %s does not exist.\n", denfilename);
		RESOLVE = 1;
	}
	if (LAPLACE == 0)
	{
		sprintf(denimgfilename, "%sDensityimgN%dk%d.txt", (ABS_MODE==1)?"Const_":"Rel_", (int)(N), (int)(WAVE));
		if ((fpdensityimg = fopen(denimgfilename, "r"))==NULL)
		{
//			printf("Open file %s error.\n", denimgfilename);
			printf("File %s does not exist.\n", denimgfilename);
			RESOLVE = 1;
		}
	}

	if (RESOLVE == 1)
	{
		begin = time(NULL);
		if (LAPLACE == 0)
		{
			if ( (counter < OMP_SIZE_THRESHOLD) || (counter > FASTSIZEBOUND) )
			{
				printf("Regular.\n");
				BiCGcomplex(counter, epsilon);
			}
			else
			{
				printf("Parallel.\n");
				parallel_BiCGcomplex(counter, epsilon, chunk);
			}
		}
//		BiCGSTABcomplex(counter);
		else
			BiCGSTAB(counter);

		end = time(NULL);
		seconds = end - begin;
		printf("BiCG took %ld minutes %ld seconds.\n", seconds/60, seconds%60);
	}
	else
	{
		for ( i = 0; i < counter; i++)
		{
			fscanf(fpdensity, "%lf", &Density[i]);
			if (LAPLACE == 0)
				fscanf(fpdensityimg, "%lf", &Densityimg[i]);
		}
		fclose(fpdensity);
		if (LAPLACE == 0)
			fclose(fpdensityimg);
	}


	if (DIRICHLET == 1)
		print_kernel(counter, samplex, sampley, samplez);
	
	NeumannTest(counter, 0., epsilon, samplex, sampley, samplez);

	free(Density);
	free(Densityimg);

	}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////



	if (TEST_NEW_ESTIMATOR == 1)
	{
		printf("Combo test begin.\n\n");
		if (POLY_TEST == 1)
		{
			if ( (fpnewPolydensity = fopen( newPolydenfilename, "r"))==NULL)
			{
				if (LAPLACE == 1)
					printf("File %s does not exist.\n", newPolydenfilename);
				fpnewPolydensity = fopen( newPolydenfilename, "w");
				RESOLVE_NEW = 1;
			}
			if ( (fpPolyResult = fopen( PolyResultfilename, "r")) == NULL)
			{
				printf("File %s does not exist.\n", PolyResultfilename);
				fpPolyResult = fopen( PolyResultfilename, "w");
				RESOLVE_NEW = 1;
			}
		}
		if (SINE_TEST == 1)
		{
			if ( (fpnewSinedensity = fopen( newSinedenfilename, "r"))==NULL)
			{
				if (LAPLACE == 1)
					printf("File %s does not exist.\n", newSinedenfilename);
				fpnewSinedensity = fopen( newSinedenfilename, "w");
				RESOLVE_NEW = 1;
			}
			if ( (fpSineResult = fopen( SineResultfilename, "r")) == NULL)
			{
				printf("File %s does not exist.\n", SineResultfilename);
				fpSineResult = fopen( SineResultfilename, "w");
				RESOLVE_NEW = 1;
			}
		}
		if ( (fpnewpointdata = fopen( newpointdatafilename, "r")) == NULL)
		{
			printf("File %s does not exist.\n", newpointdatafilename);
			fpnewpointdata = fopen(newpointdatafilename, "w");
			RESOLVE_NEW = 1;
		}
	
	}
	if (LAPLACE == 0)
	{

		if (TEST_NEW_ESTIMATOR == 1)
		{
			if (POLY_TEST == 1)
			{
				sprintf(newPolydenimgfilename, "%snewPolyDensityimgN%dk%ddel%.1lfeps%.1lf.txt", \
				(ABS_MODE==1)?"Const_":"Rel_", (int)(N), (int)(WAVE), DELTAFACTOR, EPSFACTOR);
				if ( (fpnewPolydensityimg = fopen(newPolydenimgfilename, "r")) == NULL)
				{
					printf("File %s does not exist.\n", newPolydenimgfilename);
					fpnewPolydensityimg = fopen(newPolydenimgfilename, "w");
					RESOLVE_NEW = 1;
				}
			}
			if (SINE_TEST == 1)
			{
				sprintf(newSinedenimgfilename, "%snewSineDensityimgN%dk%ddel%.1lfeps%.1lf.txt", \
				(ABS_MODE==1)?"Const_":"Rel_", (int)(N), (int)(WAVE), DELTAFACTOR, EPSFACTOR);
				if ( (fpnewSinedensityimg = fopen(newSinedenimgfilename, "r")) == NULL)
				{
					printf("File %s does not exist.\n", newSinedenimgfilename);
					fpnewSinedensityimg = fopen(newSinedenimgfilename, "w");
					RESOLVE_NEW = 1;
				}
			}
		}
	}


	if (TEST_NEW_ESTIMATOR == 1)
	{
		newzx = (double *) malloc(newcounter * sizeof(double));
		newzy = (double *) malloc(newcounter * sizeof(double));
		newzz = (double *) malloc(newcounter * sizeof(double));
		newzstarx = (double *) malloc(newcounter * sizeof(double));
		newzstary = (double *) malloc(newcounter * sizeof(double));
		newzstarz = (double *) malloc(newcounter * sizeof(double));
		newdist = (double * ) malloc(newcounter * sizeof(double));
		newgradx = (double *) malloc(newcounter * sizeof(double));
		newgrady = (double *) malloc(newcounter * sizeof(double));
		newgradz = (double *) malloc(newcounter * sizeof(double));
		newcurv = (double *) malloc(newcounter * sizeof(double));
		newJacobian = (double *) malloc(newcounter * sizeof(double));
		newRegfactor = (double *) malloc(newcounter * sizeof(double));
		newRegfactorimg = (double *) malloc(newcounter * sizeof(double));

		if (POLY_TEST == 1)
		{
			newPolyDensity = (double *) malloc(newcounter * sizeof(double));
			newPolyDensityimg = (double *) malloc(newcounter * sizeof(double));
			PolyW = (double *) malloc(newcounter * sizeof(double));
			PolyResult = (double *) malloc(newcounter * sizeof(double));
		}
		if (SINE_TEST == 1)
		{
			newSineDensity = (double *) malloc(newcounter * sizeof(double));
			newSineDensityimg = (double *) malloc(newcounter * sizeof(double));
			SineW = (double *) malloc(newcounter * sizeof(double));
			SineResult = (double *) malloc(newcounter * sizeof(double));
		}
		if (RESOLVE_NEW == 1)
		{
			begin = time(NULL);
			if (POLY_TEST == 1)
			{
				newPolyKernelH = (double **) malloc(newcounter * sizeof(double *));
				newPolyKernelHimg = (double **) malloc(newcounter * sizeof(double *));
				for ( i = 0; i < newcounter; i++)
				{
					newPolyKernelH[i] = (double *) malloc(newcounter * sizeof(double));
					newPolyKernelHimg[i] = (double *) malloc(newcounter * sizeof(double));
				}
			}
			if (SINE_TEST == 1)
			{
				newSineKernelH = (double **) malloc(newcounter * sizeof(double *));
				newSineKernelHimg = (double **) malloc(newcounter * sizeof(double *));
				for ( i = 0; i < newcounter; i++)
				{
					newSineKernelH[i] = (double *) malloc(newcounter * sizeof(double));
					newSineKernelHimg[i] = (double *) malloc(newcounter * sizeof(double));
				}
			}
			if (LAPLACE == 0)
			{

				if (COMBO == 0)
				{
				new_HelmholtzKernel(newcounter, delta0, epsilon0, newzx, newzy, newzz, newzstarx, \
				newzstary, newzstarz, newdist, newgradx, newgrady, newgradz, newcurv, newJacobian, \
				PolyW, SineW, PolyResult, SineResult, newPolyKernelH, newPolyKernelHimg, newSineKernelH, \
				newSineKernelHimg);
				}
				else
				{
				new_HelmholtzKernelCombo(newcounter, delta0, epsilon0, DIST, DIST_DX, DIST_DY, DIST_DZ, \
				MCURV, GCURV, newzx, newzy, newzz, newzstarx, newzstary, newzstarz, newdist, newgradx, \
				newgrady, newgradz, newJacobian, newRegfactor, newRegfactorimg, \
				PolyW, SineW, PolyResult, SineResult, newPolyKernelH, newPolyKernelHimg, newSineKernelH, \
				newSineKernelHimg, WAVE, eta);
				}

////////////////////////////////////////////////////////////////////////////////////////
/*
				fptestA = fopen("A.txt", "w");
				fptestb = fopen("b.txt", "w");
				fptestAimg = fopen("Aimg.txt", "w");
				fptestbimg = fopen("bimg.txt", "w");
				for (int i = 0; i < newcounter; i++)
				{
					for (int j = 0; j < newcounter; j++)
					{
						fprintf(fptestA, "%.12lf\t", newPolyKernelH[i][j]);
						fprintf(fptestAimg, "%.12lf\t", newPolyKernelHimg[i][j]);
					}
					fprintf(fptestb, "%.12lf\n", newBoundaryValue[i]);
					fprintf(fptestbimg, "%.12lf\n", newBoundaryValueimg[i]);
					fprintf(fptestA, "\n");
					fprintf(fptestAimg, "\n");
				}
				fclose(fptestA);
				fclose(fptestAimg);
				fclose(fptestb);
				fclose(fptestbimg);
*/
//////////////////////////////////////////////////////////////////////////////////////////

				for ( i = 0; i < newcounter; i++)
				{
					fprintf(fpnewpointdata, "%.12lf\t%.12lf\t%.12lf\t%.12lf\t%.12lf\t%.12lf\t%.12lf\t%.12lf\n", \
					newdist[i], newzstarx[i], newzstary[i], newzstarz[i], newgradx[i], newgrady[i], \
					newgradz[i], newRegfactor[i] );
				}
				fclose(fpnewpointdata);

				if (POLY_TEST == 1)
				{
					if (newcounter < OMP_SIZE_THRESHOLD)
					{
						CompBiCG(newcounter, newPolyDensity, newPolyDensityimg, newPolyKernelH, \
						newPolyKernelHimg, newBoundaryValue, newBoundaryValueimg);
					}
					else
					{
						omp_CompBiCG(newcounter, newPolyDensity, newPolyDensityimg, newPolyKernelH, \
						newPolyKernelHimg, newBoundaryValue, newBoundaryValueimg, chunk);
					}
					for ( i = 0; i < newcounter; i++)
					{
						fprintf(fpnewPolydensity, "%.12lf\n", newPolyDensity[i]);
						fprintf(fpnewPolydensityimg, "%.12lf\n", newPolyDensityimg[i]);
						fprintf(fpPolyResult, "%.12lf\n", PolyResult[i]);
						free(newPolyKernelH[i]);
						free(newPolyKernelHimg[i]);
					}
					fclose(fpnewPolydensity);
					fclose(fpnewPolydensityimg);
					fclose(fpPolyResult);

					free(newPolyKernelH);
					free(newPolyKernelHimg);
				}
				if (SINE_TEST == 1)
				{
					if (newcounter < OMP_SIZE_THRESHOLD)
					{
						CompBiCG(newcounter, newSineDensity, newSineDensityimg, newSineKernelH, \
						newSineKernelHimg, newBoundaryValue, newBoundaryValueimg);
					}
					else
					{
						omp_CompBiCG(newcounter, newSineDensity, newSineDensityimg, newSineKernelH, \
						newSineKernelHimg, newBoundaryValue, newBoundaryValueimg, chunk);
					}
					for ( i = 0; i < newcounter; i++)
					{
						fprintf(fpnewSinedensity, "%.12lf\n", newSineDensity[i]);
						fprintf(fpnewSinedensityimg, "%.12lf\n", newSineDensityimg[i]);	
						fprintf(fpSineResult, "%.12lf\n", SineResult[i]);
						free(newSineKernelH[i]);
						free(newSineKernelHimg[i]);
					}
					fclose(fpnewSinedensity);
					fclose(fpnewSinedensityimg);
					fclose(fpSineResult);
					free(newSineKernelH);
					free(newSineKernelHimg);
				}

//				fptestx = fopen("x.txt", "w");
//				fptestximg = fopen("ximg.txt", "w");

//					fprintf(fptestx, "%.12lf\n", newPolyDensity[i]);
//					fprintf(fptestximg, "%.12lf\n", newPolyDensityimg[i]);
//				fclose(fptestx);
//				fclose(fptestximg);

			}
			else
			{
				printf("Not now.\n");
				exit(0);
//				newBiCGSTAB
			}
			end = time(NULL);
			seconds = end-begin;
			printf("Solving new matrix took %ld minutes %ld seconds.\n", seconds/60, seconds%60);
		}
		else
		{
			for ( i = 0; i < newcounter; i++)
			{
				if (POLY_TEST == 1)
				{
					fscanf(fpnewPolydensity, "%lf", &newPolyDensity[i]);
					if (LAPLACE == 0)
						fscanf(fpnewPolydensityimg, "%lf", &newPolyDensityimg[i]);
					fscanf(fpPolyResult, "%lf", &PolyResult[i]);
				}
				if (SINE_TEST == 1)
				{
					fscanf(fpnewSinedensity, "%lf", &newSineDensity[i]);
					if (LAPLACE == 0)
						fscanf(fpnewSinedensityimg, "%lf", &newSineDensityimg[i]);
					fscanf(fpSineResult, "%lf", &SineResult[i]);
				}
				fscanf(fpnewpointdata, "%lf%lf%lf%lf%lf%lf%lf%lf", &newdist[i], &newzstarx[i], \
					&newzstary[i], &newzstarz[i], &newgradx[i], &newgrady[i], &newgradz[i], \
					&newRegfactor[i]);
			}
			if (POLY_TEST == 1)
			{
				fclose(fpnewPolydensity);
				fclose(fpPolyResult);
				if (LAPLACE == 0)
					fclose(fpnewPolydensityimg);
			}
			if (SINE_TEST == 1)
			{
				fclose(fpnewSinedensity);
				fclose(fpSineResult);
				if (LAPLACE == 0)
					fclose(fpnewSinedensityimg);
			}
			fclose(fpnewpointdata);
//			printf("Here.\n");
		}
	}
	if (TEST_NEW_ESTIMATOR == 1)
	{
		if (COMBO == 1)
			Combo_NeumannTest(newcounter, samplex, sampley, samplez, eta);
	}

	if (POLY_TEST == 1)
	{
		free(newPolyDensity);
		free(newPolyDensityimg);
	}
	if (SINE_TEST == 1)
	{
		free(newSineDensity);
		free(newSineDensityimg);
	}
	


/*
	if (DRAW==1)
	{
		printf("Drawing slice with z coordinate %lf.\n", ZSLICE);
		print_apruz(counter);
		if (argc > 4)
		{
			printf("Drawing slice with x cooridnate %lf.\n", XSLICE);
			print_aprux(counter);
			if (argc > 5)
			{
				printf("Drawing slice with y cooridnate %lf.\n", YSLICE);
				print_apruy(counter);
			}
		}
	}
*/

	return 0;
}

int print_f(int size, char *filename, char *filenameimg, double delta, double epsilon, int mode)
{
	double f[size], fimg[size], tempvalue[2];
	double zx, zy, zz, zstarx, zstary, zstarz;
	double thisdist;

	int i, j, k;
	int counter = 0;

	FILE *fpf, *fpvdl;
	FILE *fpfimg, *fpvdlimg;

	if ((fpf = fopen(filename, "w"))==NULL)
	{
		printf("Open file %s error.\n", filename);
//		printf("Open file 3Drealb.txt error.\n");
		return 0;
	}
	fprintf(fpf, "%d\n", size);
	if (LAPLACE == 0)
	{
		if ((fpfimg = fopen(filenameimg, "w"))==NULL)
		{
			printf("Open file %s error.\n", filenameimg);
//			printf("Open file large3Dimgb.txt error.\n");
			return 0;
		}
		fprintf(fpfimg, "%d\n", size);
	}


	for ( i = 0; i < SIZE; i++)
	{
		zy = INT_L + H * (double)(i);
	for ( j = 0; j < SIZE; j++)
	{
		zx = INT_L + H * (double)(j);
	for ( k = 0; k < SIZE; k++)
	{
		zz = INT_L + H * (double)(k);
//		thisdist = distance(zx, zy, zz);		
		thisdist = DIST[i][j][k];

		if (mode == 0)
		{
			if (TRAD_ONESIDED == 0)
			{
				if ( (fabs(thisdist) > epsilon) || (fabs(thisdist) < delta) )
					continue;
			}
			else
			{
				if ( (thisdist > epsilon) || (thisdist < delta) )
					continue;
			}
		}
		else
		{
			if (COMBO_ONESIDED == 0)
			{
				if ( ( fabs(thisdist) > epsilon) || (fabs(thisdist) < delta) )
					continue;
			}
			else
			{
				if (COMBO_INTERIOR == 1)
				{
					if ( (thisdist > epsilon) || (thisdist < delta) )
						continue;
				}
				else
				{
					if ( (thisdist < -1.*epsilon) || (thisdist > -1.*delta) )
						continue;
				}
			}
		}

//		zstarx = starx(zx, zy, zz);
//		zstary = stary(zx, zy, zz);
//		zstarz = starz(zx, zy, zz);
		zstarx = zx - thisdist * DIST_DX[i][j][k];
		zstary = zy - thisdist * DIST_DY[i][j][k];
		zstarz = zz - thisdist * DIST_DZ[i][j][k];

		if (fabs(cal_dist(DIST, zstarx, zstary, zstarz, N, H))> 0.01)
		{
			printf("Att (%lf, %lf, %lf), grad (%lf, %lf, %lf), dist %lf\n", zx, zy, zz, \
				DIST_DX[i][j][k], DIST_DY[i][j][k], DIST_DZ[i][j][k], thisdist);
			printf("pf Projected to (%lf, %lf, %lf), dist %lf.\n", zstarx, zstary, zstarz, \
				cal_dist(DIST, zstarx, zstary, zstarz, N, H));
			exit(0);
		}


		if (DIRICHLET == 1)
			f[counter] = cal_u(zstarx, zstary, zstarz);
		else
			f[counter] = cal_partialu(zstarx, zstary, zstarz);

		if (LAPLACE==0)
		{
			if (DIRICHLET == 1)
				fimg[counter] = cal_uimg(zstarx, zstary, zstarz);
			else
				fimg[counter] = cal_partialuimg(zstarx, zstary, zstarz);
		}

		if (NONHOM_MODE==1)
		{
			cal_nonhom(zstarx, zstary, zstarz, tempvalue);
			f[counter] -= tempvalue[0];	//	integral of psi and phi for inhomogeneous mode	//
			if (LAPLACE == 0)
				fimg[counter] -= tempvalue[1];
		}

		fprintf(fpf, "%.15lf\n", f[counter]);
		if (LAPLACE == 0)
			fprintf(fpfimg, "%.15lf\n", fimg[counter]);

		if (mode == 1)
		{
			newBoundaryValue[counter] = f[counter];
			newBoundaryValueimg[counter] = fimg[counter];
		}

		counter++;
	}
	}
	}
	fclose(fpf);
	if (LAPLACE==0)
		fclose(fpfimg);

	return 0;
}


int print_kernel(int size, double samplex, double sampley, double samplez)
{
	FILE *fpk, *fpkimg;
	FILE *fpkreg, *fpkimgreg;

	int i, j, k, indicator;
	int counti = 0;
	double zx, zy, zz, zstarx, zstary, zstarz, delta, thisdist, thisJ;
	double partial, partial_reg, temp, partialimg, partialimg_reg;
	double partialx, partialy, partialz, partialximg, partialyimg, partialzimg;
	double gradx, grady, gradz;

	static double regfactor, term1, term2, term3, term4, regfactorimg;
	static double epsilon;
	static double tau, threshold;
	epsilon = FACTOR * H;
	tau = H * TAU_FACTOR;
	threshold = tau * tau;

	if ((fpk = fopen("large3Dsamplekernel.txt", "w+"))==NULL)
	{
		printf("Open file 3Dsamplekernel.txt error.\n");
		return 0;
	}
	if ((fpkreg = fopen("large3Dsamplekernelreg.txt", "w+"))==NULL)
	{
		printf("Open file 3Dsamplekernelreg.txt error.\n");
		return 0;
	}

	if (LAPLACE == 0)
	{
		if ((fpkimg = fopen("large3Dsamplekernelimg.txt", "w+"))==NULL)
		{
			printf("Open file 3Dsamplekernelimg.txt error.\n");
			return 0;
		}
		if ((fpkimgreg = fopen("large3Dsamplekernelimgreg.txt", "w+"))==NULL)
		{
			printf("Open file samplekernelimgreg.txt error.\n");
			return 0;
		}
		fprintf(fpkimg, "%d\n", size);
		fprintf(fpkimgreg, "%d\n", size);
		fprintf(fpkimg, "%.15lf %.15lf %.15lf\n", samplex, sampley, samplez);
		fprintf(fpkimgreg, "%.15lf %.15lf %.15lf\n", samplex, sampley, samplez);
		if (INTERIOR_MODE == 0)		//	EXTERIOR MODE, u IS OUTSIDE	//
		{
			if (distance(samplex, sampley, samplez) <= 0.0)	//	OUTSIDE		//
			{
				fprintf(fpkimg, "%.15lf\n", cal_uimg(samplex, sampley, samplez));
				fprintf(fpkimgreg, "%.15lf\n", cal_uimg(samplex, sampley, samplez));
			}
			else					//	INSIDE		//
			{
				fprintf(fpkimg, "%.15lf\n", cal_vdlimg(samplex, sampley, samplez));
				fprintf(fpkimgreg, "%.15lf\n", cal_vdlimg(samplex, sampley, samplez));
			}
		}
		else				//	INTERIOR MODE, u IS INSIDE	//
		{
			if (distance(samplex, sampley, samplez) >= 0.0)	//	INSIDE		//
			{
				fprintf(fpkimg, "%.15lf\n", cal_uimg(samplex, sampley, samplez));
				fprintf(fpkimgreg, "%.15lf\n", cal_uimg(samplex, sampley, samplez));
			}
			else					//	OUTSIDE		//
			{
				fprintf(fpkimg, "%.15lf\n", cal_vdlimg(samplex, sampley, samplez));
				fprintf(fpkimgreg, "%.15lf\n", cal_vdlimg(samplex, sampley, samplez));
			}
		}
	}



	
	fprintf(fpk, "%d\n", size);
	fprintf(fpk, "%.15lf %.15lf %.15lf\n", samplex, sampley, samplez);
	fprintf(fpkreg, "%d\n", size);
	fprintf(fpkreg, "%.15lf %.15lf %.15lf\n", samplex, sampley, samplez);


	if (INTERIOR_MODE == 0)		//	EXTERIOR MODE, u IS OUTSIDE	//
	{
		if (distance(samplex, sampley, samplez) <= 0.0)	//	OUTSIDE		//
		{
			fprintf(fpk, "%.15lf\n", cal_u(samplex, sampley, samplez));
			fprintf(fpkreg, "%.15lf\n", cal_u(samplex, sampley, samplez));
		}
		else					//	INSIDE		//
		{
			fprintf(fpk, "%.15lf\n", cal_vdl(samplex, sampley, samplez));
			fprintf(fpkreg, "%.15lf\n", cal_vdl(samplex, sampley, samplez));
		}
	}
	else				//	INTERIOR MODE, u IS INSIDE	//
	{
		if (distance(samplex, sampley, samplez) >= 0.0)	//	INSIDE		//
		{
			fprintf(fpk, "%.15lf\n", cal_u(samplex, sampley, samplez));
			fprintf(fpkreg, "%.15lf\n", cal_u(samplex, sampley, samplez));
		}
		else					//	OUTSIDE		//
		{
			fprintf(fpk, "%.15lf\n", cal_vdl(samplex, sampley, samplez));
			fprintf(fpkreg, "%.15lf\n", cal_vdl(samplex, sampley, samplez));
		}
	}


	term1 = 1.0/(4.0*PI*tau*radius);		//	1/8pi*tau(2/r)	//
	term2 = (-5.0*tau)/(256.0*pow(radius,3.0)*PI);	//	-(1/pi)*(5/512)*(2/r^3)*tau	//
	term3 = (-50.0*tau)/(1536.0*pow(radius,3.0)*PI);//	-(1/pi)*(25/1536)*(2/r^3)*tau	//
	term4 = (pow(WAVE,2.0)*tau)/(24.0*PI*radius);	//	for Helmholtz 3D		//
	
	if (LAPLACE==1)
	{
		if (INTERIOR_MODE==1)
			regfactor = term1 + term2 + term3;
		else
			regfactor = -1.0 * (term1 + term2 + term3);
	}
	else
	{
		if (INTERIOR_MODE==1)
			regfactor = term1 + term2 + term3 + term4;
		else
			regfactor = -1.0 * (term1 + term2 + term3 + term4);
	}

	for ( i = 0; i < SIZE; i++)
	{
		zy = INT_L + H * (double)(i);
		for ( j = 0; j < SIZE; j++)		
		{
			zx = INT_L + H*(double)(j);
			for ( k = 0; k < SIZE; k++)
			{
				zz = INT_L + H * (double)(k);
				thisdist = DIST[i][j][k];

				if (fabs(thisdist) > epsilon)
					continue;
				counti++;
//				zstarx = starx(zx, zy, zz);
//				zstary = stary(zx, zy, zz);
//				zstarz = starz(zx, zy, zz);
				zstarx = zx - thisdist * DIST_DX[i][j][k];				
				zstary = zy - thisdist * DIST_DY[i][j][k];				
				zstarz = zz - thisdist * DIST_DZ[i][j][k];				

//				gradx = cal_interpgrad(1, zstarx, zstary, zstarz);
//				grady = cal_interpgrad(2, zstarx, zstary, zstarz);
//				gradz = cal_interpgrad(3, zstarx, zstary, zstarz);
				indicator = cal_interpgrad(N, H, zstarx, zstary, zstarz, DIST_DX, DIST_DY, DIST_DZ, \
					&gradx, &grady, &gradz);

				if (indicator == 1)
				{
					printf("At (%lf, %lf, %lf), grad (%lf, %lf, %lf), dist %lf\n",
					zx, zy, zz, DIST_DX[i][j][k], DIST_DY[i][j][k], DIST_DZ[i][j][k], thisdist);
					printf("Projected to (%lf, %lf, %lf), dist %lf\n", \
					zstarx, zstary, zstarz, cal_dist(DIST, zstarx, zstary, zstarz, N, H));
				}

				partialx = cal_partial(1, samplex, sampley, samplez, zstarx, zstary, zstarz);
				partialy = cal_partial(2, samplex, sampley, samplez, zstarx, zstary, zstarz);
				partialz = cal_partial(3, samplex, sampley, samplez, zstarx, zstary, zstarz);

				if (LAPLACE == 0)
				{
					partialximg = cal_partial(4, samplex, sampley, samplez, zstarx, zstary, zstarz);
					partialyimg = cal_partial(5, samplex, sampley, samplez, zstarx, zstary, zstarz);
					partialzimg = cal_partial(6, samplex, sampley, samplez, zstarx, zstary, zstarz);
				}


				//	for interior problem	//
				if (INTERIOR_MODE==1)
				{
					partial = -1.0*( gradx*partialx + grady*partialy + gradz*partialz);	//	outer normal
					if (LAPLACE == 0)
						partialimg = -1.0*( gradx*partialximg + grady*partialyimg + gradz*partialzimg );
				}
				
				//	for exterior problem	//
				else
				{
					partial = gradx*partialx + grady*partialy + gradz*partialz;		//	inner normal
					if (LAPLACE == 0)
						partialimg = gradx * partialximg + grady * partialyimg + gradz * partialzimg;
				}

				if ( ( (samplex-zstarx)*(samplex-zstarx) + (sampley-zstary)*(sampley-zstary) +(samplez-zstarz)*(samplez-zstarz)) < threshold )
				{
					if (INTERIOR_MODE==1)
					{
						partial_reg = regfactor;
						if (LAPLACE == 0)
							partialimg_reg = regfactorimg;
					}
					else
					{
						partial_reg = regfactor;
						if (LAPLACE == 0)
							partialimg_reg = regfactorimg;
					}
				}
				else
				{
					partial_reg = partial;
					if (LAPLACE == 0)
						partialimg_reg = partialimg;
				}

				delta = cal_delta(epsilon, distance(zx, zy, zz));
				thisJ = 1. + 2. * thisdist * MCURV[i][j][k] + thisdist * thisdist * GCURV[i][j][k];

				temp = HCUBED * partial * delta * thisJ;
				fprintf(fpk, "%.15lf ", temp);
				temp = HCUBED * partial_reg * delta * thisJ;
				fprintf(fpkreg, "%.15lf ", temp);
				if (LAPLACE == 0)
				{
					temp = HCUBED * partialimg * delta * thisJ;
					fprintf(fpkimg, "%.15lf ", temp);
					temp = HCUBED * partialimg_reg * delta * thisJ;
					fprintf(fpkimgreg, "%.15lf ", temp);
				}

				if (counti==size)
				{
					fprintf(fpk, "\n");
					fprintf(fpkreg, "\n");
					counti = 0;
				}
			}
		}
	}
	fclose(fpk);
	fclose(fpkreg);
	if (LAPLACE == 0)
	{
		fclose(fpkimg);
		fclose(fpkimgreg);
	}
	return 0;
}


int NeumannTest(int size, double delta0, double epsilon0, double samplex, double sampley, double samplez)
{
	FILE *fplog;
	char logfilename[50];
	int i, j, k;
	double zx, zy, zz, zstarx, zstary, zstarz, delta, thisJ, thisphi, thisphiimg, temp, thisnorm, result;
	double gradx, grady, gradz, gradnorm;
	double samplexstar, sampleystar, samplezstar, samplephi, samplephiimg, samplenorm, sampleu, sampleuimg, sampleunorm;
	double sampledx, sampledy, sampledz, sampledist;
	double nearx, neary, nearz, neardist, nearu, nearuimg, nearphi, nearphiimg, nearnorm, nearunorm;
	double farx, fary, farz, fardist, faru, faruimg, farphi, farphiimg, farnorm, farunorm;

	double fd_u, fd_uimg, fd_unorm;
	double tempvalue[2];
	double exactu, exactuimg, exactunorm, erroru, erroruimg, error, newerroru, newerroruimg, newerror;
	double exactnearu, exactnearuimg, exactnearunorm, exactfaru, exactfaruimg, exactfarunorm;
	double exactsampleu, exactsampleuimg, exactsampleunorm, sampleuerror, sampleuimgerror, sampleunormerror;

	int counti;
	double thisdist, thisgradx, thisgrady, thisgradz;
	double dxx, dyy, dzz, dxy, dxz, dyz, cur_meancurv, cur_gcurv;
	double tau, threshold, term1, term2, term3, term4, regfactor, regfactorimg;


	tau = H * TAU_FACTOR;
	threshold = tau * tau;

//	printf("Beginning H = %lf.\n", H);

	sprintf(logfilename, "%sResult_N%dk%d.txt", (ABS_MODE==1)?"Const_":"Rel_", N, (int)WAVE);

	exactsampleu = cal_u(samplex, sampley, samplez);
	exactsampleuimg = cal_uimg(samplex, sampley, samplez);
	exactsampleunorm = sqrt( exactsampleu*exactsampleu + exactsampleuimg*exactsampleuimg );

	if (BOUNDARY_UTEST == 1)
	{

		sampledist = cal_dist(DIST, samplex, sampley, samplez, N, H);
		cal_interpgrad(N, H, samplex, sampley, samplez, DIST_DX, DIST_DY, DIST_DZ, &sampledx, &sampledy, &sampledz);
		samplexstar = samplex - sampledist * sampledx;
		sampleystar = sampley - sampledist * sampledy;
		samplezstar = samplez - sampledist * sampledz;


		exactu = cal_u(samplexstar, sampleystar, samplezstar);
		exactuimg = cal_uimg(samplexstar, sampleystar, samplezstar);
		exactunorm = sqrt( exactu*exactu + exactuimg*exactuimg );


	cal_interpgrad(N, H, samplexstar, sampleystar, samplezstar, DIST_DX, DIST_DY, DIST_DZ, &gradx, &grady, &gradz);

		gradnorm = sqrt( gradx*gradx + grady*grady + gradz*gradz );

//		printf("sample point (%lf, %lf), projected to (%lf, %lf)\n", samplex, sampley, samplexstar, sampleystar);

		if (ABS_FDTEST == 1)
		{
			neardist = ABS_NEARDIST;
			fardist = ABS_FARDIST;	
		}
		else
		{
			neardist = NEARFACTOR * H;
			fardist = FARFACTOR * H;
		}


		nearx = samplexstar - neardist * gradx/gradnorm;
		neary = sampleystar - neardist * grady/gradnorm;
		nearz = samplezstar - neardist * gradz/gradnorm;

		farx = samplexstar - fardist * gradx/gradnorm;
		fary = sampleystar - fardist * grady/gradnorm;
		farz = samplezstar - fardist * gradz/gradnorm;

//		printf("Near point (%lf, %lf), far point (%lf, %lf)\n", nearx, neary, farx, fary);
		exactnearu = cal_u( nearx, neary, nearz );
		exactnearuimg = cal_uimg( nearx, neary, nearz);
		exactnearunorm = sqrt( exactnearu*exactnearu + exactnearuimg*exactnearuimg );

		exactfaru = cal_u( farx, fary, farz );
		exactfaruimg = cal_uimg( farx, fary, farz );
		exactfarunorm = sqrt( exactfaru*exactfaru + exactfaruimg*exactfaruimg );
	}


	counti = 0;
	nearu = 0.;
	nearuimg = 0.;
	faru = 0.;
	faruimg = 0.;
	sampleu = 0.;
	sampleuimg = 0.;

//	if (ABS_MODE != 1)
//		epsilon = FACTOR * h;
//	else
//		epsilon

	for ( i = 0; i < SIZE; i++)
	{
		zy = INT_L + H * i;
	for ( j = 0; j < SIZE; j++)		
	{
		zx = INT_L + H * j;
	for ( k = 0; k < SIZE; k++)
	{
		zz = INT_L + H * k;

//		printf("i, j, k = (%d, %d, %d), h = %lf\n", i, j, k, h);

//		if (ABS_MODE == 1)
//		{
//			epsilon = ABS_EPSILON;
//		}
//		else if (ABS_MODE == 0)
//		{
//
//			thisgradx = gradientdx(zx, zy, zz);
//			thisgrady = gradientdy(zx, zy, zz);
//			thisgradz = gradientdz(zx, zy, zz);
//			epsilon = ( fabs(thisgradx) + fabs(thisgrady) + fabs(thisgradz) ) * FACTOR * h;
//		}

//		thisdist = distance(zx, zy, zz);
		thisdist = DIST[i][j][k];		

		if (TRAD_ONESIDED == 0)
		{
			if (fabs(thisdist) > epsilon0)
				continue;
		}
		else
		{
			if ( (thisdist > epsilon0) || (thisdist < delta0) )
				continue;
		}

//		zstarx = starx(zx, zy, zz);
//		zstary = stary(zx, zy, zz);
//		zstarz = starz(zx, zy, zz);
		zstarx = zx - thisdist * DIST_DX[i][j][k];
		zstary = zy - thisdist * DIST_DY[i][j][k];
		zstarz = zz - thisdist * DIST_DZ[i][j][k];

		if (TRAD_ONESIDED == 0)
			delta = cal_delta(epsilon0, thisdist);
		else
			delta = 2. * cal_delta(epsilon0, thisdist);

//		dxx = secondpartialxx(zx, zy, zz);
//		dyy = secondpartialyy(zx, zy, zz);
//		dzz = secondpartialzz(zx, zy, zz);
//		dxy = secondpartialxy(zx, zy, zz);
//		dxz = secondpartialxz(zx, zy, zz);
//		dyz = secondpartialyz(zx, zy, zz);
//		cur_meancurv = -0.5 * (dxx + dyy + dzz);
//		cur_gcurv = dxx*dyy + dyy*dzz + dxx*dzz - dxy*dxy - dxz*dxz - dyz*dyz;
		thisJ = 1 + 2. * thisdist * MCURV[i][j][k] + thisdist*thisdist * GCURV[i][j][k];

//		thisJ = cal_J(J_M, zx, zy, zz);

		result = HCUBED * delta * thisJ;
//		if ((delta != delta) || (thisJ != thisJ) )
//		{
//			printf("Delta = %lf, J = %lf, point = (%lf, %lf, %lf), distance %lf.\n", \
//				delta, thisJ, zx, zy, zz, thisdist);
//			exit(0);
//		}


		if (BOUNDARY_UTEST == 1)
		{
			nearnorm = sqrt( (zstarx-nearx)*(zstarx-nearx) + (zstary-neary)*(zstary-neary) + \
					(zstarz - nearz)*(zstarz-nearz));
			farnorm = sqrt( (zstarx-farx)*(zstarx-farx) + (zstary-fary)*(zstary-fary) + \
					(zstarz - nearz)*(zstarz - nearz));
			if (nearnorm > tau)
			{
				nearphi = cal_phi( nearx, neary, nearz, zstarx, zstary, zstarz);
				nearphiimg = phiimg( nearx, neary, nearz, zstarx, zstary, zstarz);
			}
			else
			{
				nearphi = 0.;
				nearphiimg = 0.;
			}
			if (farnorm > tau)
			{
				farphi = cal_phi( farx, fary, farz, zstarx, zstary, zstarz);
				farphiimg = phiimg( farx, fary, farz, zstarx, zstary, zstarz);
			}
			else
			{
				farphi = 0.;
				farphiimg = 0.;
			}
			nearu += ( ( nearphi * Density[counti] - nearphiimg * Densityimg[counti] ) * result );
			nearuimg += ( ( nearphi * Densityimg[counti] + nearphiimg * Density[counti] ) * result );
			faru += ( ( farphi * Density[counti] - farphiimg * Densityimg[counti] ) * result );
			faruimg += ( ( farphi * Densityimg[counti] + farphiimg * Density[counti] ) * result );
		}

		samplenorm = sqrt( (zstarx-samplex)*(zstarx-samplex) + (zstary-sampley)*(zstary-sampley) + \
				(zstarz-samplez)*(zstarz-samplez) );
		if (samplenorm > tau)
		{
			samplephi = cos(WAVE * samplenorm)/(-4.*PI*samplenorm);
			samplephiimg = sin(WAVE * samplenorm)/(-4.*PI*samplenorm);
		}
		else
		{
			samplephi = 0.;
			samplephiimg = 0.;
		}

		sampleu += ( (samplephi * Density[counti] - samplephiimg * Densityimg[counti]) * result );
		sampleuimg += ( (samplephi * Densityimg[counti] + samplephiimg * Density[counti]) * result );


		if ( (nearu != nearu) || (faru != faru) || (sampleu != sampleu) )
//		if (counti > 1000)
		{
			printf("At %d, Density = (%lf, %lf), nearphi (%lf, %lf), farphi(%lf, %lf), result %lf\n", \
			counti, Density[counti], Densityimg[counti], nearphi, nearphiimg, farphi, farphiimg, \
			result);
			printf("(%lf, %lf, %lf), has distance %lf, delta = %lf, J = %lf, N = %d, H = %lf\n", \
			zx, zy, zz, thisdist, delta, thisJ, N, H);
			printf("sample point (%lf, %lf, %lf), phi (%lf, %lf), u = %lf, uimg = %lf.\n", \
			samplex, sampley, samplez, samplephi, samplephiimg, sampleu, sampleuimg);
//			exit(0);
		}

		counti++;
	}
	}
	}

	if (BOUNDARY_UTEST == 1)
	{
		nearunorm = sqrt(nearu*nearu + nearuimg*nearuimg);
		farunorm = sqrt(faru*faru + faruimg*faruimg);

		fd_u = nearu + neardist * (nearu - faru)/(fardist-neardist);
		fd_uimg = nearuimg + neardist * (nearuimg - faruimg)/(fardist-neardist);

		fd_unorm = sqrt(fd_u*fd_u + fd_uimg*fd_uimg);

		erroru = fabs(fd_u - exactu)/fabs(exactu);
		erroruimg = fabs(fd_uimg - exactuimg)/fabs(exactuimg);
		error = fabs(fd_unorm - exactunorm)/fabs(exactunorm);
	
//	printf("FarR: exact = %lf, appr = %lf, error %.2E\n", exactfaru, faru, fabs(exactfaru-faru)/fabs(exactfaru));
//	printf("FaIm: exact = %lf, appr = %lf, error %.2E\n", exactfaruimg, faruimg, \
//								fabs(exactfaruimg-faruimg)/fabs(exactfaruimg));
//	printf("FarN: exact = %lf, appr = %lf, error %.2E\n", exactfarunorm, farunorm, \
//							fabs(exactfarunorm-farunorm)/fabs(exactfarunorm));
//
//	printf("NeRe: exact = %lf, appr = %lf, error %.2E\n", exactnearu, nearu, fabs(exactnearu-nearu)/fabs(exactnearu));
//	printf("NeIm: exact = %lf, appr = %lf, error %.2E\n", exactnearuimg, nearuimg, \
//								fabs(exactnearuimg-nearuimg)/fabs(exactnearuimg));
//	printf("NeaN: exact = %lf, appr = %lf, error %.2E\n", exactnearunorm, nearunorm, \
//					fabs(exactnearunorm-nearunorm)/fabs(exactnearunorm));

//	printf("Real: exact = %lf, fd_u = %lf, error %.2E\n", exactu, fd_u, erroru);
//	printf("Imag: exact = %lf, fd_u = %lf, error %.2E\n", exactuimg, fd_uimg, erroruimg);
//	printf("Norm: exact = %lf, fd_u = %lf, error %.2E\n", exactunorm, fd_unorm, error);
	
	}

	sampleunorm = sqrt(sampleu*sampleu + sampleuimg*sampleuimg);
	sampleuerror = fabs(sampleu - exactsampleu)/fabs(exactsampleu);
	sampleuimgerror = fabs(sampleuimg - exactsampleuimg)/fabs(exactsampleuimg);
	sampleunormerror = fabs(sampleunorm - exactsampleunorm)/fabs(exactsampleunorm);
	printf("Real: exact = %lf, trad_approxu = %lf, error %.2E\n", exactsampleu, sampleu, sampleuerror);
	printf("Imag: exact = %lf, trad_approxu = %lf, error %.2E\n", exactsampleuimg, sampleuimg, sampleuimgerror);
	printf("Norm: exact = %lf, trad_approxu = %lf, error %.2E\n", exactsampleunorm, sampleunorm, sampleunormerror);

	if ((fplog = fopen(logfilename, "w"))==NULL)
	{
		printf("Open %s error.\n", logfilename);
		exit(0);
	}
	else
	{
	fprintf(fplog, "Real: exact = %lf, trad_approxu = %lf, error %.2E\n", exactsampleu, sampleu, sampleuerror);
	fprintf(fplog,"Imag: exact = %lf, trad_approxu = %lf, error %.2E\n", exactsampleuimg, sampleuimg, sampleuimgerror);
	fprintf(fplog,"Norm: exact = %lf, trad_approxu = %lf, error %.2E\n", exactsampleunorm, sampleunorm, sampleunormerror);
	fclose(fplog);
	}

	return 0;
}



///	EVALUATE U AT A POINT OFF OR ON INTERFACE.
//	IF A POINT IS ON INTERFACE, BOTH DOUBLE AND SINGLE LAYER USE THE DYNAMIC ESTIMATOR (X CHANGES WITH Y IN INTEGRAL)
//	IF A POINT IS OFF INTERFACE, DO NOT USE DYNAMIC ESTIMATOR (X STAYS FIXED FOR ANY Y IN INTEGRAL)

int Combo_NeumannTest(int size, double samplex, double sampley, double samplez, double eta)
{
	FILE *fplog;
	char logfilename[50];
	int i, j, k;
	int regcount;
	double zx, zy, zz, zstarx, zstary, zstarz, delta, thisJ, thisphi, thisphiimg, temp, thisnorm, result;
	double samplexstar, sampleystar, samplezstar, samplegradx, samplegrady, samplegradz, samplegradnorm;
	double *newx1, *newx2, *newx3, *xminusy1, *xminusy2, *xminusy3, *xminusynorm, *samplekernel, *samplekernelimg;
	double threshold = 0.01*H*H;
	double polysum, polysumimg, sinesum, sinesumimg, polyu, polyuimg, polyunorm, sineu, sineuimg, sineunorm;
	double polyerroru, polyerroruimg, polyerrorunorm, sineerroru, sineerroruimg, sineerrorunorm;
	double exactu, exactuimg, exactunorm, exactustar, exactustarimg, exactustarnorm;


	double term1, term2, term3, term4, regfactor, regfactorimg, tau;

	newx1 = (double *) malloc(size * sizeof(double));
	newx2 = (double *) malloc(size * sizeof(double));
	newx3 = (double *) malloc(size * sizeof(double));
	xminusy1 = (double *) malloc(size * sizeof(double));
	xminusy2 = (double *) malloc(size * sizeof(double));
	xminusy3 = (double *) malloc(size * sizeof(double));
	samplekernel = (double *) malloc(size * sizeof(double));
	samplekernelimg = (double *) malloc(size * sizeof(double));

//	threshold = 0.0001*h*h;
//	tau = 0.01*h;

	sprintf(logfilename, "%sResult_N%dk%dCombo.txt", (ABS_MODE==1)?"Const_":"Rel_", N, (int)WAVE);

	tau = TAU_FACTOR * H;
	threshold = tau * tau;


//	term1 = 1.0/(4.0*PI*tau*radius);		//	1/8pi*tau(2/r)	//
//	term2 = (-5.0*tau)/(256.0*pow(radius,3.0)*PI);	//	-(1/pi)*(5/512)*(2/r^3)*tau	//
//	term3 = (-50.0*tau)/(1536.0*pow(radius,3.0)*PI);//	-(1/pi)*(25/1536)*(2/r^3)*tau	//
//	term4 = (pow(WAVE,2.0)*tau)/(24.0*PI*radius);	//	Helmholtz 3D, (k^2/48)*(k1+k2)*tau	//

//	if (INTERIOR_MODE==1)
//		regfactor = term1 + term2 + term3 + term4;
//	else
//		regfactor = -1.0*(term1 + term2 + term3 + term4);
//	regfactor = 0.0;		
//	regfactorimg = 0.0;




//	samplexstar = starx(samplex, sampley, samplez);
//	sampleystar = stary(samplex, sampley, samplez);
//	samplezstar = starz(samplex, sampley, samplez);
//
//	samplegradx = cal_interpgrad(1, samplexstar, sampleystar, samplezstar);
//	samplegrady = cal_interpgrad(2, samplexstar, sampleystar, samplezstar);
//	samplegradz = cal_interpgrad(3, samplexstar, sampleystar, samplezstar);
//	samplegradnorm = sqrt( samplegradx*samplegradx + samplegrady*samplegrady + samplegradz*samplegradz);
//	samplegradx = samplegradx / samplegradnorm;
//	samplegrady = samplegrady / samplegradnorm;
//	samplegradz = samplegradz / samplegradnorm;


	exactu = cal_u(samplex, sampley, samplez);
	exactuimg = cal_uimg(samplex, sampley, samplez);
	exactunorm = sqrt( exactu*exactu + exactuimg*exactuimg );

//	exactustar = cal_u(samplexstar, sampleystar, samplezstar);
//	exactustarimg = cal_uimg(samplexstar, sampleystar, samplezstar);
//	exactustarnorm = sqrt( exactu*exactu + exactuimg*exactuimg );

	//	OFF INTERFACE, u(x) = [ dG/dny(x,y*) + i * eta * G(x,y*) ] W(y) J(y) beta(y) H^2	//
	regcount = 0;
	#pragma omp parallel for private(i) schedule(static, 16) reduction(+:regcount)
	for ( i = 0; i < size; i++)
	{
		double xmy1, xmy2, xmy3;
		double partialx, partialy, partialz, partialximg, partialyimg, partial, partialimg;
		double CosineValue, SineValue, inner_pd, fourpinorm, fourpinorm2, parameter;
		double singlereal, singleimg, doublereal, doubleimg;

		xmy1 = samplex - newzstarx[i];
		xmy2 = sampley - newzstary[i];
		xmy3 = samplez - newzstarz[i];

		thisnorm = sqrt( xmy1*xmy1 + xmy2*xmy2 + xmy3*xmy3 );
		parameter = WAVE * thisnorm;

		//	cos(k|x-y*|)	//
		CosineValue = cos(parameter);
		//	sin(k|x-y*|)	//
		SineValue = sin(parameter);

		inner_pd = xmy1 * newgradx[i] + xmy2 * newgrady[i] + xmy3 * newgradz[i];
		fourpinorm = 4. * PI * thisnorm;
		fourpinorm2 = fourpinorm * thisnorm;

		//	dG/dny = -1 * Grad G dot grad d = 
		//	REAL: 1/(4pi|x-y|^2) * (cos(k|x-y|)/|x-y| + ksin(k|x-y|)) * [(x-y) dot grad d]	//
		//	IMAG: 1/(4pi|x-y|^2) * (sin(k|x-y|)/|x-y| - kcos(k|x-y|)) * [(x-y) dot grad d]	//
		//	There's a negative from H0' = -H1, one from ny = -grad d, one from d|x-y*|/dy* = -(x-y*)/|x-y*|

		if (thisnorm > tau)
		{
			partial = (CosineValue/thisnorm + WAVE * SineValue) * inner_pd / fourpinorm2;
			partialimg = (SineValue/thisnorm - WAVE * CosineValue) * inner_pd / fourpinorm2;
		}
		else
		{
			if (DIRICHLET == 1)
				partial = newRegfactor[j];
			else
				partial = newRegfactor[i];
			partialimg = 0.0;
			regcount += 1;
		}


		//	-i * eta * G(x,y) = -i eta * (-i/4)*(J0 + iY0) = 0.25*eta*(J0+iY0)	//
		//	-i * eta * G(x,y) = -i eta * -1.0 * ( cos(k|x-y|) + isin(k|x-y|) )/ 4PI|x-y|	//
		//			  = eta * ( -sin(k|x-y|) + i cos(k|x-y|) ) / 4PI|x-y|	//
		singlereal = eta * SineValue / fourpinorm;
		singleimg = -1. * eta * CosineValue / fourpinorm;
//		singlereal = -1. * eta * SineValue / fourpinorm;
//		singleimg = eta * CosineValue / fourpinorm;

		samplekernel[i] = partial + singlereal;
		samplekernelimg[i] = partialimg + singleimg;
	}

	polysum = 0.0;
	polysumimg = 0.0;
	sinesum = 0.0;
	sinesumimg = 0.0;

	#pragma omp parallel for private(i) schedule(static, 16) reduction(+:polysum,polysumimg,sinesum,sinesumimg)
	for ( i = 0 ; i < size; i++)
	{
		if (POLY_TEST == 1)
		{
			polysum += ( (samplekernel[i] * newPolyDensity[i] - samplekernelimg[i] * newPolyDensityimg[i] ) * \
				PolyResult[i] );
			polysumimg += ( (samplekernelimg[i] * newPolyDensity[i] + samplekernel[i] * newPolyDensityimg[i]) * \
				PolyResult[i] );
		}
		if (SINE_TEST == 1)
		{
			sinesum += ( (samplekernel[i] * newSineDensity[i] - samplekernelimg[i] * newSineDensityimg[i] ) * \
				SineResult[i] );
			sinesumimg += ( (samplekernelimg[i] * newSineDensity[i] + samplekernel[i] * newSineDensityimg[i]) * \
				SineResult[i] );
		}
	}

	if (POLY_TEST == 1)
	{
		polyu = polysum;
		polyuimg = polysumimg;
		polyunorm = sqrt( polyu * polyu + polyuimg * polyuimg );
		polyerroru = fabs(polyu - exactu)/fabs(exactu);
		polyerroruimg = fabs(polyuimg - exactuimg)/fabs(exactuimg);
		polyerrorunorm = fabs(polyunorm - exactunorm)/fabs(exactunorm);
	}
	if (SINE_TEST == 1)
	{
		sineu = sinesum;
		sineuimg = sinesumimg;
		sineunorm = sqrt( sineu * sineu + sineuimg * sineuimg );
		sineerroru = fabs(sineu - exactu)/fabs(exactu);
		sineerroruimg = fabs(sineuimg - exactuimg)/fabs(exactuimg);
		sineerrorunorm = fabs(sineunorm - exactunorm)/fabs(exactunorm);
	}

	printf("Regcount = %d.\n", regcount);
//	printf("Poly:\n");
		
	printf("Real: exact = %lf, ", exactu);
	if (POLY_TEST == 1)
		printf("polyu %lf, error %.2E, ", polyu, polyerroru);
	if (SINE_TEST == 1)
		printf("sineu %lf, error %.2E\n", sineu, sineerroru);
	printf("Imag: exact = %lf, ", exactuimg);
	if (POLY_TEST == 1)
		printf("polyu %lf, error %.2E, ", polyuimg, polyerroruimg);
	if (SINE_TEST == 1)
		printf("sineu %lf, error %.2E\n", sineuimg, sineerroruimg);
	printf("Norm: exact = %lf, ", exactunorm);
	if (POLY_TEST == 1)
		printf("polyu %lf, error %.2E, ", polyunorm, polyerrorunorm);
	if (SINE_TEST == 1)
		printf("sineu %lf, error %.2E\n", sineunorm, sineerrorunorm);

	free(newx1);
	free(newx2);
	free(newx3);
	free(xminusy1);
	free(xminusy2);
	free(xminusy3);
	free(samplekernel);
	free(samplekernelimg);

	if ((fplog = fopen(logfilename, "w"))==NULL)
	{
		printf("Open %s error.\n", logfilename);
		exit(0);
	}
	else
	{
		fprintf(fplog, "Regcount = %d.\n", regcount);
	//	printf("Poly:\n");
		
		fprintf(fplog, "Real: exact = %lf, ", exactu);
		if (POLY_TEST == 1)
			fprintf(fplog, "polyu %lf, error %.2E, ", polyu, polyerroru);
		if (SINE_TEST == 1)
			fprintf(fplog, "sineu %lf, error %.2E\n", sineu, sineerroru);
		fprintf(fplog, "Imag: exact = %lf, ", exactuimg);
		if (POLY_TEST == 1)
			fprintf(fplog, "polyu %lf, error %.2E, ", polyuimg, polyerroruimg);
		if (SINE_TEST == 1)
			fprintf(fplog, "sineu %lf, error %.2E\n", sineuimg, sineerroruimg);
		fprintf(fplog, "Norm: exact = %lf, ", exactunorm);
		if (POLY_TEST == 1)
			fprintf(fplog, "polyu %lf, error %.2E, ", polyunorm, polyerrorunorm);
		if (SINE_TEST == 1)
			fprintf(fplog, "sineu %lf, error %.2E\n", sineunorm, sineerrorunorm);
		fclose(fplog);
	}


	return 0;
		

/*
		
		//	nx dot ny	//
		inner_pd1 = samplegradx * newgradx + samplegrady * newgrady;
		//	nx dot (x-y*)	//
		inner_pd2 = samplegradx * xmy1 + samplegrady * xmy2;
		//	ny dot (x-y*)	//
		inner_pd3 = xmy1 * newgradx + xmy2 * newgrady;

		//	k|x-y*|		//
		parameter = WAVE * thisnorm;
		commonterm1 = inner_pd1 / thisnorm;
		commonterm2 = inner_pd2 * inner_pd3 / (thisnorm * thisnorm);
		BesselY1value = BesselY1( parameter );
		BesselJ1value = BesselJ1( parameter );
		//	first term
	}

	for (int i = 0; i < size; i++)
	{
		newx1[i] = samplexstar - samplegradx * fabs(newdist[i]);
		newx2[i] = sampleystar - samplegrady * fabs(newdist[i]);

		xminusy1[i] = newx1[i] - newzstarx[i];
		xminusy2[i] = newx2[i] - newzstary[i];
		thisnorm = sqrt( xminusy1[i]*xminusy1[i] + xminusy2[i]*xminusy2[i] );

		//	nx dot ny	//
		inner_pd1 = samplegradx * newgradx[i] + samplegrady * newgrady[i];
		//	(x-y*) dot nx	//
		inner_pd2 = samplegradx * xminusy1[i] + samplegrady * xminusy2[i];
		//	(x-y*) dot ny	//
		inner_pd3 = newgradx[i] * xminusy1[i] + newgrady[i] * xminusy2[i];

		parameter = thisnorm * WAVE;
		commonterm1 = inner_pd1 / thisnorm;
		commonterm2 = inner_pd2 * inner_pd3 / (thisnorm * thisnorm);
		BesselY1value = BesselY1( parameter );
		BesselJ1value = BesselJ1( parameter );

		//	first term, H1(k|x-y*|) * [nx dot ny] / |x-y*|	(divided by k/2) 	//
		term1R = BesselY1value * commonterm1;
		term1I = BesselJ1value * commonterm1;
		//	second term, k * H0(k|x-y*|) * [(x-y*) dot nx] * [(x-y*) dot ny]/ |x-y*|^2	//
		term2R = BesselY0(parameter) * WAVE * commonterm2;
		term2I = BesselJ0(parameter) * WAVE * commonterm2;
		//	third term, H1(k|x-y*|) * [(x-y*) dot nx] * [(x-y*) dot ny] / |x-y*|^3
		term3R = BesselY1value * commonterm2 / xminusynorm;
		term3I = BesselJ1value * commonterm2 / xminusynorm;

		


		if (thisnorm > threshold)
		{
			singlereal = BesselY0(WAVE * thisnorm)/4.;
			singleimg = BesselY0(WAVE * thisnorm)/(-4.);

			doublereal = WAVE * (term1R + term2R + term3R)/2.;
			doubleimg = WAVE * (term1I + term2I + term3I)/(-2.);
			

			samplekernel[i] = BesselY0(WAVE * thisnorm)/4.;
			samplekernelimg[i] = BesselJ0(WAVE * thisnorm)/(-4.);
		}
		else
		{
			samplekernel[i] = 0.;
			samplekernelimg[i] = 0.;
		}
	}

	polysum = 0.;
	polysumimg = 0.;
	sinesum = 0.;
	sinesumimg = 0.;
	for (int i = 0; i < size; i++)
	{
		polysum +=  (( samplekernel[i] * newPolyDensity[i] - samplekernelimg[i] * newPolyDensityimg[i] )  * \
			PolyResult[i] );
		polysumimg += ( ( samplekernelimg[i] * newPolyDensity[i] + samplekernel[i] * newPolyDensityimg[i] ) * \
			PolyResult[i] );
		sinesum += ( ( samplekernel[i] * newSineDensity[i] - samplekernelimg[i] * newSineDensityimg[i] ) * \
			SineResult[i] );
		sinesumimg += ( ( samplekernelimg[i] * newSineDensity[i] + samplekernel[i] * newSineDensityimg[i] ) * \
			SineResult[i] );
//		printf("%d, Den %lf, Denimg %lf, PolyResult %lf\n",\
//			 i, newPolyDensity[i], newPolyDensityimg[i], PolyResult[i]);
	}
	polyu = polysum;
	polyuimg = polysumimg;
	polyunorm = sqrt( polyu * polyu + polyuimg * polyuimg );
	sineu = sinesum;
	sineuimg = sinesumimg;
	sineunorm = sqrt( sineu * sineu + sineuimg * sineuimg );

	polyerroru = fabs(polyu - exactu)/fabs(exactu);
	polyerroruimg = fabs(polyuimg - exactuimg)/fabs(exactuimg);
	polyerrorunorm = fabs(polyunorm - exactunorm)/fabs(exactunorm);

	sineerroru = fabs(sineu - exactu)/fabs(exactu);
	sineerroruimg = fabs(sineuimg - exactuimg)/fabs(exactuimg);
	sineerrorunorm = fabs(sineunorm - exactunorm)/fabs(exactunorm);
	
//	printf("Poly:\n");
	printf("Real: exact = %lf, polyu %lf, error %.2E, sineu %lf, error %.2E\n", \
	exactu, polyu, polyerroru, sineu, sineerroru);
	printf("Imag: exact = %lf, polyu %lf, error %.2E, sineu %lf, error %.2E\n", \
	exactuimg, polyuimg, polyerroruimg, sineuimg, sineerroruimg);
	printf("Norm: exact = %lf, polyu %lf, error %.2E, sineu %lf, error %.2E\n", \
	exactunorm, polyunorm, polyerrorunorm, sineunorm, sineerrorunorm);

	free(newx1);
	free(newx2);
	free(xminusy1);
	free(xminusy2);
	free(samplekernel);
	free(samplekernelimg);

	return 0;
*/

}


int BiCGSTAB(int size)
{
	int i, j, k, counter;
	FILE *fp3Dx, *fprealb;
	FILE *fp3DC;		//	DEBUG	//	

	double r[size], rhat[size], b[size], x[size], v[size], p[size], s[size], t[size];
	double error[size];
	double alpha = 1.0, rhoi = 1.0, rhoi1 = 1.0, omega = 1.0;
	double beta;
	double errornorm, preerrornorm, bnorm;
	double sum = 0.0;
	double accpower;
	double accuracy = DEFAULTACC;

	int countrun = 0;
	int checksize;
//	static int flag = 0;

	time_t loop_initime, loop_endtime, seconds;

	if ((fp3Dx = fopen("large3Dx.txt", "w+"))==NULL)
	{
		printf("Open file large3Dx.txt error.\n");
		return 0;
	}
//	fprintf(fp3Dx, "%d\n", size);

	if ((fprealb = fopen("large3Drealb.txt", "r"))==NULL)
	{
		printf("Open file large3Drealb.txt error.\n");
		return 0;
	}
	if (TEST_CONDITION==1)
	{
		if ((fp3DC = fopen("large3DC.dat", "w+"))==NULL)
		{
			printf("Open file large3DC.dat error.\n");
			return 0;
		}
	}
//	fprintf(fp3DC, "%d\n", size);

	printf("Enter the accuracy in 10^x ( 0 > x > -10)\n");
	scanf("%lf", &accpower);
	if ( (accpower > -10.0) && (accpower < 0.0) )
		accuracy = pow(10, accpower);
	printf("Using accuracy = %.8E.\n", accuracy);

	//	INITIALIZING THE THINGS NEEDED FOR BICONJUGATE GRADIENT STABILIZED	//

	fscanf(fprealb, "%d", &checksize);
	if (checksize != size)
	{
		printf("size error.\n");
		return 0;
	}
	for ( counter = 0; counter < size; counter++)
	{
		fscanf(fprealb, "%lf", &b[counter]);
		if (b[counter]!=b[counter])
			printf("b[%d] = %lf.\n", counter, b[counter]);
		r[counter] = b[counter];
		v[counter] = 0.0;
		p[counter] = 0.0;
		rhat[counter] = r[counter];
		if (rhat[counter] != r[counter])
			printf("counter = %d.\n, rhat = %lf, r = %lf.\n", counter, rhat[counter], r[counter]);
	}
	fclose(fprealb);
	///////////////////////////////////////////////////////////////////////////////////////
	bnorm = norm(size, b);

	if ( bnorm < pow(10,-13.0))	//	for zero BC	//
	{
		printf("Here.\n");
		r[0] = 0.000000001;
		rhat[0] = 0.000000001;
	}

	if ( bnorm < pow(10, -15.0) )
	{
		printf("Zero BC.\n");
		for ( i = 0; i < size; i++)
		{
			fprintf(fp3Dx, "%.15lf ", 0.0);
			
		}
		fclose(fp3Dx);
		return 0;
	}


	for ( i = 1; i <= 20000; i++)
	{
		loop_initime = time(NULL);

		rhoi1 = rhoi;					//	

		rhoi = inner(size, rhat, r);			//	rho = (r^, r)			//
//		if (rhoi!=rhoi)
//			printf("rhoi.\n");
		beta = (rhoi*alpha)/(rhoi1*omega);		//	beta = (rhoi/rhoi1) * (a/w)	//
		//		printf("beta = %lf\n", beta);
		for ( j = 0; j < size; j++)
		{
			p[j] = r[j] + beta * (p[j] - omega * v[j]);	//	p = r + beta( p - w * v )  //
//			if ( (p[j] != p[j])|| (p[j] + 1.0 ==p[j]) )
//				printf("j = %d, r[j] = %lf, v[j] = %lf, beta = %lf.\n", j, r[j], v[j], beta);
		}

		//////	FIRST MATRIX MULTIPLICATION	////////////
		for ( j = 0; j < size; j++)
		{
			sum = 0.0;
			for ( k = 0; k < size; k++)
			{
				if (i==1)
				{
					if (TEST_CONDITION==1)
						fprintf(fp3DC, "%.12lf ", get_C(size, j, k));
				}
				sum += (get_C(size, j, k) * p[k]);
			}
			if (i==1)
			{
				if (TEST_CONDITION==1)
					fprintf(fp3DC, "\n");
			}
			v[j] = sum;				//	v = A * p	//
			if ((v[j] != v[j])||(v[j] +1 == v[j]))
				printf("v wrong, j = %d, v[j] = %lf.\n", j, v[j]);
		}
		loop_endtime = time(NULL);
		seconds = loop_endtime - loop_initime;
		printf("Made to 1st. Took %ld minutes %ld seconds.\n", seconds/60, seconds % 60);
		/////	FIRST MATRIX MULTIPLICATION ENDS	////


		alpha = rhoi/inner(size, rhat, v);	//	a = rho/ (r^, v)	//
		for ( j = 0; j < size; j++)		
		{
			s[j] = r[j] - alpha * v[j];	//	s = r - a * v	//
			if ((s[j] !=s[j])||(s[j] + 1 ==s[j]))
				printf("s wrong, j = %d, alpha = %lf, r[j] = %lf", j, alpha, r[j]);
		}

		//////	SECOND MATRIX MULTIPLICATION	/////////////
		for ( j = 0; j < size; j++)
		{
			sum = 0.0;
			for ( k = 0; k < size; k++)
			{
				sum += (get_C(size, j, k) * s[k]);
			}
			t[j] = sum;				//	t = A * s		//
			if (t[j]!=t[j])
				printf("t wrong, j = %d.\n", j);
		}
		///////	SECOND MATRIX MULTIPLICATION ENDS	//////

		omega = inner(size, t, s)/inner(size, t, t);	//	w = (t, s)/(t, t)	//
		for ( j = 0; j < size; j++)
		{
			x[j] = x[j] + alpha * p[j] + omega * s[j];	//	x = x + a * p + w * s	//
		}
	
		///////	THIRD MATRIX MULTIPLICATION	/////////////
//		for (int j = 0; j < size; j++)
//		{
//			sum = 0.0;
//			for (int k = 0; k < size; k++)
//			{
//				sum += (get_C(size, j, k) * x[k]);
//			}
//			error[j] = sum - b[j];			//	err = b - A * x		//
//			if (error[j] != error[j])
//				printf("error wrong, j = %d, x = %lf, b = %lf, omega = %lf.", j, x[j], b[j], omega);
//		}

		for ( j = 0; j < size; j++)
		{
			r[j] = s[j] - omega * t[j];		//	r = s - w * t		//
		}

		preerrornorm = errornorm;
		errornorm = norm(size, r)/bnorm;
		printf("Run %d, error = %.8E\n", i, errornorm);

		loop_endtime = time(NULL);
		seconds = loop_endtime - loop_initime;
		printf("This run took %ld minutes and %ld seconds.\n", seconds/60, seconds % 60);
		countrun++;

		if (errornorm < accuracy)
			break;
		if (errornorm != errornorm)
			break;
		if (fabs(preerrornorm - errornorm) < pow(10, -3) * accuracy)
			break;
		

	}

	printf("After %d runs, the error is %.10E.\n", countrun, errornorm);
	
	for ( i = 0; i < size; i++)
		fprintf(fp3Dx, "%.15lf ", x[i]);
	
	fclose(fp3Dx);
//	fclose(fp3DC);
	return 0;
}


int BiCGcomplex(int size, double epsilon)
{
	int i, j, k, counter;
	FILE *fp3Dx, *fprealb, *fp3Dximg, *fpimgb;
//	FILE *fpverifyb, *fpverifybimg;
	FILE *fp3DC;		//	DEBUG	//	
	FILE *fp3DCimg;		//	DEBUG	//

	int compsize = 2 * size;
	double nonhomogeneous = 0.0;
	double nonhomogeneousimg = 0.0;
	double zx, zy;
	double zstarx[size], zstary[size];
//	double epsilon = FACTOR * h;
	int count = 0;;

	double tempC[2];
//	double r[size], rhat[size], x[size], v[size], p[size], s[size], t[size];
//	double rimg[size], rhatimg[size], ximg[size], vimg[size], pimg[size], simg[size], timg[size];
	double *x, *ximg, *r, *rimg, *rtilde, *rtildeimg;
	double *p, *pimg, *ptilde, *ptildeimg, *pold, *poldimg, *ptildeold, *ptildeoldimg;
	double *Akp, *Akpimg, *Akptilde, *Akptildeimg;

//	double b[3] = TESTB;
//	double bimg[3] = TESTBIMG;
	double *b, *bimg;




	double *error, *errorimg;
	double alpha = 1.0, rhoi = 1.0, rhoi1 = 1.0, omega = 1.0;
	double alphaimg = 0.0, rhoiimg = 0.0, rhoi1img = 0.0, omegaimg = 0.0;
	double beta, betaimg;
	double errornorm, preerrornorm;
	double comptemp, comptempnum, comptempnumimg, comptempden, comptempa, comptempb, comptempc, comptempd;
	double innerrold, innerroldimg;

	double sum = 0.0, regsum = 0.0, sumimg = 0.0, regsumimg = 0.0;
	double normb = 0.0;
	double BC_sum = 0.0;		//	for zero BC	//
	double accpower;
	double accuracy = DEFAULTACC;
	char denfilename[50], denimgfilename[50], filenameb[50], filenamebimg[50];

	int countrun = 0;
	int checksize, checksizeimg;

//	static int flag = 0;

	time_t loop_initime, loop_endtime, seconds;

	sprintf(denfilename, "%sDensityN%dk%d.txt", (ABS_MODE==1)?"Const_":"Rel_", (int)(N), (int)(WAVE));
	sprintf(denimgfilename, "%sDensityimgN%dk%d.txt", (ABS_MODE==1)?"Const_":"Rel_", (int)(N), (int)(WAVE));

	sprintf(filenameb, "%sbN%dk%d.txt", (ABS_MODE==1)?"Const_":"Rel_", N, (int)WAVE);
	sprintf(filenamebimg, "%sbimgN%dk%d.txt", (ABS_MODE==1)?"Const_":"Rel_", N, (int)WAVE);


	x = (double *) malloc(size * sizeof(double));
	b = (double *) malloc(size * sizeof(double));
	p = (double *) malloc(size * sizeof(double));
	pold = (double *) malloc(size * sizeof(double));
	r = (double *) malloc(size * sizeof(double));
	ptilde = (double *) malloc(size * sizeof(double));
	ptildeold = (double *) malloc(size * sizeof(double));
	rtilde = (double *) malloc(size * sizeof(double));
	Akp = (double *) malloc(size * sizeof(double));
	Akptilde = (double *) malloc(size * sizeof(double));
	error = (double *) malloc(size * sizeof(double));
	ximg = (double *) malloc(size * sizeof(double));
	bimg = (double *) malloc(size * sizeof(double));
	pimg = (double *) malloc(size * sizeof(double));
	poldimg = (double *) malloc(size * sizeof(double));
	rimg = (double *) malloc(size * sizeof(double));
	ptildeimg = (double *) malloc(size * sizeof(double));
	ptildeoldimg = (double *) malloc(size * sizeof(double));
	rtildeimg = (double *) malloc(size * sizeof(double));
	Akpimg = (double *) malloc(size * sizeof(double));
	Akptildeimg = (double *) malloc(size * sizeof(double));
	errorimg = (double *) malloc(size * sizeof(double));


	

	if ((fp3Dx = fopen(denfilename, "w+"))==NULL)
	{
		printf("Open file %s error.\n", denfilename);
		return 0;
	}
	if ((fp3Dximg = fopen(denimgfilename, "w+"))==NULL)
	{
		printf("Open file %s error.\n", denimgfilename);
		return 0;
	}
//	fprintf(fp3Dx, "%d\n", size);

//	if ((fpverifyb = fopen("verifyb.txt", "w+"))==NULL)
//	{
//		printf("Open file verifyb.txt error.\n");
//		exit(0);
//	}
//	if ((fpverifybimg = fopen("verifybimg.txt", "w+"))==NULL)
//	{
//		printf("Open file verifybimg.txt error.\n");
//		exit(0);
//	}

	if (TEST_CONDITION == 1)
	{
		if ((fp3DC = fopen("testC.txt", "w+"))==NULL)
		{
			printf("Open file testC.txt error.\n");
			return 0;
		}


//		fprintf(fp3DC, "%d\n", size);
		if ((fp3DCimg = fopen("testCimg.txt", "w+"))==NULL)
		{
			printf("Open file testCimg.txt error.\n");
			return 0;
		}
//		fprintf(fp3DCimg, "%d\n", size);
	}

	if (INTERACTIVE_ACCURACY == 1)
	{
		printf("Enter the accuracy in 10^x ( 0 > x > -10)\n");
		scanf("%lf", &accpower);
		if ( (accpower > -10.0) && (accpower < 0.0) )
			accuracy = pow(10, accpower);
	}
	printf("Using accuracy = %.8E.\n", accuracy);

	//	INITIALIZING THE THINGS NEEDED FOR BICONJUGATE GRADIENT STABILIZED	//


	if ((fprealb = fopen(filenameb, "r"))==NULL)
	{
		printf("Open file %s error.\n", filenameb);
		return 0;
	}
	if ((fpimgb = fopen(filenamebimg, "r"))==NULL)
	{
		printf("Open file %s error.\n", filenamebimg);
		return 0;
	}

	fscanf(fprealb, "%d", &checksize);
	fscanf(fpimgb, "%d", &checksizeimg);
	if ((checksize != size) || (checksizeimg!=size))
	{
		printf("size error. Realsize = %d, realpart = %d, imgpart = %d\n", size, checksize, checksizeimg);
		return 0;
	}
	for ( counter = 0; counter < size; counter++)
	{
		fscanf(fprealb, "%lf", &b[counter]);
		fscanf(fpimgb, "%lf", &bimg[counter]);

		if ((b[counter]!=b[counter])||(b[size+counter]!=b[size+counter]))
		{
			printf("b[%d] = %lf.\n", counter, b[counter]);
			printf("bimg[%d] = %lf.\n", counter, bimg[counter]);
			exit(0);
		}
	}
	fclose(fprealb);
	fclose(fpimgb);


////////////////////////////////////////////////////////////////////////////
//	TESTB;	
//	bimg[3] = TESTBIMG;

//	fprintf(fpverifyb, "%d\n", size);
//	fprintf(fpverifybimg, "%d\n", size);
//	for ( i = 0; i < size ; i++)
//	{
//		fprintf(fpverifyb, "%lf ", b[i]);
//		fprintf(fpverifybimg, "%lf ", bimg[i]);
//	}
////////////////////////////////////////////////////////////////////////////


	normb = sqrt(inner(size, b, b) + inner(size, bimg, bimg));
	printf("b norm is %lf\n", normb);

	for ( counter = 0; counter < size; counter++)
	{
		x[counter] = 0.0;
		ximg[counter] = 0.0;
		r[counter] = b[counter];
		rimg[counter] = bimg[counter];
		p[counter] = r[counter];
		pimg[counter] = rimg[counter];

		rtilde[counter] = r[counter];
		rtildeimg[counter] = (-1.0)*rimg[counter];
		ptilde[counter] = rtilde[counter];
		ptildeimg[counter] = rtildeimg[counter];

		if ( (rtilde[counter] != rtilde[counter]) || (rtildeimg[counter] != rtildeimg[counter]) )
		{
			printf("counter = %d.\n, rtilde = %lf, r = %lf.\n", counter, rtilde[counter], r[counter]);
			printf("rtildeimg = %lf, rimg = %lf.\n", rtildeimg[counter], rimg[counter]);
		}
	}
	//	(r, r^)	//
	innerrold = inner(size, r, rtilde) + inner(size, rimg, rtildeimg);
	innerroldimg = inner(size, rimg, rtilde) - inner(size, r, rtildeimg);

	///////////////////////////////////////////////////////////////////////////////////////


	for ( i = 0; i < size; i++)
	{
		BC_sum += fabs(r[i]);
		BC_sum += fabs(rimg[i]);
	}
	if ( BC_sum < pow(10,-13.0))	//	for zero BC	//
	{
		printf("Here.\n");
		r[0] = 0.000000001;
		rimg[0] = 0.000000001;
		rtilde[0] = 0.000000001;
		rtildeimg[0] = 0.000000001;
	}

	if ( BC_sum < pow(10, -15.0) )
	{
		printf("Zero BC.\n");
		for ( i = 0; i < size; i++)
		{
			fprintf(fp3Dx, "%.15lf ", 0.0);
			fprintf(fp3Dximg, "%.15lf ", 0.0);
		}
		fclose(fp3Dx);
		fclose(fp3Dximg);
		return 0;
	}


	for ( i = 1; i <= 2000; i++)
	{
		if (i == 1)
			loop_initime = time(NULL);

		//////	FIRST MATRIX MULTIPLICATION	////////////
		for ( j = 0; j < size; j++)
		{
			sum = 0.0;
			sumimg = 0.0;
			for ( k = 0; k < size; k++)
			{
//				printf("Here. i = %d, j = %d, k = %d\n", i, j, k);
				get_Ccomp(size, j, k, epsilon, tempC);
//				printf("After. i, j, k = (%d, %d, %d)\n", i, j, k);
				if (i==1)
				{
					if (TEST_CONDITION==1)
					{
						fprintf(fp3DC, "%.8lf ", tempC[0]);
						fprintf(fp3DCimg, "%.8lf ", tempC[1]);
					}
				}
				sum += (tempC[0] * p[k] - tempC[1] * pimg[k]);		//	A * p	//
				sumimg += (tempC[0] * pimg[k] + tempC[1] * p[k]);	//	A * p	//
			}
			if (i==1)
			{
				if (TEST_CONDITION==1)
				{
					fprintf(fp3DC, "\n");
					fprintf(fp3DCimg, "\n");
				}
			}
			Akp[j] = sum;				//	v = A * p	//
			Akpimg[j] = sumimg;

			if ((Akp[j] != Akp[j])||(Akp[j] +1 == Akp[j])||(Akpimg[j] != Akpimg[j])||(Akpimg[j]+1==Akpimg[j]))
			{
				printf("Akp wrong, j = %d, Akp[j] = %lf, Akpimg[j] = %lf.\n", j, Akp[j], Akpimg[j]);
				exit(0);
			}
//			printf("v[%d] = %.15lf\n", j, v[j]);
		}
		loop_endtime = time(NULL);
		seconds = loop_endtime - loop_initime;
		if (i == 1)
		{
			printf("Made to 1st. Took %ld minutes %ld seconds.\n", seconds/60, seconds % 60);
			if (TEST_CONDITION == 1)
			{
				fclose(fp3DC);
				fclose(fp3DCimg);
				printf("Printed matrix for condition number.\n");
				exit(0);
			}
		}
		/////	FIRST MATRIX MULTIPLICATION ENDS	////




		//////	SECOND MATRIX MULTIPLICATION	/////////////
		for ( j = 0; j < size; j++)
		{
			sum = 0.0;
			sumimg = 0.0;
			for ( k = 0; k < size; k++)
			{
				get_Ccomp(size, k, j, epsilon, tempC);
				sum += (tempC[0]*ptilde[k] + tempC[1]*ptildeimg[k]);
				sumimg += (tempC[0]*ptildeimg[k] - tempC[1]*ptilde[k]);
			}
			Akptilde[j] = sum;				//	Akptilde = A^H * p		//
			Akptildeimg[j] = sumimg;
			if ((Akptilde[j]!=Akptilde[j])||(Akptildeimg[j]!=Akptildeimg[j]))
			{
				printf("Akptilde wrong, j = %d, Akptilde[j] = %lf, Akptildeimg[j] = %lf.\n", j, Akptilde[j], Akptildeimg[j]);
				printf("alpha = %lf, r[%d] = %lf.\n", alpha, j, r[j]);
				exit(0);
			}
		}
		///////	SECOND MATRIX MULTIPLICATION ENDS	//////


		//	DEFINE IF A = a+bi, B = c+di, then (A,B) = sigma(A*B') = (a+bi)dot(c-di)	//
		//	THAT IS, CONJUGATE THE LATTER VECTOR	//

		comptempa = inner(size, r,rtilde) + inner(size, rimg, rtildeimg);
		comptempb = inner(size, rimg,rtilde) - inner(size, r, rtildeimg);
		comptempc = inner(size, Akp, ptilde) + inner(size, Akpimg, ptildeimg);
		comptempd = inner(size, Akpimg, ptilde) - inner(size, Akp, ptildeimg);
		comptempnum = comptempa*comptempc + comptempb*comptempd;
		comptempnumimg = comptempb*comptempc - comptempa*comptempd;
		comptempden = comptempc*comptempc + comptempd*comptempd;
		alpha = comptempnum/comptempden;
		alphaimg = comptempnumimg/comptempden;

		for ( j = 0; j < size; j++)
		{
			x[j] = x[j] + alpha * p[j] - alphaimg * pimg[j];
			ximg[j] = ximg[j] + alphaimg * p[j] + alpha * pimg[j];
			r[j] = r[j] - alpha * Akp[j] + alphaimg * Akpimg[j];
			rimg[j] = rimg[j] - alphaimg * Akp[j] - alpha * Akpimg[j];

			rtilde[j] = rtilde[j] - alpha * Akptilde[j] - alphaimg * Akptildeimg[j];
			rtildeimg[j] = rtildeimg[j] + alphaimg * Akptilde[j] - alpha * Akptildeimg[j];

			if (ximg[j] != ximg[j])
			{
				printf("j = %d, r = (%lf, %lf), p = (%lf, %lf), alpha = (%lf, %lf)\n", \
				j, r[j], rimg[j], p[j], pimg[j], alpha, alphaimg);
				exit(0);
			}

//			if ( (p[j] != p[j])|| (p[j] + 1.0 ==p[j]) )
//				printf("j = %d, r[j] = %lf, v[j] = %lf, beta = %lf.\n", j, r[j], v[j], beta);
		}


		//	BEGIN beta = (r^,r)_k+1/(r^, r)k		//
		comptempc = innerrold;
		comptempd = innerroldimg;	
		comptempa = inner(size, r, rtilde) + inner(size, rimg, rtildeimg);
		comptempb = inner(size, rimg, rtilde) - inner(size, r, rtildeimg);
		comptempnum = comptempa*comptempc + comptempb*comptempd;
		comptempnumimg = comptempb*comptempc - comptempa*comptempd;
		comptempden = comptempc*comptempc + comptempd*comptempd;
		beta = comptempnum/comptempden;
		betaimg = comptempnumimg/comptempden;

		innerrold = comptempa;
		innerroldimg = comptempb;

		//	END a = rho/ (r^, v)	//

		for ( j = 0; j < size; j++)		
		{
			pold[j] = p[j];
			poldimg[j] = pimg[j];
			ptildeold[j] = ptilde[j];
			ptildeoldimg[j] = ptildeimg[j];
		}

		for ( j = 0; j < size; j++)
		{
			p[j] = r[j] + beta * pold[j] - betaimg * poldimg[j];		//	p = r + beta * p	//
			pimg[j] = rimg[j] + betaimg * pold[j] + beta * poldimg[j];

			ptilde[j] = rtilde[j] + beta * ptildeold[j] + betaimg * ptildeoldimg[j];
			ptildeimg[j] = rtildeimg[j] - betaimg * ptildeold[j] + beta * ptildeoldimg[j];


			if ((p[j] !=p[j])||(p[j] + 1 ==p[j])||(pimg[j] !=pimg[j])||(pimg[j]+1==pimg[j]))
			{
				printf("s wrong, j = %d, alpha = %lf, r[j] = %lf, rimg[j] = %lf.\n", j, alpha, r[j], rimg[j]);
				exit(0);
			}
		}



		preerrornorm = errornorm;
		errornorm = sqrt(inner(size, r, r) + inner(size, rimg, rimg))/normb;
//		printf("Run %d, error = %.8E\n", i, errornorm);


//		loop_endtime = time(NULL);
//		seconds = loop_endtime - loop_initime;
//		printf("This run took %ld minutes and %ld seconds.\n", seconds/60, seconds % 60);
		countrun++;

		if (errornorm < accuracy)
			break;
		if (errornorm != errornorm)
			break;
		if (fabs(preerrornorm - errornorm) < pow(10, -3) * accuracy)
			break;
		
	}

	printf("After %d runs, the error is %.10E.\n", countrun, errornorm);
	
	for ( i = 0; i < size; i++)
	{
		fprintf(fp3Dx, "%.15lf ", x[i]);
		fprintf(fp3Dximg, "%.15lf ", ximg[i]);
		if ( (x[i] != x[i])|| (ximg[i] != ximg[i]) )
		{
			printf("AT i = %d, x = (%lf, %lf)\n", i, x[i], ximg[i]);
			exit(0);
		}

		Density[i] = x[i];
		Densityimg[i] = ximg[i];
	}
	get_Ccomp(size, -1, -1, epsilon, tempC);


	free(x);
	free(b);
	free(r);
	free(p);
	free(rtilde);
	free(ptilde);
	free(pold);
	free(ptildeold);
	free(Akp);
	free(Akptilde);
	free(error);
	free(ximg);
	free(bimg);
	free(rimg);
	free(pimg);
	free(rtildeimg);
	free(ptildeimg);
	free(poldimg);
	free(ptildeoldimg);
	free(Akpimg);
	free(Akptildeimg);
	free(errorimg);
	
	fclose(fp3Dx);
	fclose(fp3Dximg);
//	fclose(fp3DC);
//	fclose(fp3DCimg);
	return 0;
}


int parallel_BiCGcomplex(int size, double epsilon, int chunk)
{
	int i, j, k, counter;
	FILE *fp3Dx, *fprealb, *fp3Dximg, *fpimgb;
//	FILE *fpverifyb, *fpverifybimg;
	FILE *fp3DC;		//	DEBUG	//	
	FILE *fp3DCimg;		//	DEBUG	//

	int compsize = 2 * size;
	double nonhomogeneous = 0.0;
	double nonhomogeneousimg = 0.0;
	double zx, zy;
//	double zstarx[size], zstary[size];
	int count = 0;;

	double tempA[2];
	double *x, *ximg, *r, *rimg, *rtilde, *rtildeimg;
	double *p, *pimg, *ptilde, *ptildeimg, *pold, *poldimg, *ptildeold, *ptildeoldimg;
	double *Akp, *Akpimg, *Akptilde, *Akptildeimg;
	double **matrix, **matriximg;

	double *b, *bimg;




	double *error, *errorimg;
	double alpha = 1.0, rhoi = 1.0, rhoi1 = 1.0, omega = 1.0;
	double alphaimg = 0.0, rhoiimg = 0.0, rhoi1img = 0.0, omegaimg = 0.0;
	double beta, betaimg;
	double errornorm, preerrornorm;
	double comptemp, comptempnum, comptempnumimg, comptempden, comptempa, comptempb, comptempc, comptempd;
	double innerrold, innerroldimg;

	double sum = 0.0, regsum = 0.0, sumimg = 0.0, regsumimg = 0.0;
	double normb = 0.0;
	double BC_sum = 0.0;		//	for zero BC	//
	double accpower;
	double accuracy = DEFAULTACC;
	char denfilename[50], denimgfilename[50], filenameb[50], filenamebimg[50];

	int countrun = 0;
	int checksize, checksizeimg;


	time_t loop_initime, loop_endtime, seconds;

	sprintf(denfilename, "%sDensityN%dk%d.txt", (ABS_MODE==1)?"Const_":"Rel_", (int)(N), (int)(WAVE));
	sprintf(denimgfilename, "%sDensityimgN%dk%d.txt", (ABS_MODE==1)?"Const_":"Rel_", (int)(N), (int)(WAVE));

	sprintf(filenameb, "%sbN%dk%d.txt", (ABS_MODE==1)?"Const_":"Rel_", (int)(N), (int)(WAVE));
	sprintf(filenamebimg, "%sbimgN%dk%d.txt", (ABS_MODE==1)?"Const_":"Rel_", (int)(N), (int)(WAVE));


	x = (double *) malloc(size * sizeof(double));
	b = (double *) malloc(size * sizeof(double));
	p = (double *) malloc(size * sizeof(double));
	pold = (double *) malloc(size * sizeof(double));
	r = (double *) malloc(size * sizeof(double));
	ptilde = (double *) malloc(size * sizeof(double));
	ptildeold = (double *) malloc(size * sizeof(double));
	rtilde = (double *) malloc(size * sizeof(double));
	Akp = (double *) malloc(size * sizeof(double));
	Akptilde = (double *) malloc(size * sizeof(double));
	error = (double *) malloc(size * sizeof(double));
	ximg = (double *) malloc(size * sizeof(double));
	bimg = (double *) malloc(size * sizeof(double));
	pimg = (double *) malloc(size * sizeof(double));
	poldimg = (double *) malloc(size * sizeof(double));
	rimg = (double *) malloc(size * sizeof(double));
	ptildeimg = (double *) malloc(size * sizeof(double));
	ptildeoldimg = (double *) malloc(size * sizeof(double));
	rtildeimg = (double *) malloc(size * sizeof(double));
	Akpimg = (double *) malloc(size * sizeof(double));
	Akptildeimg = (double *) malloc(size * sizeof(double));
	errorimg = (double *) malloc(size * sizeof(double));

	matrix = (double **) malloc(size * sizeof(double *));
	matriximg = (double **) malloc(size * sizeof(double *));
	for (i = 0; i < size; i++)
	{
		matrix[i] = (double * ) malloc(size * sizeof(double));
		matriximg[i] = (double * ) malloc(size * sizeof(double));
	}

	

	if ((fp3Dx = fopen(denfilename, "w"))==NULL)
	{
		printf("Open file %s error.\n", denfilename);
		return 0;
	}
	if ((fp3Dximg = fopen(denimgfilename, "w"))==NULL)
	{
		printf("Open file %s error.\n", denimgfilename);
		return 0;
	}
//	fprintf(fp3Dx, "%d\n", size);

//	if ((fpverifyb = fopen("verifyb.txt", "w+"))==NULL)
//	{
//		printf("Open file verifyb.txt error.\n");
//		exit(0);
//	}
//	if ((fpverifybimg = fopen("verifybimg.txt", "w+"))==NULL)
//	{
//		printf("Open file verifybimg.txt error.\n");
//		exit(0);
//	}




	if (TEST_CONDITION == 1)
	{
		if ((fp3DC = fopen("testC.txt", "w+"))==NULL)
		{
			printf("Open file testC.txt error.\n");
			return 0;
		}


//		fprintf(fp3DC, "%d\n", size);
		if ((fp3DCimg = fopen("testCimg.txt", "w+"))==NULL)
		{
			printf("Open file testCimg.txt error.\n");
			return 0;
		}
//		fprintf(fp3DCimg, "%d\n", size);
	}

	if (INTERACTIVE_ACCURACY == 1)
	{
		printf("Enter the accuracy in 10^x ( 0 > x > -10)\n");
		scanf("%lf", &accpower);
		if ( (accpower > -10.0) && (accpower < 0.0) )
			accuracy = pow(10, accpower);
	}
	printf("Using accuracy = %.8E.\n", accuracy);

	//	INITIALIZING THE THINGS NEEDED FOR BICONJUGATE GRADIENT STABILIZED	//


	if ((fprealb = fopen(filenameb, "r"))==NULL)
	{
		printf("Open file %s error.\n", filenameb);
		return 0;
	}
	if ((fpimgb = fopen(filenamebimg, "r"))==NULL)
	{
		printf("Open file %s error.\n", filenamebimg);
		return 0;
	}

	fscanf(fprealb, "%d", &checksize);
	fscanf(fpimgb, "%d", &checksizeimg);
	if ((checksize != size) || (checksizeimg!=size))
	{
		printf("size error. Realsize = %d, realpart = %d, imgpart = %d\n", size, checksize, checksizeimg);
		return 0;
	}

	for ( counter = 0; counter < size; counter++)
	{
		fscanf(fprealb, "%lf", &b[counter]);
		fscanf(fpimgb, "%lf", &bimg[counter]);

		if ((b[counter]!=b[counter])||(b[size+counter]!=b[size+counter]))
		{
			printf("b[%d] = %lf.\n", counter, b[counter]);
			printf("bimg[%d] = %lf.\n", counter, bimg[counter]);
			exit(0);
		}
	}
	fclose(fprealb);
	fclose(fpimgb);



////////////////////////////////////////////////////////////////////////////

//	fprintf(fpverifyb, "%d\n", size);
//	fprintf(fpverifybimg, "%d\n", size);
//	for ( i = 0; i < size ; i++)
//	{
//		fprintf(fpverifyb, "%lf ", b[i]);
//		fprintf(fpverifybimg, "%lf ", bimg[i]);
//	}

////////////////////////////////////////////////////////////////////////////


	normb = sqrt(omp_inner(size, b, b, chunk) + omp_inner(size, bimg, bimg, chunk ));
	printf("b norm is %lf\n", normb);

	#pragma omp parallel for private(counter) schedule(static, chunk)
	for ( counter = 0; counter < size; counter++)
	{
//		if (counter== 0)
//		{
//			printf("%d threads in parallel BiCG.\n", omp_get_num_threads());
//		}
		x[counter] = 0.0;
		ximg[counter] = 0.0;
		r[counter] = b[counter];
		rimg[counter] = bimg[counter];
		p[counter] = b[counter];
		pimg[counter] = bimg[counter];

		rtilde[counter] = b[counter];
		rtildeimg[counter] = (-1.0)*bimg[counter];
		ptilde[counter] = b[counter];
		ptildeimg[counter] = (-1.0)*bimg[counter];

//		if ( (rtilde[counter] != rtilde[counter]) || (rtildeimg[counter] != rtildeimg[counter]) )
//		{
//			printf("counter = %d.\n, rtilde = %lf, r = %lf.\n", counter, rtilde[counter], r[counter]);
//			printf("rtildeimg = %lf, rimg = %lf.\n", rtildeimg[counter], rimg[counter]);
//		}
	}
	//	(r, r^)	//
	innerrold = omp_inner(size, r, rtilde, chunk) + omp_inner(size, rimg, rtildeimg, chunk);
	innerroldimg = omp_inner(size, rimg, rtilde, chunk) - omp_inner(size, r, rtildeimg, chunk);

	///////////////////////////////////////////////////////////////////////////////////////

	#pragma omp parallel for private(i) schedule(static, chunk) reduction(+:BC_sum)
	for ( i = 0; i < size; i++)
	{
		BC_sum += fabs(r[i]) + fabs(rimg[i]);
	}
	if ( BC_sum < pow(10,-13.0))	//	for zero BC	//
	{
		printf("Here.\n");
		r[0] = 0.000000001;
		rimg[0] = 0.000000001;
		rtilde[0] = 0.000000001;
		rtildeimg[0] = 0.000000001;
	}

	if ( BC_sum < pow(10, -15.0) )
	{
		printf("Zero BC.\n");
		for ( i = 0; i < size; i++)
		{
			fprintf(fp3Dx, "%.15lf ", 0.0);
			fprintf(fp3Dximg, "%.15lf ", 0.0);
		}
		fclose(fp3Dx);
		fclose(fp3Dximg);
		return 0;
	}


	for ( i = 1; i <= 2000; i++)
	{
		if (i == 1)
			loop_initime = time(NULL);

		#pragma omp parallel for private(j) schedule(static, chunk)
		for (j = 0; j < size; j++)
		{
			Akp[j] = 0.0;
			Akpimg[j] = 0.0;
			Akptilde[j] = 0.0;
			Akptildeimg[j] = 0.0;
		}

		//////	FIRST FIRST MATRIX MULTIPLICATION	////////////
		if (i == 1)
		{
			for ( j = 0; j < size; j++)
			{
			for ( k = 0; k < size; k++)
			{
				omp_get_Ccomp(size, j, k, epsilon, tempA);
				matrix[j][k] = tempA[0];
				matriximg[j][k] = tempA[1];
				if (TEST_CONDITION==1)
				{
					fprintf(fp3DC, "%.8lf ", tempA[0]);
					fprintf(fp3DCimg, "%.8lf ", tempA[1]);
				}
				Akp[j] += (tempA[0] * p[k] - tempA[1] * pimg[k]);		//	A * p	//
				Akpimg[j] += (tempA[0] * pimg[k] + tempA[1] * p[k]);	//	A * p	//
			}	//	k

				if (TEST_CONDITION==1)
				{
					fprintf(fp3DC, "\n");
					fprintf(fp3DCimg, "\n");
				}
			}	//	j
			loop_endtime = time(NULL);
			seconds = loop_endtime - loop_initime;
			printf("Made to 1st. Took %ld minutes %ld seconds.\n", seconds/60, seconds % 60);
			if (TEST_CONDITION == 1)
			{
				fclose(fp3DC);
				fclose(fp3DCimg);
				printf("Printed matrix for condition number.\n");
				exit(0);
			}
			omp_get_Ccomp(size, -1, -1, epsilon, tempA);
		}
		else
		{
			///////	FIRST FIRST MATRIX MULTIPLICATION ENDS	//////////////
			#pragma omp parallel for private(j,k) schedule(static, chunk)
			for ( j = 0; j < size; j++)
			{
//				double tempC[2];
			
				for ( k = 0; k < size; k++)
				{
//					get_Ccomp(size, j, k, epsilon, tempC);
//				if (TEST_CONDITION==1)
//				{
//					fprintf(fp3DC, "%.8lf ", tempC[0]);
//					fprintf(fp3DCimg, "%.8lf ", tempC[1]);
//				}
					Akp[j] += (matrix[j][k] * p[k] - matriximg[j][k] * pimg[k]);		//	A * p	//
					Akpimg[j] += (matrix[j][k] * pimg[k] + matriximg[j][k] * p[k]);	//	A * p	//
				}
//			if ((Akp[j] != Akp[j])||(Akp[j] +1 == Akp[j])||(Akpimg[j] != Akpimg[j])||(Akpimg[j]+1==Akpimg[j]))
//			{
//				printf("Akp wrong, j = %d, Akp[j] = %lf, Akpimg[j] = %lf.\n", j, Akp[j], Akpimg[j]);
//				exit(0);
//			}
//			printf("v[%d] = %.15lf\n", j, v[j]);
			}				
		}
		/////	FIRST MATRIX MULTIPLICATION ENDS	////





		//////	SECOND MATRIX MULTIPLICATION	/////////////
		#pragma omp parallel for private(j,k) schedule(static, chunk)
		for ( j = 0; j < size; j++)
		{
//			double tempC[2];
			for ( k = 0; k < size; k++)
			{
//				get_Ccomp(size, k, j, epsilon, tempC);
				Akptilde[j] += (matrix[k][j]*ptilde[k] + matriximg[k][j]*ptildeimg[k]);
				Akptildeimg[j] += (matrix[k][j]*ptildeimg[k] - matriximg[k][j]*ptilde[k]);
			}
//			if ((Akptilde[j]!=Akptilde[j])||(Akptildeimg[j]!=Akptildeimg[j]))
//			{
//				printf("Akptilde wrong, j = %d, Akptilde[j] = %lf, Akptildeimg[j] = %lf.\n", j, Akptilde[j], Akptildeimg[j]);
//				printf("alpha = %lf, r[%d] = %lf.\n", alpha, j, r[j]);
//				exit(0);
//			}
		}
//	printf("here.\n");
//	exit(0);
		///////	SECOND MATRIX MULTIPLICATION ENDS	//////


		//	DEFINE IF A = a+bi, B = c+di, then (A,B) = sigma(A*B') = (a+bi)dot(c-di)	//
		//	THAT IS, CONJUGATE THE LATTER VECTOR	//


		comptempa = omp_inner(size, r,rtilde, chunk) + omp_inner(size, rimg, rtildeimg, chunk);
		comptempb = omp_inner(size, rimg,rtilde, chunk) - omp_inner(size, r, rtildeimg, chunk);
		comptempc = omp_inner(size, Akp, ptilde, chunk) + omp_inner(size, Akpimg, ptildeimg, chunk);
		comptempd = omp_inner(size, Akpimg, ptilde, chunk) - omp_inner(size, Akp, ptildeimg, chunk);
		comptempnum = comptempa*comptempc + comptempb*comptempd;
		comptempnumimg = comptempb*comptempc - comptempa*comptempd;
		comptempden = comptempc*comptempc + comptempd*comptempd;
		alpha = comptempnum/comptempden;
		alphaimg = comptempnumimg/comptempden;


		#pragma omp parallel sections private(j)
		{
			#pragma omp section
			{
				for ( j = 0; j < size; j++)
					x[j] = x[j] + alpha * p[j] - alphaimg * pimg[j];
			}
			#pragma omp section
			{
				for ( j = 0; j < size; j++)
					ximg[j] = ximg[j] + alphaimg * p[j] + alpha * pimg[j];
			}
			#pragma omp section
			{
				for ( j = 0; j < size; j++)
					r[j] = r[j] - alpha * Akp[j] + alphaimg * Akpimg[j];
			}
			#pragma omp section
			{
				for ( j = 0; j < size; j++)
					rimg[j] = rimg[j] - alphaimg * Akp[j] - alpha * Akpimg[j];
			}
			#pragma omp section
			{
				for ( j = 0; j < size; j++)
					rtilde[j] = rtilde[j] - alpha * Akptilde[j] - alphaimg * Akptildeimg[j];
			}
			#pragma omp section
			{
				for ( j = 0; j < size; j++)
					rtildeimg[j] = rtildeimg[j] + alphaimg * Akptilde[j] - alpha * Akptildeimg[j];
			}
//			if (ximg[j] != ximg[j])
//			{
//				printf("j = %d, r = (%lf, %lf), p = (%lf, %lf), alpha = (%lf, %lf)\n", \
//				j, r[j], rimg[j], p[j], pimg[j], alpha, alphaimg);
//				exit(0);
//			}

//			if ( (p[j] != p[j])|| (p[j] + 1.0 ==p[j]) )
//				printf("j = %d, r[j] = %lf, v[j] = %lf, beta = %lf.\n", j, r[j], v[j], beta);
		}


		//	BEGIN beta = (r^,r)_k+1/(r^, r)k		//
		comptempc = innerrold;
		comptempd = innerroldimg;	
		comptempa = omp_inner(size, r, rtilde, chunk) + omp_inner(size, rimg, rtildeimg, chunk);
		comptempb = omp_inner(size, rimg, rtilde, chunk) - omp_inner(size, r, rtildeimg, chunk);
		comptempnum = comptempa*comptempc + comptempb*comptempd;
		comptempnumimg = comptempb*comptempc - comptempa*comptempd;
		comptempden = comptempc*comptempc + comptempd*comptempd;
		beta = comptempnum/comptempden;
		betaimg = comptempnumimg/comptempden;

		innerrold = comptempa;
		innerroldimg = comptempb;

		//	END a = rho/ (r^, v)	//
		#pragma omp parallel for private(j) schedule(static, chunk)
		for ( j = 0; j < size; j++)		
		{
			pold[j] = p[j];
			poldimg[j] = pimg[j];
			ptildeold[j] = ptilde[j];
			ptildeoldimg[j] = ptildeimg[j];
		}

		#pragma omp parallel for private(j) schedule(static, chunk)
		for ( j = 0; j < size; j++)
		{
			p[j] = r[j] + beta * pold[j] - betaimg * poldimg[j];		//	p = r + beta * p	//
			pimg[j] = rimg[j] + betaimg * pold[j] + beta * poldimg[j];

			ptilde[j] = rtilde[j] + beta * ptildeold[j] + betaimg * ptildeoldimg[j];
			ptildeimg[j] = rtildeimg[j] - betaimg * ptildeold[j] + beta * ptildeoldimg[j];
//			if ((p[j] !=p[j])||(p[j] + 1 ==p[j])||(pimg[j] !=pimg[j])||(pimg[j]+1==pimg[j]))
//			{
//				printf("s wrong, j = %d, alpha = %lf, r[j] = %lf, rimg[j] = %lf.\n",\
//					 j, alpha, r[j], rimg[j]);
//				exit(0);
//			}
		}




		preerrornorm = errornorm;
		errornorm = sqrt(omp_inner(size, r, r, chunk) + omp_inner(size, rimg, rimg, chunk))/normb;
//		printf("Run %d, error = %.8E\n", i, errornorm);


//		loop_endtime = time(NULL);
//		seconds = loop_endtime - loop_initime;
//		printf("This run took %ld minutes and %ld seconds.\n", seconds/60, seconds % 60);
		countrun++;

		if (errornorm < accuracy)
			break;
		if (errornorm != errornorm)
			break;
		if (fabs(preerrornorm - errornorm) < pow(10, -3) * accuracy)
			break;
		
	}

	printf("After %d runs, the error is %.10E.\n", countrun, errornorm);
	
	#pragma omp parallel for private(i) schedule(static, chunk)
	for (i = 0; i < size; i++)
	{
		Density[i] = x[i];
		Densityimg[i] = ximg[i];
	}


	#pragma omp parallel sections private(i)
	{
		#pragma omp section
		{
			for ( i = 0; i < size; i++)
			{
				fprintf(fp3Dx, "%.15lf ", x[i]);
				if (x[i] != x[i])
				{
					printf("AT i = %d, x = %lf\n", i, x[i]);
					exit(0);
				}
			}
		}
		#pragma omp section
		{
			for ( i = 0; i < size; i++)
			{
				fprintf(fp3Dximg, "%.15lf ", ximg[i]);
				if (ximg[i] != ximg[i])
				{
					printf("AT i = %d, ximg = %lf\n", i, ximg[i]);
					exit(0);
				}
			}
		}


	}
//	get_Ccomp(size, -1, -1, epsilon, tempA);



	for ( i = 0; i < size; i++)
	{
		free(matrix[i]);
		free(matriximg[i]);
	}
	free(matrix);
	free(matriximg);
	free(x);
	free(b);
	free(r);
	free(p);
	free(rtilde);
	free(ptilde);
	free(pold);
	free(ptildeold);
	free(Akp);
	free(Akptilde);
	free(error);
	free(ximg);
	free(bimg);
	free(rimg);
	free(pimg);
	free(rtildeimg);
	free(ptildeimg);
	free(poldimg);
	free(ptildeoldimg);
	free(Akpimg);
	free(Akptildeimg);
	free(errorimg);
	
	fclose(fp3Dx);
	fclose(fp3Dximg);
//	fclose(fp3DC);
//	fclose(fp3DCimg);
	return 0;
}




double get_C(int size, int pi, int pj)
{
	int i, j, k;
	int counti = 0;
//	int countj = 0;
	int counter = 0;
	static double *zx, *zy, *zz;
	static double *zstarx, *zstary, *zstarz;
	static double *J, *delta;
	static double *gradx, *grady, *gradz;
	static double *result;
	static double regfactor;
	static double tau, threshold;
	static double epsilon;
	static double term1, term2, term3, term4;
	static int flag = 0;
	static int ini = 0;
	int indicator;


	static double *sumpartial, *sumpartialreg;
	double tempsum = 0.0, tempsumreg = 0.0;
//	static int fillpartialcount = 0, fillpartial = 0;

	static double **matrix;


	double tempzx, tempzy, tempzz, tstarx, tstary, tstarz;
	double starxi, staryi, starzi, starxj, staryj, starzj;
	double thisdist, thisJ;
	double C, C_reg;
	double partialx, partialy, partialz, partial, partial_reg;

	static unsigned long int fastflag = 0;
	static unsigned long int fillinflag = 0;
//	double h = (INT_U - INT_L)/N;


	if (size <= FASTSIZEBOUND)
	{
		if ( fastflag >= size*size )
			return matrix[pi][pj];
		fastflag++;
	}



	if (ini == 0)	//	FIRST TIME, THE ONLY TIME DOING INITIALIZATION FOR X, Y, Z POINTS //
	{
		tau = H * TAU_FACTOR;
		threshold = tau * tau;
		epsilon = FACTOR * H;

		term1 = 1.0/(4.0*PI*tau*radius);		//	1/8pi*tau(2/r)	//
		term2 = (-5.0*tau)/(256.0*pow(radius,3.0)*PI);	//	-(1/pi)*(5/512)*(2/r^3)*tau	//
		term3 = (-50.0*tau)/(1536.0*pow(radius,3.0)*PI);//	-(1/pi)*(25/1536)*(2/r^3)*tau	//
		

		if (INTERIOR_MODE==1)
			regfactor = term1 + term2 + term3;
		else
			regfactor = -1.0 * (term1 + term2 + term3);

		zx = (double *) malloc(size * sizeof(double));
		zy = (double *) malloc(size * sizeof(double));
		zz = (double *) malloc(size * sizeof(double));
		zstarx = (double *) malloc(size * sizeof(double));
		zstary = (double *) malloc(size * sizeof(double));
		zstarz = (double *) malloc(size * sizeof(double));
		delta = (double *) malloc(size * sizeof(double));
		J = (double *) malloc(size * sizeof(double));
		gradx = (double *) malloc(size * sizeof(double));
		grady = (double *) malloc(size * sizeof(double));
		gradz = (double *) malloc(size * sizeof(double));
		result= (double *) malloc(size * sizeof(double));

		sumpartial = (double *) malloc(size * sizeof(double));
		sumpartialreg = (double *) malloc(size * sizeof(double));
		

		if ( (zx==NULL)||(zy==NULL)||(zz==NULL)|| (zstarx==NULL) || (zstary==NULL) || (zstarz==NULL) )
			printf("Memory allocation error....points.\n");
		if ( (delta==NULL)||(J==NULL)||(gradx==NULL)|| (grady==NULL) || (gradz==NULL) )
			printf("Memory allocation error....calcs.\n");
		if (result==NULL)
			printf("Memory allocation error....result.\n");

		for ( i = 0; i < SIZE; i++)
		{
			tempzy = INT_L + H * (double)(i);
			for ( j = 0; j < SIZE; j++)
			{
				tempzx = INT_L + H * (double)(j);
				for( k = 0; k < SIZE; k++)
				{
					tempzz = INT_L + H * (double) (k);
					thisdist = DIST[i][j][k];
					//	if not in tubular region, don't calculate.	//
					if (fabs(thisdist) > epsilon)
						continue;

					zx[counti] = tempzx;
					zy[counti] = tempzy;
					zz[counti] = tempzz;
//					tstarx = starx(tempzx, tempzy, tempzz);
//					tstary = stary(tempzx, tempzy, tempzz);
//					tstarz = starz(tempzx, tempzy, tempzz);

					tstarx = tempzx - thisdist * DIST_DX[i][j][k];
					tstary = tempzy - thisdist * DIST_DY[i][j][k];
					tstarz = tempzz - thisdist * DIST_DZ[i][j][k];

					zstarx[counti] = tstarx;
					zstary[counti] = tstary;
					zstarz[counti] = tstarz;

					delta[counti] = cal_delta(epsilon, thisdist);
//					J[counti] = cal_J(J_M, tempzx, tempzy, tempzz);
					J[counti] = 1. + 2.*thisdist * MCURV[i][j][k] + thisdist * thisdist * GCURV[i][j][k];
//					gradx[counti] = cal_interpgrad(1, tstarx, tstary, tstarz);
//					grady[counti] = cal_interpgrad(2, tstarx, tstary, tstarz);
//					gradz[counti] = cal_interpgrad(3, tstarx, tstary, tstarz);

					indicator = cal_interpgrad(N, H, tstarx, tstary, tstarz, DIST_DX, DIST_DY, DIST_DZ, \
						&gradx[counti], &grady[counti], &gradz[counti]);

					if (indicator == 1)
					{
						printf("At (%lf, %lf, %lf), grad (%lf, %lf, %lf), dist %lf\n",
						tempzx, tempzy, tempzz, DIST_DX[i][j][k], DIST_DY[i][j][k], \
						DIST_DZ[i][j][k], thisdist);
						printf("Projected to (%lf, %lf, %lf), dist %lf\n", \
						tstarx, tstary, tstarz, cal_dist(DIST, tstarx, tstary, tstarz, N, H));
					}

					result[counti] = HCUBED * delta[counti] * J[counti];

					counti++;	// in the end should total up to size	//
				}
			}
		}
		if (size <= FASTSIZEBOUND)
		{
			matrix = (double **) malloc(size * sizeof(double *));
			for ( i = 0; i < size; i++)
			{
				matrix[i] = (double *) malloc(size * sizeof(double));
				if (matrix[i] == NULL)
				{
					printf("Memory allocation error for matrix[%d].\n", i);
					exit(0);
				}
			}
		}
		else
		{
			matrix = (double **) malloc(FASTSIZEBOUND * sizeof(double *));
			for ( i = 0; i < FASTSIZEBOUND; i++)
			{
				matrix[i] = (double *) malloc(FASTSIZEBOUND *sizeof(double));
				if (matrix[i] == NULL)
				{
					printf("Memory allocation error for matrix[%d].\n", i);
					exit(0);
				}
			}
		}
		ini = 1;	//	DENOTING DATA INITIALIZED	//
	}
//	if (counti==size)
//		printf("It's a pass.\n");

	starxi = zstarx[pi];
	staryi = zstary[pi];
	starzi = zstarz[pi];
	starxj = zstarx[pj];
	staryj = zstary[pj];
	starzj = zstarz[pj];

	partialx = cal_partial(1, starxi, staryi, starzi, starxj, staryj, starzj);// x coordinate	//
	partialy = cal_partial(2, starxi, staryi, starzi, starxj, staryj, starzj);// y coordinate	//
	partialz = cal_partial(3, starxi, staryi, starzi, starxj, staryj, starzj);// z coordinate	//

	if (INTERIOR_MODE==1)
		partial = (-1.0) * ( partialx * gradx[pj] + partialy * grady[pj] + partialz * gradz[pj] );
	else
		partial = partialx * gradx[pj] + partialy * grady[pj] + partialz * gradz[pj];

	
	if ( ( (starxi-starxj)*(starxi-starxj) + (staryi-staryj)*(staryi-staryj) + (starzi-starzj)*(starzi-starzj) ) < threshold )
	{
		partial_reg = regfactor;
	}
	else
	{
		partial_reg = partial;
	}
					
	//	calculate matrix elements	//

	C = result[pj] * partial ;		//	( h^3*deltaj*Jj ) * partial	//
	C_reg = result[pj] * partial_reg;


/*	if (ini==1)
	{
		for (int i = 0; i < size; i++)
		{
			starxi = zstarx[i];
			staryi = zstary[i];
			starzi = zstarz[i];
			tempsum = 0.0;
			tempsumreg = 0.0;
			for (int j = 0; j < size; j++)
			{

				starxj = zstarx[j];
				staryj = zstary[j];
				starzj = zstarz[j];

				partialx = cal_partial(1, starxi, staryi, starzi, starxj, staryj, starzj);// x coordinate	//
				partialy = cal_partial(2, starxi, staryi, starzi, starxj, staryj, starzj);// y coordinate	//
				partialz = cal_partial(3, starxi, staryi, starzi, starxj, staryj, starzj);// z coordinate	//

				partial = (-1.0) * ( partialx * gradx[j] + partialy * grady[j] + partialz * gradz[j] );

	
				if ( ( (starxi-starxj)*(starxi-starxj) + (staryi-staryj)*(staryi-staryj) + (starzi-starzj)*(starzi-starzj) ) < threshold )
				{
					partial_reg = regfactor;
				}
				else
				{
					partial_reg = partial;
				}
					
				//	calculate matrix elements	//

				tempsum += result[j] * partial ;		//	( h^3*deltaj*Jj ) * partial	//
				tempsumreg += result[j] * partial_reg;
			}
			sumpartial[i] = tempsum;
			sumpartialreg[i] = tempsumreg;
		}
		ini++;
	}
*/


	//	C = B + 1/2 I	//
	if (pi == pj)
	{
		C += 0.5;	//(delta*cal_J(J_M, wx, wy));
		C_reg += 0.5;	//(delta*cal_J(J_M, wx, wy));

//		C += 1.0;	//(delta*cal_J(J_M, wx, wy));
//		C_reg += 1.0;	//(delta*cal_J(J_M, wx, wy));
//		C -= sumpartial[pi];
//		C_reg -= sumpartialreg[pi];
	}

	if (size <= FASTSIZEBOUND)
	{
		if (REG_MODE==0)
			matrix[pi][pj] = C;
		else
			matrix[pi][pj] = C_reg;
	}
	//	TRY TO SAVE MORE	//
	else if (fillinflag < (long int) (FASTSIZEBOUND) * (long int) (FASTSIZEBOUND) )
	{
		if ( (pi < FASTSIZEBOUND) && (pj < FASTSIZEBOUND) )
		{
			if (REG_MODE != 1)
				matrix[pi][pj] = C;
			else
				matrix[pi][pj] = C_reg;
			fillinflag++;
		}
	}

	if ((C_reg!=C_reg)||(C_reg+1==C_reg))
	{
		printf("i = %d, j = %d, partial = %.10lf, C_reg = %.10lf.\n", pi,  pj, partial_reg, C_reg);
	}
	if ((C!=C)||(C+1==C))
	{
		printf("i = %d, j = %d, partial = %.10lf, C = %.10lf.\n", pi, pj, partial, C);
	}

	if (REG_MODE != 1)
		return C;
	else
		return C_reg;
}


double get_Ccomp(int size, int pi, int pj, double epsilon, double tempC[])
{
	FILE *fp, *fp2;

	int i, j, k, indicator;
	int counti = 0;
	int countj = 0;

	static int errorflag = 0;
	static double *zx, *zy, *zz, *wx, *wy, *wz;
	static double *zstarx, *zstary, *zstarz, *wstarx, *wstary, *wstarz;
	static double *J, *delta;
	static double *gradx, *grady, *gradz;
	static double *result, *regfactor, *regfactorimg;
//	static double epsilon;
	static double tau, threshold;
	static double term1, term2, term3, term4, thisdist, cur_meancurv, cur_gcurv;
	double regterm1, regterm2, regterm3, regterm4, gradnorm;

	static double **matrix;
	static double **matriximg;
	static double *temp;

	static const double recipfourpi = 1./(4.*PI);
	double xmy1,xmy2, xmy3, xmynorm, recipxmynorm;
	double SineValue, CosValue, extSineValue, extCosValue;
	double commonterm1, commonterm2, parameter, denom, inner_pd;


	double tempzx, tempzy, tempzz, tstarx, tstary, tstarz;
	double starxi, staryi, starzi, starxj, staryj, starzj;
	double C, C_reg;
	double partial, partial_reg, partialimg, partialimg_reg;
	double partialx, partialy, partialz, partialximg, partialyimg, partialzimg;
//	double h = (INT_U - INT_L)/N;
//	double epsilon = FACTOR * h;
//	double tau = h;
//	time_t begin, middle, end, seconds;

	static int ini = 0;
	static int flag = 0;
	static long int fastflag = 0;
	static long int fillinflag = 0;
	static int regcount = 0;
	static int printregcount = 0;



	if ( (pi == -1) && (pj == -1) )
	{
		if (size <= FASTSIZEBOUND)
		{
			for ( i = 0; i < size; i++)
			{
				free(matrix[i]);
				free(matriximg[i]);
			}
			printf("regcount = %d, on average %d per point.\n", regcount, regcount/size);
			free(matrix);
			free(matriximg);
		}
		else
		{
			for (i = 0; i < FASTSIZEBOUND; i++)
			{
				free(matrix[i]);
				free(matriximg[i]);
			}
			free(matrix);
			free(matriximg);
		}

		free(zx);
		free(zy);
		free(zz);
		free(zstarx);
		free(zstary);
		free(zstarz);
		free(delta);
		free(J);
		free(result);
		free(gradx);
		free(grady);
		free(gradz);
		free(regfactor);
		free(regfactorimg);
		return 0.;
	}

	if (size <= FASTSIZEBOUND)
	{
		if (fastflag >= (long int ) size * (long int) size )
		{
			tempC[0] = matrix[pi][pj];
			tempC[1] = matriximg[pi][pj];
			if (printregcount == 0)
			{
//				printf("The regcount is %d.\n", regcount);
				printregcount = 1;
			}
			return 0;
		}
		fastflag++;
	}
	else if ((pi < FASTSIZEBOUND) && (pj < FASTSIZEBOUND) && (fillinflag >= \
			(long int)(FASTSIZEBOUND) * (long int)(FASTSIZEBOUND)))
	{
		tempC[0] = matrix[pi][pj];
		tempC[1] = matriximg[pi][pj];
		return 0;
	}

	if (ini == 0)	//	FIRST TIME, THE ONLY TIME DOING INITIALIZATION FOR X, Y, Z POINTS //
	{
		tau = H * TAU_FACTOR;
		threshold = tau * tau;


//		term1 = 1.0/(4.0*PI*tau*radius);		//	1/8pi*tau(2/r)	//
//		term2 = (-5.0*tau)/(256.0*pow(radius,3.0)*PI);	//	-(1/pi)*(5/512)*(2/r^3)*tau	//
//		term3 = (-50.0*tau)/(1536.0*pow(radius,3.0)*PI);//	-(1/pi)*(25/1536)*(2/r^3)*tau	//
//		term4 = (pow(WAVE,2.0)*tau)/(24.0*PI*radius);	//	Helmholtz 3D, (k^2/48)*(k1+k2)*tau	//
//
//		if (INTERIOR_MODE==1)
//			regfactor = term1 + term2 + term3 + term4;
//		else
//			regfactor = -1.0*(term1 + term2 + term3 + term4);
//
//		regfactorimg = 0.0;

		zx = (double *) malloc(size * sizeof(double));
		zy = (double *) malloc(size * sizeof(double));
		zz = (double *) malloc(size * sizeof(double));
		zstarx = (double *) malloc(size * sizeof(double));
		zstary = (double *) malloc(size * sizeof(double));
		zstarz = (double *) malloc(size * sizeof(double));
		delta = (double *) malloc(size * sizeof(double));
		J = (double *) malloc(size * sizeof(double));
		gradx = (double *) malloc(size * sizeof(double));
		grady = (double *) malloc(size * sizeof(double));
		gradz = (double *) malloc(size * sizeof(double));
		result= (double *) malloc(size * sizeof(double));
		regfactor = (double * ) malloc(size * sizeof(double));
		regfactorimg = (double *) malloc(size * sizeof(double));
	
//		printf("Here after initialization.\n");

		if ( (zx==NULL)||(zy==NULL)||(zz==NULL) || (zstarx==NULL) || (zstary==NULL) ||(zstarz==NULL))
			printf("Memory allocation error....points.\n");
		if ( (delta==NULL)||(J==NULL)||(gradx==NULL)|| (grady==NULL)||(gradz==NULL) )
			printf("Memory allocation error....calcs.\n");
		if (result==NULL)
			printf("Memory allocation error....result.\n");

		for ( i = 0; i < SIZE; i++)
		{
			tempzy = INT_L + H * (double)(i);
		for ( j = 0; j < SIZE; j++)
		{
			tempzx = INT_L + H * (double)(j);
		for ( k = 0; k < SIZE; k++)
		{
			tempzz = INT_L + H * (double)(k);

//			thisdist = distance(tempzx, tempzy, tempzz);
			thisdist = DIST[i][j][k];			
			//	if not in tubular region, don't calculate.	//

			if (TRAD_ONESIDED == 0)
			{
				if (fabs(thisdist) > epsilon)
					continue;
			}
			else
			{
				if ( (thisdist > epsilon) || (thisdist < 0.) )
					continue;
			}

			zx[counti] = tempzx;
			zy[counti] = tempzy;
			zz[counti] = tempzz;
			tstarx = tempzx - thisdist * DIST_DX[i][j][k];
			tstary = tempzy - thisdist * DIST_DY[i][j][k];
			tstarz = tempzz - thisdist * DIST_DZ[i][j][k];
	
			zstarx[counti] = tstarx;
			zstary[counti] = tstary;
			zstarz[counti] = tstarz;

			if (TRAD_ONESIDED == 0)
				delta[counti] = cal_delta(epsilon, thisdist);
			else
				delta[counti] = 2. * cal_delta(epsilon, thisdist);
				
//			J[counti] = cal_J(J_M, tempzx, tempzy, tempzz);
			J[counti] = 1. + 2. * thisdist * MCURV[i][j][k] + thisdist*thisdist * GCURV[i][j][k];
//			gradx[counti] = cal_interpgrad(1, tstarx, tstary, tstarz);
//			grady[counti] = cal_interpgrad(2, tstarx, tstary, tstarz);
//			gradz[counti] = cal_interpgrad(3, tstarx, tstary, tstarz);

			indicator = cal_interpgrad(N,H,tstarz, tstary, tstarz, DIST_DX, DIST_DY, DIST_DZ, \
					&gradx[counti], &grady[counti], &gradz[counti]);
			if (indicator == 1)
			{
				printf("At (%lf, %lf, %lf), grad (%lf, %lf, %lf), dist %lf\n",
				tempzx, tempzy, tempzz, DIST_DX[i][j][k], DIST_DY[i][j][k], DIST_DZ[i][j][k], thisdist);
				printf("Ccomp Projected to (%lf, %lf, %lf), dist %lf\n", \
				tstarx, tstary, tstarz, cal_dist(DIST, tstarx, tstary, tstarz, N, H));
			}
			gradnorm =sqrt(gradx[counti]*gradx[counti]+grady[counti]*grady[counti]+gradz[counti]*gradz[counti]);
			gradx[counti] = gradx[counti]/gradnorm;
			grady[counti] = grady[counti]/gradnorm;
			gradz[counti] = gradz[counti]/gradnorm;

//			gradnorm =sqrt(gradx[counti]*gradx[counti]+grady[counti]*grady[counti]+gradz[counti]*gradz[counti]);
			if (fabs(gradnorm-1.) > 0.3)
			{
//				printf("point (%lf, %lf, %lf), grad (%lf, %lf, %lf), norm %lf, dist %lf\n", \
//				tempzx, tempzy, tempzz, DIST_DX[i][j][k], DIST_DY[i][j][k], DIST_DZ[i][j][k], \
//				sqrt(DIST_DX[i][j][k]*DIST_DX[i][j][k] + DIST_DY[i][j][k]*DIST_DY[i][j][k] \
//				+ DIST_DZ[i][j][k]*DIST_DZ[i][j][k]), thisdist);
//
				printf("At (%lf, %lf,%lf), grad (%lf, %lf, %lf) norm %lf.\n", \
				tstarx, tstary, tstarz, gradx[counti], grady[counti], gradz[counti], gradnorm);
//				errorflag++;
//				if (errorflag > 100)
//					exit(0);
			}
			///	REGULARIZATION		////////////////

//			cur_meancurv = cal_curv(DIST, tstarx, tstary, tstarz, N, H);
//			cur_gcurv = cal_GaussianCurv(DIST, tstarx, tstary, tstarz, N, H);

			cal_BothCurvs( tstarx, tstary, tstarz, H, MCURV, GCURV, &cur_meancurv, &cur_gcurv);

			regterm1 = cur_meancurv/(4.*PI*tau);		//	(k1+k2)/8pitau

			//	(5tau/512pi)) (k1+k2)(k1^2 + k1k2 + k2^2)	//
			regterm2 = -5.*tau*(2.*cur_meancurv)*(4.*cur_meancurv*cur_meancurv-cur_gcurv)/(512.*PI);	
			//	(25tau/1536pi) (k1+k2) * k1k2	//
			regterm3 = -25.*tau*(2.*cur_meancurv)*cur_gcurv/(1536.*PI);
			//	k^2 tau / (48pi) * (k1+k2)
			regterm4 = WAVE*WAVE * tau * 2.* cur_meancurv / (48.*PI);

			if (INTERIOR_MODE == 1)
				regfactor[counti] = regterm1 + regterm2 + regterm3 + regterm4;
			else
				regfactor[counti] = -1. * (regterm1 + regterm2 + regterm3 + regterm4);
			regfactorimg[counti] = 0.0;


			result[counti] = HCUBED * delta[counti] * J[counti];
	
			counti++;	// in the end should total up to size	//
		}
		}
		}

		if (size <= FASTSIZEBOUND)
		{
			matrix = (double **) malloc(size * sizeof(double *));
			matriximg = (double **) malloc(size * sizeof(double *));
			for ( i = 0; i < size; i++)
			{
				matrix[i] = (double *) malloc(size * sizeof(double));
				matriximg[i] = (double *) malloc(size * sizeof(double));
				if ( (matrix[i]==NULL) || (matriximg[i]==NULL) )
				{
					printf("Memory allocation error for matrix[%d].\n", i);
					exit(0);
				}
			}
		}
		else
		{
			matrix = (double **) malloc(FASTSIZEBOUND * sizeof(double *));
			matriximg = (double **) malloc(FASTSIZEBOUND * sizeof(double *));
			for (i = 0; i < FASTSIZEBOUND; i++)
			{
				matrix[i] = (double *) malloc(FASTSIZEBOUND * sizeof(double));
				matriximg[i] = (double *) malloc(FASTSIZEBOUND * sizeof(double));
				if ( (matrix[i]==NULL) || (matriximg[i]==NULL) )
				{
					printf("Memory allocation error for matrix[%d].\n", i);
					exit(0);
				}
			}
		}

		ini = 1;	//	DENOTING DATA INITIALIZED	//
	}


	starxi = zstarx[pi];
	staryi = zstary[pi];
	starzi = zstarz[pi];
	starxj = zstarx[pj];
	staryj = zstary[pj];
	starzj = zstarz[pj];

	xmy1 = starxi - starxj;
	xmy2 = staryi - staryj;
	xmy3 = starzi - starzj;
	xmynorm = sqrt( xmy1*xmy1 + xmy2*xmy2 + xmy3*xmy3 );
	
	if (xmynorm > tau)
	{
		recipxmynorm = 1./xmynorm;
		parameter = WAVE * xmynorm;		//	k|x-y|
		SineValue = sin(parameter);		//	sin(k|x-y|)
		CosValue = cos(parameter);		//	cos(k|x-y|)
		extSineValue = SineValue * recipxmynorm;	//	sin(k|x-y|)/|x-y|
		extCosValue = CosValue * recipxmynorm;		//	cos(k|x-y|)/|x-y|
		denom = recipfourpi*recipxmynorm*recipxmynorm;	//	1./4pi|x-y|^2
		
		commonterm1 = -1.* (WAVE * SineValue + extCosValue) * denom;
		commonterm2 = (WAVE * CosValue - extSineValue) * denom;
		
		if (DIRICHLET == 1)
			inner_pd = gradx[pj] * xmy1 + grady[pj] * xmy2 + gradz[pj] * xmy3;
		else
			inner_pd = gradx[pi] * xmy1 + grady[pi] * xmy2 + gradz[pi] * xmy3;

//		partialx = commonterm1 * (xmy1);
//		partialy = commonterm1 * (xmy2);
//		partialz = commonterm1 * (xmy3);
//		partialximg = commonterm2 * (xmy1);
//		partialyimg = commonterm2 * (xmy2);
//		partialzimg = commonterm2 * (xmy3);


//	partialx = cal_partial(1, starxi, staryi, starzi, starxj, staryj, starzj);// x coordinate real	//
//	partialy = cal_partial(2, starxi, staryi, starzi, starxj, staryj, starzj);// y coordinate real	//
//	partialz = cal_partial(3, starxi, staryi, starzi, starxj, staryj, starzj);// z coordinate real	//
//	partialximg = cal_partial(4, starxi, staryi, starzi, starxj, staryj, starzj);// x coordinate imaginary //
//	partialyimg = cal_partial(5, starxi, staryi, starzi, starxj, staryj, starzj);// y coordinate imaginary //
//	partialzimg = cal_partial(6, starxi, staryi, starzi, starxj, staryj, starzj);// z coordinate imaginary //


		//	for interior problem	//
		if (INTERIOR_MODE==1)	//	outer normal	//
		{
			partial = (-1.0) * commonterm1 * inner_pd;
			partialimg = (-1.0) * commonterm2 * inner_pd;
		}
		//	for exterior problem	//
		else	
		{
			partial = commonterm1 * inner_pd;		
			partialimg = commonterm2 * inner_pd;
		}
	}
	else
	{
		regcount += 1;
		if (DIRICHLET == 1)
			partial = regfactor[j];
		else
			partial = regfactor[i];
		partialimg = 0.0;
	}

	//	calculate matrix elements	//

	tempC[0] = result[pj] * partial ;		//	( h^3*deltaj*Jj ) * partial	//
	tempC[1] = result[pj] * partialimg;
//	if (REG_MODE==1)
//	{
//		tempC[0] = result[pj] * partial_reg;
//		tempC[1] = result[pj] * partialimg_reg;
//	}

	//	C = B + 1/2 I	//
	if (pi == pj)
	{
		if (DIRICHLET == 1)
			tempC[0] += 0.5;	//(delta*cal_J(J_M, wx, wy));
		else
			tempC[0] += 0.5;
//		if (REG_MODE==1)
//		{
//			C_reg += 0.5;	//(delta*cal_J(J_M, wx, wy));
//		}
//		C += 1.0;	//	+1 is for Helmholtz, + 0.5 is for Laplace->which is wrong.
//		C_reg += 1.0;	//	+1 is for Helmholtz, + 0.5 is for Laplace. Both should be + 0.5
	}

	if (size <= FASTSIZEBOUND)
	{
		matrix[pi][pj] = tempC[0];
		matriximg[pi][pj] = tempC[1];
	}
	else if ( (pi < FASTSIZEBOUND) && (pj < FASTSIZEBOUND) && ((long int)fillinflag < \
			(long int) (FASTSIZEBOUND) * (long int)(FASTSIZEBOUND)) )
	{
		matrix[pi][pj] = tempC[0];
		matriximg[pi][pj] = tempC[1];
		fillinflag++;
	}


//	if ( (matrix[pi][pj] != matrix[pi][pj]) || (matrix[pi][pj]+1 == matrix[pi][pj]) ||
//		(matriximg[pi][pj] != matriximg[pi][pj]) || (matriximg[pi][pj]+1 == matriximg[pi][pj]) )
//	{
//		printf("i = %d, j = %d, partial = %lf, matrix[pi][pj] = (%lf, %lf), fillflag = %ld, fastflag = %ld.\n", \
//		pi, pj, partial, matrix[pi][pj], matriximg[pi][pj], fillinflag, fastflag);
//		exit(0);
//	}
	if ((tempC[0]!=tempC[0])||(tempC[0]+1==tempC[0])||(tempC[1]!=tempC[1])||(tempC[1]+1==tempC[1]))
	{
		printf("i = %d, j = %d, partial = %.10lf, C = %.10lf, Cimg = %.10lf.\n", \
			pi, pj, partial, tempC[0], tempC[1]);
		exit(0);
	}

	return tempC[0];

//	if (REG_MODE==0)
//		return C;
//	else
//		return C_reg;
}


double omp_get_Ccomp(int size, int pi, int pj, double epsilon, double tempC[])
{
	FILE *fp, *fp2;

	int i, j, k, indicator;
	int counti = 0;
	int countj = 0;

	static double *zx, *zy, *zz, *wx, *wy, *wz;
	static double *zstarx, *zstary, *zstarz, *wstarx, *wstary, *wstarz;
	static double *J, *delta;
	static double *gradx, *grady, *gradz;
	static double *result;
	static double *regfactor, *regfactorimg;
	static double tau, threshold;
	static double term1, term2, term3, term4, thisdist;
	double regterm1, regterm2, regterm3, regterm4, gradnorm;

	static double *temp;


	double tempzx, tempzy, tempzz, tstarx, tstary, tstarz;
	double starxi, staryi, starzi, starxj, staryj, starzj;
	double C, C_reg;
	double partial, partial_reg, partialimg, partialimg_reg;
	double partialx, partialy, partialz, partialximg, partialyimg, partialzimg;
	double dxx, dyy, dzz, dxy, dxz, dyz, cur_meancurv, cur_gcurv;

	static const double recipfourpi = 1./(4.*PI);
	double xmy1, xmy2, xmy3, xmynorm, recipxmynorm, parameter, CosValue, SineValue, extCosValue, extSineValue;
	double commonterm1, commonterm2, denom, inner_pd;

	static int ini = 0;
	static int flag = 0;
	static int regcount = 0;
	static int printregcount = 0;
	static int errorflag = 0;

	if ( (pi == -1) && (pj == -1) )
	{
		printf("regcount = %d, on average %d per point.\n", regcount, regcount/size);
		free(zx);
		free(zy);
		free(zz);
		free(zstarx);
		free(zstary);
		free(zstarz);
		free(delta);
		free(J);
		free(result);
		free(gradx);
		free(grady);
		free(gradz);
		free(regfactor);
		free(regfactorimg);
		return 0.;
	}
	if (ini == 0)	//	FIRST TIME, THE ONLY TIME DOING INITIALIZATION FOR X, Y, Z POINTS //
	{
		tau = H * TAU_FACTOR;
		threshold = tau * tau;


//		term1 = 1.0/(4.0*PI*tau*radius);		//	1/8pi*tau(2/r)	//
//		term2 = (-5.0*tau)/(256.0*pow(radius,3.0)*PI);	//	-(1/pi)*(5/512)*(2/r^3)*tau	//
//		term3 = (-50.0*tau)/(1536.0*pow(radius,3.0)*PI);//	-(1/pi)*(25/1536)*(2/r^3)*tau	//
//		term4 = (pow(WAVE,2.0)*tau)/(24.0*PI*radius);	//	Helmholtz 3D, (k^2/48)*(k1+k2)*tau	//
//
//		if (INTERIOR_MODE==1)
//			regfactor = term1 + term2 + term3 + term4;
//		else
//			regfactor = -1.0*(term1 + term2 + term3 + term4);
//
//		regfactorimg = 0.0;

		zx = (double *) malloc(size * sizeof(double));
		zy = (double *) malloc(size * sizeof(double));
		zz = (double *) malloc(size * sizeof(double));
		zstarx = (double *) malloc(size * sizeof(double));
		zstary = (double *) malloc(size * sizeof(double));
		zstarz = (double *) malloc(size * sizeof(double));
		delta = (double *) malloc(size * sizeof(double));
		J = (double *) malloc(size * sizeof(double));
		gradx = (double *) malloc(size * sizeof(double));
		grady = (double *) malloc(size * sizeof(double));
		gradz = (double *) malloc(size * sizeof(double));
		regfactor = (double *) malloc(size * sizeof(double));
		regfactorimg = (double *) malloc(size * sizeof(double));
		result= (double *) malloc(size * sizeof(double));
	
//		printf("Here after initialization.\n");

		if ( (zx==NULL)||(zy==NULL)||(zz==NULL) || (zstarx==NULL) || (zstary==NULL) ||(zstarz==NULL))
			printf("Memory allocation error....points.\n");
		if ( (delta==NULL)||(J==NULL)||(gradx==NULL)|| (grady==NULL)||(gradz==NULL) )
			printf("Memory allocation error....calcs.\n");
		if (result==NULL)
			printf("Memory allocation error....result.\n");

		for ( i = 0; i < SIZE; i++)
		{
			tempzy = INT_L + H * (double)(i);
		for ( j = 0; j < SIZE; j++)
		{
			tempzx = INT_L + H * (double)(j);
		for ( k = 0; k < SIZE; k++)
		{
			tempzz = INT_L + H * (double)(k);

//			thisdist = distance(tempzx, tempzy, tempzz);
			thisdist = DIST[i][j][k];			
			//	if not in tubular region, don't calculate.	//

			if (TRAD_ONESIDED == 0)
			{
				if (fabs(thisdist) > epsilon)
					continue;
			}
			else
			{
				if ( (thisdist > epsilon) || (thisdist < 0.) )
					continue;
			}

			zx[counti] = tempzx;
			zy[counti] = tempzy;
			zz[counti] = tempzz;

			tstarx = tempzx - thisdist * DIST_DX[i][j][k];
			tstary = tempzy - thisdist * DIST_DY[i][j][k];
			tstarz = tempzz - thisdist * DIST_DZ[i][j][k];
	
			zstarx[counti] = tstarx;
			zstary[counti] = tstary;
			zstarz[counti] = tstarz;

			if (TRAD_ONESIDED == 0)
				delta[counti] = cal_delta(epsilon, thisdist);
			else
				delta[counti] = 2. * cal_delta(epsilon, thisdist);
				
//			J[counti] = cal_J(J_M, tempzx, tempzy, tempzz);

//			dxx = secondpartialxx(tempzx, tempzy, tempzz);
//			dyy = secondpartialyy(tempzx, tempzy, tempzz);
//			dzz = secondpartialzz(tempzx, tempzy, tempzz);
//			dxy = secondpartialxy(tempzx, tempzy, tempzz);
//			dxz = secondpartialxz(tempzx, tempzy, tempzz);
//			dyz = secondpartialyz(tempzx, tempzy, tempzz);
//			cur_meancurv = -0.5 * (dxx + dyy + dzz);
//			cur_gcurv = dxx*dyy + dyy*dzz + dxx*dzz - dxy*dxy - dxz*dxz - dyz*dyz;
			J[counti] = 1 + 2. * thisdist * MCURV[i][j][k] + thisdist*thisdist * GCURV[i][j][k];


//			gradx[counti] = cal_interpgrad(1, tstarx, tstary, tstarz);
//			grady[counti] = cal_interpgrad(2, tstarx, tstary, tstarz);
//			gradz[counti] = cal_interpgrad(3, tstarx, tstary, tstarz);
			indicator = cal_interpgrad(N, H, tstarx, tstary, tstarz, DIST_DX, DIST_DY, DIST_DZ, \
					&gradx[counti], &grady[counti], &gradz[counti]);
			if (indicator == 1)
			{
				printf("At (%lf, %lf, %lf), grad (%lf, %lf, %lf), dist %lf\n",
				tempzx, tempzy, tempzz, DIST_DX[i][j][k], DIST_DY[i][j][k], DIST_DZ[i][j][k], thisdist);
				printf("ompProjected to (%lf, %lf, %lf), dist %lf\n", \
				tstarx, tstary, tstarz, cal_dist(DIST, tstarx, tstary, tstarz, N, H));
			}
			gradnorm =sqrt(gradx[counti]*gradx[counti]+grady[counti]*grady[counti]+gradz[counti]*gradz[counti]);
			gradx[counti] = gradx[counti]/gradnorm;
			grady[counti] = grady[counti]/gradnorm;
			gradz[counti] = gradz[counti]/gradnorm;

//			gradnorm =sqrt(gradx[counti]*gradx[counti]+grady[counti]*grady[counti]+gradz[counti]*gradz[counti]);
			if (fabs(gradnorm-1.) > 0.1)
			{
				printf("point (%lf, %lf, %lf), grad (%lf, %lf, %lf), norm %lf, dist %lf\n", \
				tempzx, tempzy, tempzz, DIST_DX[i][j][k], DIST_DY[i][j][k], DIST_DZ[i][j][k], \
				sqrt(DIST_DX[i][j][k]*DIST_DX[i][j][k]+DIST_DY[i][j][k]*DIST_DY[i][j][k] \
				+ DIST_DZ[i][j][k]*DIST_DZ[i][j][k]), thisdist);

				printf("At (%lf, %lf,%lf), grad (%lf, %lf, %lf) norm %lf.\n", \
				tstarx, tstary, tstarz, gradx[counti], grady[counti], gradz[counti], gradnorm);
				errorflag++;
				if (errorflag > 100)
					exit(0);
			}
			///	REGULARIZATION		////////////////

//			cur_meancurv = cal_curv(DIST, tstarx, tstary, tstarz, N, H);
//			cur_gcurv = cal_GaussianCurv(DIST, tstarx, tstary, tstarz, N, H);
			cal_BothCurvs(tstarx, tstary, tstarz, H, MCURV, GCURV, &cur_meancurv, &cur_gcurv);

			regterm1 = cur_meancurv/(4.*PI*tau);		//	(k1+k2)/8pitau
	
			//	(5tau/512pi)) (k1+k2)(k1^2 + k1k2 + k2^2)	//
			regterm2 = -5.*tau*(2.*cur_meancurv)*(4.*cur_meancurv*cur_meancurv-cur_gcurv)/(512.*PI);	
			//	(25tau/1536pi) (k1+k2) * k1k2	//
			regterm3 = -25.*tau*(2.*cur_meancurv)*cur_gcurv/(1536.*PI);
			//	k^2 tau / (48pi) * (k1+k2)
			regterm4 = WAVE*WAVE * tau * 2.* cur_meancurv / (48.*PI);

			if (INTERIOR_MODE == 1)	
				regfactor[counti] = regterm1 + regterm2 + regterm3 + regterm4;
			else
				regfactor[counti] = -1. * (regterm1 + regterm2 + regterm3 + regterm4);
			regfactorimg[counti] = 0.0;
	
			result[counti] = HCUBED * delta[counti] * J[counti];
	
			counti++;	// in the end should total up to size	//
		}	//	k
		}	//	j	
		}	//	i
		ini = 1;
		if (size != counti)
		{
			printf("Size %d, counti%d.\n", size, counti);
			exit(0);
		}
	}



	starxi = zstarx[pi];
	staryi = zstary[pi];
	starzi = zstarz[pi];
	starxj = zstarx[pj];
	staryj = zstary[pj];
	starzj = zstarz[pj];


	xmy1 = starxi - starxj;
	xmy2 = staryi - staryj;
	xmy3 = starzi - starzj;
	xmynorm = sqrt( xmy1*xmy1 + xmy2*xmy2 + xmy3*xmy3 );
	
	if (xmynorm > tau)
	{
		recipxmynorm = 1./xmynorm;
		parameter = WAVE * xmynorm;		//	k|x-y|
		SineValue = sin(parameter);		//	sin(k|x-y|)
		CosValue = cos(parameter);		//	cos(k|x-y|)
		extSineValue = SineValue * recipxmynorm;	//	sin(k|x-y|)/|x-y|
		extCosValue = CosValue * recipxmynorm;		//	cos(k|x-y|)/|x-y|
		denom = recipfourpi*recipxmynorm*recipxmynorm;	//	1./4pi|x-y|^2
		
		commonterm1 = -1.* (WAVE * SineValue + extCosValue) * denom;
		commonterm2 = (WAVE * CosValue - extSineValue) * denom;
		

//		partialx = commonterm1 * (xmy1);
//		partialy = commonterm1 * (xmy2);
//		partialz = commonterm1 * (xmy3);
//		partialximg = commonterm2 * (xmy1);
//		partialyimg = commonterm2 * (xmy2);
//		partialzimg = commonterm2 * (xmy3);

		if (DIRICHLET == 1)
			inner_pd = xmy1 * gradx[pj] + xmy2 * grady[pj] + xmy3 * gradz[pj];
		else
			inner_pd = xmy1 * gradx[pi] + xmy2 * grady[pi] + xmy3 * gradz[pi];


//	partialx = cal_partial(1, starxi, staryi, starzi, starxj, staryj, starzj);// x coordinate real	//
//	partialy = cal_partial(2, starxi, staryi, starzi, starxj, staryj, starzj);// y coordinate real	//
//	partialz = cal_partial(3, starxi, staryi, starzi, starxj, staryj, starzj);// z coordinate real	//
//	partialximg = cal_partial(4, starxi, staryi, starzi, starxj, staryj, starzj);// x coordinate imaginary //
//	partialyimg = cal_partial(5, starxi, staryi, starzi, starxj, staryj, starzj);// y coordinate imaginary //
//	partialzimg = cal_partial(6, starxi, staryi, starzi, starxj, staryj, starzj);// z coordinate imaginary //


		//	for interior problem	//
		if (INTERIOR_MODE==1)	//	outer normal	//
		{
			partial = (-1.0) * commonterm1 * inner_pd;
			partialimg = (-1.0) * commonterm2 * inner_pd;
		}
		//	for exterior problem	//
		else	
		{
			partial = commonterm1 * inner_pd;
			partialimg = commonterm2 * inner_pd;
		}
	}
	else
	{
		regcount += 1;
		if (DIRICHLET == 1)
			partial = regfactor[j];
		else
			partial = regfactor[i];
		partialimg = 0.0;
	}
	
//	if ( ( (starxi-starxj)*(starxi-starxj) + (staryi-staryj)*(staryi-staryj) + \
//		(starzi-starzj)*(starzi-starzj) ) < threshold )
//	{
//		partial_reg = regfactor;
//		partialimg_reg = regfactorimg;
//		regcount += 1;
//	}
//	else
//	{
//		partial_reg = partial;
//		partialimg_reg = partialimg;
//	}

	//	calculate matrix elements	//

	tempC[0] = result[pj] * partial ;		//	( h^3*deltaj*Jj ) * partial	//
	tempC[1] = result[pj] * partialimg;
//	if (REG_MODE==1)
//	{
//		tempC[0] = result[pj] * partial_reg;
//		tempC[1] = result[pj] * partialimg_reg;
//	}

	//	C = B + 1/2 I	//
	if (pi == pj)
	{
		if (DIRICHLET == 1)
			tempC[0] += 0.5;	//(delta*cal_J(J_M, wx, wy));
		else
			tempC[0] += 0.5;
	}


	if ((tempC[0]!=tempC[0])||(tempC[0]+1==tempC[0])||(tempC[1]!=tempC[1])||(tempC[1]+1==tempC[1]))
	{
		printf("i = %d, j = %d, partial = %.10lf, C = %.10lf, Cimg = %.10lf.\n", \
			pi, pj, partial, tempC[0], tempC[1]);
		exit(0);
	}

	return tempC[0];

}



double cal_u(double x, double y, double z)
{
	int l, m;
	const static double a[4][4] = COEFFICIENTA;
	const static double b[4][4] = COEFFICIENTB;
	static double Lambda[6][6] = ROOTBESSEL;
	static double Lambda2[6][6] = ROOTBESSEL;
	static int lambda_flag = 0;

	double r, temptheta, theta, tempphi, phi, templ;
	double lsum = 0.0;
	double sum = 0.0;
//	double a[4] = COEFFICIENTA;
//	double b[4] = COEFFICIENTB;
	static int flag = 0;
	static int incidental_flag = 0;

	if (INCIDENTALWAVE == 1)
	{
		if (incidental_flag == 0)
			printf("Incidental test.\n");
		incidental_flag = 1;
		return 0.0;
	}


//	return 3.0;			///	DEBUG

	
	if ( (LAPLACE==0) && (NONHOM_MODE==1) )		//	Total wave Dirichlet hard condition.	//
		return 0.0;


	r = sqrt(pow(x-Cx,2.0) + pow(y-Cy,2.0) + pow(z-Cz, 2.0));	// note that not necessarily radius!

	temptheta = acos((z-Cz)/r);
	theta = temptheta;
	tempphi = atan((y-Cy)/(x-Cx));	// note the range of inverse trig functions is tricky	//
	if (tempphi!=tempphi)
		tempphi = 0.0;


	if ((x >= Cx) && (y < Cy))	//	4th quadrant
		phi = tempphi + 2.0*PI;
	else if ((x >= Cx) && (y >= Cy))	//	1st quadrant
		phi = tempphi;
	else
		phi = tempphi + PI;




	if (LAPLACE==0)
	{

		if (INTERIOR_MODE == 1)
		{
			return sin(WAVE*r)/(WAVE*r);		//	besselj(0, WAVE*r);
		}
		else
		{
//			return sin(WAVE*r)/(WAVE*r);
			return -1.0*cos(WAVE*r)/(WAVE*r);	//	bessely(0, WAVE*r);
		}
	}
	else
	{
		if (INTERIOR_MODE==0)
		{
			for ( l = 0; l <= 3; l++)
			{
				templ = pow(r, -1.0*(double)(l)-1.0);
				for ( m = 0; m <= l; m++)
					sum += templ * ( a[l][m]*cos(m*phi) + b[l][m]*sin(m*phi) ) * cal_f(m,l,cos(theta));
			}
			return sum;
		}
		else
		{
			//	THE FORMULA FOR Ue = sigma r^l	sigma Fourrier*Legendre(cos(theta))	//
			for ( l = 0; l <= 3; l++)
			{
				templ = pow(r, (double)(l));
				for ( m = 0; m <= l; m++)
					sum += templ * ( a[l][m]*cos(m*phi) + b[l][m]*sin(m*phi) ) * cal_f(m,l,cos(theta));
			}
				return sum;
		}
	}


//	if (flag < 30)
//	{	
//		printf("r = %.10E, theta = %.10E, sum = %.10E.\n", r, theta, sum);
//		flag++;
//	}
	return sum;
}

double cal_uimg(double x, double y, double z)
{
	const static double a[4][4] = COEFFICIENTA;
	const static double b[4][4] = COEFFICIENTB;
	static double Lambda[6][6] = ROOTBESSEL;
	static double Lambda2[6][6] = ROOTBESSEL;
	static int lambda_flag = 0;

	double r, temptheta, theta, tempphi, phi, templ;
	double lsum = 0.0;
	double sum = 0.0;
//	double a[4] = COEFFICIENTA;
//	double b[4] = COEFFICIENTB;
	static int flag = 0;

//	return 3.0;			///	DEBUG

	if (INCIDENTALWAVE==1)
		return 0.0;
	
	if ( (LAPLACE==0) && (NONHOM_MODE==1) )		//	Total wave Dirichlet hard condition.	//
		return 0.0;


	r = sqrt(pow(x-Cx,2.0) + pow(y-Cy,2.0) + pow(z-Cz, 2.0));	// note that not necessarily radius!

	temptheta = acos((z-Cz)/r);
	theta = temptheta;
	tempphi = atan((y-Cy)/(x-Cx));	// note the range of inverse trig functions is tricky	//
	if (tempphi!=tempphi)
		tempphi = 0.0;


	if ((x >= Cx) && (y < Cy))	//	4th quadrant
		phi = tempphi + 2.0*PI;
	else if ((x >= Cx) && (y >= Cy))	//	1st quadrant
		phi = tempphi;
	else
		phi = tempphi + PI;

//	if (flag==0)
//	{
//		printf("uimg has not been completed yet.\n");
//		flag++;
//	}
	return -1.0 * sin(WAVE * r)/(WAVE * r);
	return -1.0*cos(WAVE*r)/(WAVE*r);	//	bessely(0, WAVE*r);
	return 0.0;
}


double cal_partialu(double x, double y, double z)
{ 
	double r, r2, temptheta, theta, phi, tempphi;


	if (INCIDENTALWAVE == 1)
	{
		double Dx, Dy, Dz, Dnorm, inner_pd, inner_pd2;
		cal_interpgrad(N, H, x, y, z, DIST_DX, DIST_DY, DIST_DZ, &Dx, &Dy, &Dz);
		Dnorm = sqrt(Dx*Dx + Dy*Dy + Dz*Dz);
		Dx = Dx/Dnorm;
		Dy = Dy/Dnorm;
		Dz = Dz/Dnorm;
		inner_pd = x * PLANEWAVEX + y * PLANEWAVEY + z * PLANEWAVEZ;
		inner_pd2 = -1. * (PLANEWAVEX * Dx + PLANEWAVEY * Dy + PLANEWAVEZ * Dz);
		return 10. * WAVE * inner_pd2 * sin(WAVE * inner_pd);
	}


	r = sqrt(pow(x-Cx,2.0) + pow(y-Cy,2.0) + pow(z-Cz, 2.0));	// note that not necessarily radius!

	temptheta = acos((z-Cz)/r);
	theta = temptheta;
	tempphi = atan((y-Cy)/(x-Cx));	// note the range of inverse trig functions is tricky	//
	if (tempphi!=tempphi)
		tempphi = 0.0;

	if ((x >= Cx) && (y < Cy))	//	4th quadrant
		phi = tempphi + 2.0*PI;
	else if ((x >= Cx) && (y >= Cy))	//	1st quadrant
		phi = tempphi;
	else
		phi = tempphi + PI;

	if (LAPLACE == 1)
	{
		printf("Shouldn't need to use this program to solve Laplace Neumann Problem yet.\n");
		exit(0);
	}
	else
	{
		if (INTERIOR_MODE == 1)
		{
			printf("Not yet.\n");
			return BesselJ1(WAVE*r)/BesselY0(WAVE*radius);	
		}
		else
		{
			return ( sin(WAVE*r) / r ) + ( cos(WAVE * r)/(WAVE * r * r) );
		}
	}
	printf("Why here in cal_partialu?\n");
	exit(0);

}


double cal_partialuimg(double x, double y, double z)
{
	int indicator;
	double r, r2, temptheta, theta, phi, tempphi;

	if (INCIDENTALWAVE == 1)
	{
		double Dx, Dy, Dz, Dnorm, inner_pd, inner_pd2;
		indicator = cal_interpgrad(N, H, x, y, z, DIST_DX, DIST_DY, DIST_DZ, &Dx, &Dy, &Dz);
		if (indicator == 1)
		{
			printf("At (%lf, %lf, %lf), dist %lf\n",
			x, y, z, cal_dist(DIST, x, y, z, N, H));
//			printf("Projected to (%lf, %lf, %lf), dist %lf\n", \
//			zstarx, zstary, zstarz, cal_dist(DIST, zstarx, zstary, zstarz, N, H));
		}
		Dnorm = sqrt(Dx*Dx + Dy*Dy + Dz*Dz);
		Dx = Dx/Dnorm;
		Dy = Dy/Dnorm;
		Dz = Dz/Dnorm;
		inner_pd = x * PLANEWAVEX + y * PLANEWAVEY + z * PLANEWAVEZ;
		inner_pd2 = -1. * (PLANEWAVEX * Dx + PLANEWAVEY * Dy + PLANEWAVEZ * Dz);
		return -10. * WAVE * inner_pd2 * cos(WAVE * inner_pd);
	}

	r = sqrt(pow(x-Cx,2.0) + pow(y-Cy,2.0) + pow(z-Cz, 2.0));	// note that not necessarily radius!

	temptheta = acos((z-Cz)/r);
	theta = temptheta;
	tempphi = atan((y-Cy)/(x-Cx));	// note the range of inverse trig functions is tricky	//
	if (tempphi!=tempphi)
		tempphi = 0.0;

	if ((x >= Cx) && (y < Cy))	//	4th quadrant
		phi = tempphi + 2.0*PI;
	else if ((x >= Cx) && (y >= Cy))	//	1st quadrant
		phi = tempphi;
	else
		phi = tempphi + PI;

	if (LAPLACE == 1)
	{
		printf("Shouldn't need to use this program to solve Laplace Neumann Problem yet.\n");
		exit(0);
	}
	else
	{
		if (INTERIOR_MODE == 1)
		{
			printf("Not yet.\\n");
		}
		else
		{
			return (sin(WAVE*r)/(WAVE * r * r)) - (cos(WAVE*r)/r );
		}
	}
	printf("Why here in cal_partialuimg?\n");
	exit(0);
}





double cal_vdl(double x, double y, double z)
{
	int l, m;
	const static double a[4][4] = COEFFICIENTA;
	const static double b[4][4] = COEFFICIENTB;
	double r, temptheta, theta, tempphi, phi, templ;
	double lsum = 0.0;
	double sum = 0.0;
//	double a[4] = COEFFICIENTA;
//	double b[4] = COEFFICIENTB;
	static int flag = 0;

//	return 3.0;			///	DEBUG

	
	r = sqrt(pow(x-Cx,2.0) + pow(y-Cy,2.0) + pow(z-Cz, 2.0));	// note that not necessarily radius!

	temptheta = acos((z-Cz)/r);
	theta = temptheta;
	tempphi = atan((y-Cy)/(x-Cx));	// note the range of inverse trig functions is tricky	//
	if (tempphi!=tempphi)
		tempphi = 0.0;


	if ((x >= Cx) && (y < Cy))	//	4th quadrant
		phi = tempphi + 2.0*PI;
	else if ((x >= Cx) && (y >= Cy))	//	1st quadrant
		phi = tempphi;
	else
		phi = tempphi + PI;



	if (LAPLACE==0)
	{

		if (INTERIOR_MODE == 1)
		{
			return 0.0;
			return BesselY0(WAVE*r)/BesselY0(WAVE*radius);
		}
		else
		{
			return BesselY0(WAVE*r);
		}
	}
	else
	{
		if (INTERIOR_MODE==0)
		{
			for ( l = 0; l <= 3; l++)
			{
				templ = (-1.0) * pow(radius,-2.0*(double)(l)-2.0) * pow(r, (double)(l)+1.0);
				for ( m = 0; m <= l; m++)
					sum += templ * ( a[l][m]*cos(m*phi) + b[l][m]*sin(m*phi) ) * cal_f(m,l,cos(theta));
			}
			return sum;
		}
		else
		{
			//	THE FORMULA FOR Ue = sigma r^l	sigma Fourrier*Legendre(cos(theta))	//
			for ( l = 0; l <= 3; l++)
			{
				templ = (-1.0) * pow(radius,2.0*(double)(l)+1.0) * \
					pow(r, -1.0*(double)(l)-1.0) *((double)(l)/(double)(l+1.0));
				for ( m = 0; m <= l; m++)
					sum += templ * ( a[l][m]*cos(m*phi) + b[l][m]*sin(m*phi) ) * cal_f(m,l,cos(theta));
			}
			return sum;
		}
	}


//	if (flag==0)
//	{	
//		printf("r = %.10E, theta = %.10E, sum = %.10E.\n", r, theta, sum);
//		flag++;
//	}
	return sum;
}

double cal_vdlimg(double x, double y, double z)
{
	printf("vdlimg has not been completed yet.\n");
	return 0.0;
}


