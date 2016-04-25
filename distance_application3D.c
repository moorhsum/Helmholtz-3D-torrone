
#include "distance_application3D.h"

double GaussianCurv(double ***dist, int i, int j, int k, int N, double H)
{
	double dxx, dyy, dzz, dxy, dxz, dyz;
	double crossterm, diagterm;

	dxx = XDER2(dist, i, j, k, N, H);
	dyy = YDER2(dist, i, j, k, N, H);
	dzz = ZDER2(dist, i, j, k, N, H);

	dxy = secondpartialxy(dist, i, j, k, N, H);
	dxz = secondpartialxz(dist, i, j, k, N, H);
	dyz = secondpartialyz(dist, i, j, k, N, H);

	diagterm = dxx*dyy + dyy*dzz + dxx*dzz;
	crossterm = dxy*dxy + dyz*dyz + dxz*dxz;

	return diagterm - crossterm;
}

double MeanCurv(double ***dist, int i, int j, int k, int N, double H)
{
	double dxx, dyy, dzz;

	dxx = XDER2(dist, i, j, k, N, H);
	dyy = YDER2(dist, i, j, k, N, H);
	dzz = ZDER2(dist, i, j, k, N, H);

	return -0.5 * (dxx+dyy+dzz);
}


int BothCurvs(double ***dist, int i, int j, int k, int N, double H, double *meancurv, double *gausscurv)
{
	double dxx, dyy, dzz, dxy, dxz, dyz;
	double crossterm, diagterm;

	dxx = XDER2(dist, i, j, k, N, H);
	dyy = YDER2(dist, i, j, k, N, H);
	dzz = ZDER2(dist, i, j, k, N, H);

	dxy = secondpartialxy(dist, i, j, k, N, H);
	dxz = secondpartialxz(dist, i, j, k, N, H);
	dyz = secondpartialyz(dist, i, j, k, N, H);

	diagterm = dxx*dyy + dyy*dzz + dxx*dzz;
	crossterm = dxy*dxy + dyz*dyz + dxz*dxz;

	*meancurv = -0.5 * (dxx+dyy+dzz);
	*gausscurv = diagterm - crossterm;
	return 0;
}

int ASSIGN_GRADIENT(int N, double H, double ***dist, double ***graddx, double ***graddy, double ***graddz, int weno_mode)
{
	int size;
	int i, j, k;

	size = N+1;

//	#pragma omp parallel for private(i,j,k) schedule(static, 2)
	for ( i = 0; i < size; i++)
	{
	for ( j = 0; j < size; j++)
	{
	for ( k = 0; k < size; k++)
	{
		double Dxn, Dxp, Dyn, Dyp, Dzp, Dzn;
		double aplus, aminus, bplus, bminus, cplus, cminus, dplus, dminus, eplus, eminus, fplus, fminus;
		double gradnorm;
		if (weno_mode <= 2)
		{
			FD_DX(dist, i, j, k, N, H, &Dxp, &Dxn, weno_mode);
			FD_DY(dist, i, j, k, N, H, &Dyp, &Dyn, weno_mode);
			FD_DZ(dist, i, j, k, N, H, &Dzp, &Dzn, weno_mode);
		}
		else
		{
//			printf("Not doing weno now.\n");
//			exit(0);

			WENO_DX(dist, i, j, k, N, H, &Dxp, &Dxn, weno_mode);
			WENO_DY(dist, i, j, k, N, H, &Dyp, &Dyn, weno_mode);
			WENO_DZ(dist, i, j, k, N, H, &Dzp, &Dzn, weno_mode);
		}
		aplus = max(0.0, Dxp);
		aminus = min(0.0, Dxp);
		bplus = max(0.0, Dxn);
		bminus = min(0.0, Dxn);
		cplus = max(0.0, Dyp);
		cminus = min(0.0, Dyp);
		dplus = max(0.0, Dyn);
		dminus = min(0.0, Dyn);
		eplus = max(0.0, Dzp);
		eminus = min(0.0, Dzp);
		fplus = max(0.0, Dzn);
		fminus = min(0.0, Dzn);

		if (dist[i][j][k] > 0.0)
		{
			graddx[i][j][k] = (fabs(aminus)>fabs(bplus))?aminus:bplus;
			graddy[i][j][k] = (fabs(cminus)>fabs(dplus))?cminus:dplus;
			graddz[i][j][k] = (fabs(eminus)>fabs(fplus))?eminus:fplus;
		}
		else
		{
			graddx[i][j][k] = (fabs(aplus)>fabs(bminus))?aplus:bminus;
			graddy[i][j][k] = (fabs(cplus)>fabs(dminus))?cplus:dminus;
			graddz[i][j][k] = (fabs(eplus)>fabs(fminus))?eplus:fminus;
		}
		gradnorm = sqrt(graddx[i][j][k]*graddx[i][j][k] + graddy[i][j][k]*graddy[i][j][k] + \
				graddz[i][j][k]*graddz[i][j][k]);

		if ( (N <= 50) && (fabs(gradnorm - 1.) > 0.1) && ( fabs(dist[i][j][k]) < sqrt(H) ) )
		{
//			graddx[i][j][k] = graddx[i][j][k]/gradnorm;
//			graddy[i][j][k] = graddy[i][j][k]/gradnorm;
//			graddz[i][j][k] = graddz[i][j][k]/gradnorm;
//			if (graddx[i][j][k] != graddx[i][j][k])
//			{
//			
				printf("At (%lf, %lf, %lf), (index %d, %d, %d) dist %lf\n", \
				INT_L+j*H, INT_L+i*H, INT_L+k*H, i, j, k, dist[i][j][k]);

				printf("a (%lf, %lf) b(%lf, %lf) c(%lf, %lf)\nd(%lf, %lf), e(%lf, %lf), f(%lf, %lf)\n", \
				aminus, aplus, bminus, bplus, cminus, cplus, dminus, dplus, eminus, eplus, fminus, fplus);
				printf("grad (%lf, %lf, %lf), gradnorm %lf.\n\n", \
				graddx[i][j][k], graddy[i][j][k], graddz[i][j][k], gradnorm);

//				printf("gradplus (%lf, %lf, %lf), gradminus (%lf, %lf, %lf)\n", max(aminus, bplus), \
//				max(cminus, dplus), max(eminus, fplus), max(aplus, bminus), max(cplus, dminus), \
//				max(eplus, fminus) );
//			}
		}
	}	//	end of k	//
	}	//	end of j	//
	}	//	end of i	//
	return 0;
}



