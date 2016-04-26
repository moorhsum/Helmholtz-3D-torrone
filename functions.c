
#include "functions.h"



/*
double cal_interpgrad(int mode, double xj, double yj, double zj)
{
	double interpolatex, interpolatey, interpolatez;
	double x, y, z;
//	double h = (INT_U-INT_L)/N;
	double temp_hhh, temp_hhl, temp_lhh, temp_lhl, temp_llh, temp_hlh, temp_lll, temp_hll;
	double temp_ll, temp_lh, temp_hl, temp_hh;
	double templ, temph;
	int i, j, k;

	i = (int)((xj-INT_L)*N/(INT_U-INT_L));
	j = (int)((yj-INT_L)*N/(INT_U-INT_L));
	k = (int)((zj-INT_L)*N/(INT_U-INT_L));

	x = INT_L + h*(double)(i);
	y = INT_L + h*(double)(j);
	z = INT_L + h*(double)(k);

//	if (flag < 30)
//	{
//		printf("xj = %lf, yj = %lf, zj = %lf, x = %lf, y = %lf, z = %lf.\n", xj, yj, zj, x, y, z);
//		flag++;
//	}

	if (mode==1)		//	calculating interpolated x coordinate of gradient d	//
	{
		//	Interpolation for x coordinate	//
		temp_lll = gradientdx(x,y,z);
		temp_llh = gradientdx(x,y,z+h);
		temp_lhl = gradientdx(x,y+h,z);
		temp_lhh = gradientdx(x,y+h, z+h);
		temp_hll = gradientdx(x+h,y,z);
		temp_hlh = gradientdx(x+h,y,z+h);
		temp_hhl = gradientdx(x+h,y+h,z);
		temp_hhh = gradientdx(x+h,y+h, z+h);
		temp_ll = temp_lll + (zj-z)*(temp_llh - temp_lll)/h;
		temp_lh = temp_lhl + (zj-z)*(temp_lhh - temp_lhl)/h;
		temp_hl = temp_hll + (zj-z)*(temp_hlh - temp_hll)/h;
		temp_hh = temp_hhl + (zj-z)*(temp_hhh - temp_hhl)/h;
		templ = temp_ll + (yj-y)*(temp_lh - temp_ll)/h;
		temph = temp_hl + (yj-y)*(temp_hh - temp_hl)/h;
		interpolatex = templ + (xj-x)*(temph - templ)/h;

		if ((interpolatex+1==interpolatex)||(interpolatex!=interpolatex))
		{
			printf("it's gradient d, xj = %lf, yj = %lf, zj = %lf.\n", xj, yj, zj);
			printf("x = %lf, y = %lf, z = %lf.\n", x, y, z);
			exit(0);
		}
		return interpolatex;
	}

	else if (mode==2)	//	calculating the y coordinate of gradient d	//
	{
		//	Interpolation for y coordinate	//	
		temp_lll = gradientdy(x,y,z);
		temp_llh = gradientdy(x,y,z+h);
		temp_lhl = gradientdy(x,y+h,z);
		temp_lhh = gradientdy(x,y+h, z+h);
		temp_hll = gradientdy(x+h,y,z);
		temp_hlh = gradientdy(x+h,y,z+h);
		temp_hhl = gradientdy(x+h,y+h,z);
		temp_hhh = gradientdy(x+h,y+h, z+h);
		temp_ll = temp_lll + (zj-z)*(temp_llh - temp_lll)/h;
		temp_lh = temp_lhl + (zj-z)*(temp_lhh - temp_lhl)/h;
		temp_hl = temp_hll + (zj-z)*(temp_hlh - temp_hll)/h;
		temp_hh = temp_hhl + (zj-z)*(temp_hhh - temp_hhl)/h;
		templ = temp_ll + (yj-y)*(temp_lh - temp_ll)/h;
		temph = temp_hl + (yj-y)*(temp_hh - temp_hl)/h;
		interpolatey = templ + (xj-x)*(temph - templ)/h;
		
		if ((interpolatey+1==interpolatey)||(interpolatey!=interpolatey))
		{
			printf("it's gradient d, xj = %lf, yj = %lf, zj = %lf.\n", xj, yj, zj);
			printf("x = %lf, y = %lf, z = %lf.\n", x, y, z);
			exit(0);
		}
		return interpolatey;
	}
	else if (mode==3)	//	calculating the z coordinate of gradient d	//
	{
		//	Interpolation for y coordinate	//	
		temp_lll = gradientdz(x,y,z);
		temp_llh = gradientdz(x,y,z+h);
		temp_lhl = gradientdz(x,y+h,z);
		temp_lhh = gradientdz(x,y+h, z+h);
		temp_hll = gradientdz(x+h,y,z);
		temp_hlh = gradientdz(x+h,y,z+h);
		temp_hhl = gradientdz(x+h,y+h,z);
		temp_hhh = gradientdz(x+h,y+h, z+h);
		temp_ll = temp_lll + (zj-z)*(temp_llh - temp_lll)/h;
		temp_lh = temp_lhl + (zj-z)*(temp_lhh - temp_lhl)/h;
		temp_hl = temp_hll + (zj-z)*(temp_hlh - temp_hll)/h;
		temp_hh = temp_hhl + (zj-z)*(temp_hhh - temp_hhl)/h;
		templ = temp_ll + (yj-y)*(temp_lh - temp_ll)/h;
		temph = temp_hl + (yj-y)*(temp_hh - temp_hl)/h;
		interpolatez = templ + (xj-x)*(temph - templ)/h;
		
		if ((interpolatez+1==interpolatez)||(interpolatez!=interpolatez))
		{
			printf("it's gradient d, xj = %lf, yj = %lf, zj = %lf.\n", xj, yj, zj);
			printf("x = %lf, y = %lf, z = %lf.\n", x, y, z);
			exit(0);
		}
		return interpolatez;
	}
	else
	{
		printf("Error in gradient d mode.\n");
		exit(0);
	}
}

*/
double cal_partial(int mode, double xi, double yi, double zi, double xj, double yj, double zj)
{
	double gradphix, gradphiy, gradphiz;
	static int flag = 0;
	static double tau, threshold;
	
	
	//	Center gradient of the fundamental solution on the second variable	//
//	if (REG_MODE != 2)
//	{
		switch (mode)
		{
			case 1:
				return ( cal_phi(xi, yi, zi, xj+H, yj, zj) - cal_phi(xi, yi, zi, xj-H, yj, zj) )/(2.0*H);
				break;
			case 2:
				return ( cal_phi(xi, yi, zi, xj, yj+H, zj) - cal_phi(xi, yi, zi, xj, yj-H, zj) )/(2.0*H);
				break;
			case 3:
				return ( cal_phi(xi, yi, zi, xj, yj, zj+H) - cal_phi(xi, yi, zi, xj, yj, zj-H) )/(2.0*H);
				break;
			case 4:
				return ( phiimg(xi, yi, zi, xj+H, yj, zj) - phiimg(xi, yi, zi, xj-H, yj, zj) )/(2.0*H);
				break;
			case 5:
				return ( phiimg(xi, yi, zi, xj, yj+H, zj) - phiimg(xi, yi, zi, xj, yj-H, zj) )/(2.0*H);
				break;
			case 6:
				return ( phiimg(xi, yi, zi, xj, yj, zj+H) - phiimg(xi, yi, zi, xj, yj, zj-H) )/(2.0*H);
				break;
			default:
				printf("Error in partial mode.\n");
				exit(0);
		}
//	}
/*	else
	{
		if (flag==0)
		{
			tau = h * TAU_FACTOR;
			threshold = tau * tau;
			flag++;
		}

		if ( pow(xi-xj,2.0)+pow(yi-yj,2.0)+pow(zi-zj,2.0) < threshold )
		{
			switch (mode)
			{
				case 1:
					return ( cal_phi(1,xi, yi, zi, xj+h, yj, zj) - cal_phi(1,xi, yi, zi, xj-h, yj, zj) )/(2.0*h);
					break;
				case 2:
					return ( cal_phi(1,xi, yi, zi, xj, yj+h, zj) - cal_phi(1,xi, yi, zi, xj, yj-h, zj) )/(2.0*h);
					break;
				case 3:
					return ( cal_phi(1,xi, yi, zi, xj, yj, zj+h) - cal_phi(1,xi, yi, zi, xj, yj, zj-h) )/(2.0*h);
					break;
				case 4:
					return ( phiimg(xi, yi, zi, xj+h, yj, zj) - phiimg(xi, yi, zi, xj-h, yj, zj) )/(2.0*h);
					break;
				case 5:
					return ( phiimg(xi, yi, zi, xj, yj+h, zj) - phiimg(xi, yi, zi, xj, yj-h, zj) )/(2.0*h);
					break;
				case 6:
					return ( phiimg(xi, yi, zi, xj, yj, zj+h) - phiimg(xi, yi, zi, xj, yj, zj-h) )/(2.0*h);
					break;
				default:
					printf("Error in partial mode.\n");
					exit(0);
			}
		}
		else
		{
			switch (mode)
			{
				case 1:
					return ( cal_phi(0,xi, yi, zi, xj+h, yj, zj) - cal_phi(0,xi, yi, zi, xj-h, yj, zj) )/(2.0*h);
					break;
				case 2:
					return ( cal_phi(0,xi, yi, zi, xj, yj+h, zj) - cal_phi(0,xi, yi, zi, xj, yj-h, zj) )/(2.0*h);
					break;
				case 3:
					return ( cal_phi(0,xi, yi, zi, xj, yj, zj+h) - cal_phi(0,xi, yi, zi, xj, yj, zj-h) )/(2.0*h);
					break;
				case 4:
					return ( phiimg(xi, yi, zi, xj+h, yj, zj) - phiimg(xi, yi, zi, xj-h, yj, zj) )/(2.0*h);
					break;
				case 5:
					return ( phiimg(xi, yi, zi, xj, yj+h, zj) - phiimg(xi, yi, zi, xj, yj-h, zj) )/(2.0*h);
					break;
				case 6:
					return ( phiimg(xi, yi, zi, xj, yj, zj+h) - phiimg(xi, yi, zi, xj, yj, zj-h) )/(2.0*h);
					break;
				default:
					printf("Error in partial mode.\n");
					exit(0);
			}
		}
	}
*/
}
/*
double cal_regpartial(double xi, double yi, double zi, double xj, double yj, double zj, double tau)
{
	double term1, term2, term3;
	term1 = 1.0/(4.0*PI*tau*radius);			//	1/8pi*tau(2/r)	//
	term2 = (-5.0*tau)/(256.0*pow(radius,3.0)*PI);	//	-(1/pi)*(5/512)*(2/r^3)*tau	//
	term3 = (-50.0*tau)/(1536.0*pow(radius,3.0)*PI);	//	-(1/pi)*(25/1536)*(2/r^3)*tau	//

	if ( sqrt( fabs(pow(xi-xj, 2.0) + pow(yi-yj,2.0) + pow(zi-zj,2.0))) < tau )
		return term1 + term2 + term3;
	else
		return cal_partial(xi, yi, zi, xj, yj, zj);
}
*/
/*
double starx(double zx, double zy, double zz)
{	return zx - distance(zx, zy, zz)*gradientdx(zx, zy, zz); 	}
double stary(double zx, double zy, double zz)
{	return zy - distance(zx, zy, zz)*gradientdy(zx, zy, zz);	}
double starz(double zx, double zy, double zz)
{	return zz - distance(zx, zy, zz)*gradientdz(zx, zy, zz);	}

double gradientdx(double x, double y, double z)
{	return (distance( x+h, y, z) - distance( x-h, y, z))/(2.0*h);	}
double gradientdy(double x, double y, double z)
{	return (distance( x, y+h, z) - distance( x, y-h, z))/(2.0*h);	}
double gradientdz(double x, double y, double z)
{	return (distance( x, y, z+h) - distance( x, y, z-h))/(2.0*h);	}


double secondpartialxx(double x, double y, double z)
{	return (distance( x+h, y, z) - 2.0*distance( x, y, z) + distance( x-h, y, z))/(h*h);	}
double secondpartialyy(double x, double y, double z)
{	return (distance( x, y+h, z) - 2.0*distance( x, y, z) + distance( x, y-h, z))/(h*h);	}
double secondpartialzz(double x, double y, double z)
{	return (distance( x, y, z+h) - 2.0*distance( x, y, z) + distance( x, y, z-h))/(h*h);	}


double secondpartialxy(double x, double y, double z)
{	return (distance(x+h,y+h,z) - distance(x+h,y-h,z) - distance(x-h,y+h,z) + distance(x-h,y-h,z))/(4.0*h*h);	}
double secondpartialxz(double x, double y, double z)
{	return (distance(x+h,y,z+h) - distance(x+h,y,z-h) - distance(x-h,y,z+h) + distance(x-h,y,z-h))/(4.0*h*h);	}
double secondpartialyz(double x, double y, double z)
{	return (distance(x,y+h,z+h) - distance(x,y+h,z-h) - distance(x,y-h,z+h) + distance(x,y-h,z-h))/(4.0*h*h);	}


double laplace(double x, double y, double z)
{
	double Dplusx, Dminusx, Dplusy, Dminusy, Dplusz, Dminusz;
	Dplusx = ( distance(x+h, y, z) - distance(x, y, z) )/h;
	Dminusx = ( distance(x, y, z) - distance(x-h, y, z) )/h;
	Dplusy = ( distance(x, y+h, z) - distance(x, y, z) )/h;
	Dminusy = ( distance(x, y, z) - distance(x, y-h, z) )/h;
	Dplusz = ( distance(x, y, z+h) - distance(x, y, z) )/h;
	Dminusz = ( distance(x, y, z) - distance(x, y, z-h) )/h;
	return ( Dplusx-Dminusx + Dplusy-Dminusy + Dplusz-Dminusz )/h;
}

double GaussianCurv(double x, double y, double z)
{
	int i, j;
	double A[3][3], Adj[3][3];
	double grad[3];
	double second;

	second = 0.0;
	grad[0] = gradientdx(x, y, z);
	grad[1] = gradientdy(x, y, z);
	grad[2] = gradientdz(x, y, z);

	A[0][0] = secondpartialxx(x,y,z);
	A[1][1] = secondpartialyy(x,y,z);
	A[2][2] = secondpartialzz(x,y,z);
	A[0][1] = secondpartialxy(x,y,z);
	A[1][0] = A[0][1];
	A[0][2] = secondpartialxz(x,y,z);
	A[2][0] = A[0][2];
	A[1][2] = secondpartialyz(x,y,z);
	A[2][1] = A[1][2];	
	cal_adj(A, Adj);

	for (i = 0; i < 3; i++)
		for ( j = 0; j < 3; j++)
			second += grad[i]*grad[j]*Adj[i][j];
	return second;
}

double MeanCurv(double x, double y, double z)
{	return -0.5*laplace(x,y,z);	}

double cal_J(int mode, double x, double y, double z)
{
	int i, j;
	double A[3][3], Adj[3][3];
	double grad[3];
	double norm, second = 0.0;
	if (mode==0)
		return 1.0;
	else if (mode==1)
		return 1.0 - distance(x,y,z) * laplace(x,y,z);
	else
	{
		norm = sqrt( pow(x-Cx,2.0) + pow(y-Cy,2.0) + pow(z-Cz, 2.0) );
		grad[0] = gradientdx(x, y, z);
		grad[1] = gradientdy(x, y, z);
		grad[2] = gradientdz(x, y, z);

		//	BUILD THE HESSIAN OF D(z)	//
		if (ONE_SPHERE==1)
		{
			A[0][0] = (-1.) * (pow(y-Cy,2.0) + pow(z-Cz,2.0))/pow(norm,3.0);
			A[1][1] = (-1.) * (pow(x-Cx,2.0) + pow(z-Cz,2.0))/pow(norm,3.0);
			A[2][2] = (-1.) * (pow(x-Cx,2.0) + pow(y-Cy,2.0))/pow(norm,3.0);
			A[0][1] = (x-Cx)*(y-Cy)/pow(norm,3.0);
			A[1][0] = A[0][1];
			A[0][2] = (x-Cx)*(z-Cz)/pow(norm,3.0);
			A[2][0] = A[0][2];
			A[1][2] = (y-Cy)*(z-Cz)/pow(norm,3.0);
			A[2][1] = A[1][2];
		}
		else
		{
			A[0][0] = secondpartialxx(x,y,z);
			A[1][1] = secondpartialyy(x,y,z);
			A[2][2] = secondpartialzz(x,y,z);
			A[0][1] = secondpartialxy(x,y,z);
			A[1][0] = A[0][1];
			A[0][2] = secondpartialxz(x,y,z);
			A[2][0] = A[0][2];
			A[1][2] = secondpartialyz(x,y,z);
			A[2][1] = A[1][2];	
		}
		
		cal_adj(A, Adj);

		for ( i = 0; i < 3; i++)
			for ( j = 0; j < 3; j++)
				second += grad[i]*grad[j]*Adj[i][j];

		return 1.0 - distance(x,y,z) * laplace(x,y,z) + pow(distance(x,y,z),2.0) * second;
	}
}
*/

double cal_delta(double epsilon, double x)
{
	if (strcmp(DELTA_MODE,"cosine")==0)
	{
		if (fabs(x) <= epsilon)
			return ( 1.0/(2.0*epsilon)) * (1.0 + cos(PI*x/epsilon));
		else
			return 0.0;
	}
	else if (strcmp(DELTA_MODE,"hat")==0)
	{
		if (fabs(x) <= epsilon)
			return 1.0/epsilon - fabs(x)/(epsilon*epsilon);
		else
			return 0.0;
	}
	else
		printf("Delta mode error.\n");
	return 0.0;
}

double cal_f(int m, int l , double x)
{
	double front;		//	The (-1)^m * (1-x^2)^(m/2)	factor to multiply in front	//
	if (l==0)
		return 1.0;
	if (m!=0)
		front = pow(-1.0,m) * pow( 1.0-x*x , (double)(m)/2.0);
	if ( (m > l) || (l > 3) || (l < 1))
	{
		printf("There's an error calculating Legendre.\n");
		return 0.0;
	}
	else if (m==l)
	{
		if (m==1)
			return 1.0 * front;
		if (m==2)
			return 3.0 * front;
		if (m==3)
			return 15.0 * front;
	}
	else
	{
		if (l==1)
			return x;				//	P1(x) = x		//
		if (l==2)
		{
			if (m==0)
				return 1.5*pow(x,2.0) - 0.5;	//	P2(x) = 1.5x^2 - 0.5	//
			if (m==1)
				return 3.0*x*front;
		}
		if (l==3)
		{
			if (m==0)
				return 2.5*pow(x,3.0) - 1.5*x;	//	P3(x) = 2.5x^3 - 1.5x	//
			if (m==1)
				return (7.5*pow(x,2.0) - 1.5)*front;
			if (m==2)
				return 15.0*x*front;
		}
		printf("There's error calculating Legendre because of l.\n");
	}
	return 0.0;
}

double cal_phi(double xi, double yi, double zi, double xj, double yj, double zj)
{
	double norm;
	static int flag = 0;
	norm = sqrt(fabs( (xi-xj)*(xi-xj) + (yi-yj)*(yi-yj) + (zi-zj)*(zi-zj) ));


	if (norm == 0.0)
	{
		printf("At (%lf,%lf, %lf) with (%lf, %lf, %lf), norm = 0.\n", xi, yi, zi, xj, yj, zj);
		return 0.0;
	}
///////////////////			DEBUGGING START		//////////////////////
	if (flag < 0)
	{
		printf("norm = %lf, phi = %lf.\n", norm, cos(WAVE*norm)/(-4.0*PI*norm));
		flag++;
	}
//////////////////			DEBUGGING END		//////////////////////

	if (LAPLACE==1)
		return 1.0/(-4.0*PI*norm);
	else
		return cos(WAVE*norm)/(-4.0*PI*norm);

/*
	if (REG_MODE!=2)
		return 1.0/(-4.0*PI*norm);
	if (REG_MODE==2)
	{
		if (LAPLACE==1)
		{
			if (mode==0)
				return 1.0/(-4.0*PI*norm);
			else
				return erf(norm/G_DELTA)/(-4.0*PI*norm);
		}
		else
		{
			if (mode==0)
				return cos(WAVE*norm)/(-4.0*PI*norm);
			else
				return erf(norm/G_DELTA)*cos(WAVE*norm)/(-4.0*PI*norm);	
		}
	}
*/
}

double phiimg(double xi, double yi, double zi, double xj, double yj, double zj)
{
	double norm;
	static int flag = 0;
	norm = sqrt(fabs(pow(xi-xj,2.0) + pow(yi-yj,2.0) + pow(zi-zj,2.0)));

	if (norm == 0.0)
	{
		printf("At (%lf,%lf,%lf) with (%lf,%lf,%lf), norm = 0.\n", xi, yi, zi, xj, yj, zj);
		return 0.0;
	}
///////////////////			DEBUGGING START		//////////////////////
	if (flag < 0)
	{
		printf("norm = %lf, phiimg = %lf.\n", norm, sin(WAVE*norm)/(-4.0*PI*norm) );
		flag++;
	}
//////////////////			DEBUGGING END		//////////////////////

	if (LAPLACE==1)
		return 0.0;
	else
		return sin(WAVE*norm)/(-4.0*PI*norm);
}

double psi(double xj, double yj, double zj)
{
	double norm;
	static double eta = PSI_ETA;
	static double etacube = PSI_ETA * PSI_ETA * PSI_ETA;
	static double neta = (-1.0) * PSI_ETA * PSI_ETA;
	
	norm = fabs( pow(xj-Wavex, 2.0) + pow(yj-Wavey, 2.0) + pow(zj-Wavez, 2.0) );

	if (norm > 0.01)
		return 0.0;
	else
		return exp(norm/neta)/etacube;		//	e^(-(|x|/eta)^2) / eta		//
}

double psiimg(double xj, double yj, double zj)
{
	static int flag = 0;
	if (flag ==0)
	{
		printf("psiimg has not been written yet.\n");
		flag++;
	}
	return 0.0;
}


double cal_nonhom(double zstarx, double zstary, double zstarz, double A[])
{
	int i, j, k;
	static double Waveh = (WAVEXINT_U - WAVEXINT_L)/WaveN;
	static double Wavehcube;
	double wx, wy, wz;
	double realphi, imgphi, realpsi, imgpsi;
	double sum = 0.0;
	double sumimg = 0.0;

	Wavehcube = pow(Waveh, 3.0);

	for ( i = 0; i < WaveN; i++)
	{
		wy = WAVEZINT_L + Waveh * (double) (i);
		for ( j = 0; j < WaveN; j++)
		{
			wx = WAVEYINT_L + Waveh * (double)(j);
			for ( k = 0; k < WaveN; k++)
			{
				wz = WAVEXINT_L + Waveh * (double)(k);
				if ( distance(wx, wy, wz) >= 0)
				{
					sum += 0.0;
				}
				else
				{
					realphi = cal_phi(zstarx, zstary, zstarz, wx, wy, wz);
					imgphi = phiimg(zstarx, zstary, zstarz, wx, wy, wz);
					realpsi = psi(wx, wy, wz);
					imgpsi = psiimg(wx, wy, wz);
					sum += ( realphi*realpsi - imgphi*imgpsi);
					sumimg += (imgphi*realpsi + realphi*imgpsi);
					if ( (sum != sum) || (sum+1.0==sum) )
					{
						printf("Nonhomogeneous part at (%lf,%lf,%lf) with (%lf,%lf,%lf) is %lf.\n", wx, wy, wz, zstarx, zstary, zstarz, sum);
						exit(0);
					}
				}
			}
		}
	}
	sum = sum * Wavehcube; 
	sumimg = sumimg * Wavehcube;
	A[0] = sum;
	A[1] = sumimg;
	return sum;
}


double distance(double x, double y, double z)
{	
	int i, j, k;
	if (DIST_MODE==0)
	{
//		return max(radius-sqrt(fabs((x-Cx)*(x-Cx) + (y-Cy)*(y-Cy) + (z-Cz)*(z-Cz))), \
		   radius2-sqrt(fabs((x-Cx2)*(x-Cx2) + (y-Cy2)*(y-Cy2) + (z-Cz2)*(z-Cz2))));
		return radius - sqrt( (x-Cx)*(x-Cx) + (y-Cy)*(y-Cy) + (z-Cz)*(z-Cz) );	

		return radius - sqrt( (x-Cx)*(x-Cx) + (y-Cy)*(y-Cy) + (z-Cz)*(z-Cz) );
	}
	else
	{
		static int ini = 0;
		static double ***distance;
		static double ***real;
		static double *coordinate;
		static double distance_h;
		static double **distancei, **distancei1;
		static double *distanceij, *distanceij1, *distancei1j, *distancei1j1;
		static double distanceijk;
		static double **reali;
		static double *realij;
		

		int indx, indy, indz;
		double templll, templlh, templhl, templhh, temphll, temphlh, temphhl, temphhh;
		double templl, templh, temphl, temphh, templ, temph;
		double zx, zy, zz;
		double temp;

		static int matrixsize;
		char distfilename[50];

		FILE *fpdist;

//		printf("In distance.\n");

		if (ini==0)
		{
			sprintf(distfilename, "distance%d.txt", N);
			if ((fpdist = fopen(distfilename, "r"))==NULL)
			{
				printf("Open file %s error.\n", distfilename);
				exit(0);
			}
			fscanf(fpdist, "%d", &matrixsize);
			distance_h = (INT_U-INT_L)/matrixsize;
//			double matrix[matrixsize][matrixsize];
			real = (double ***) malloc((int)(SIZE) * sizeof(double **));
			distance = (double ***) malloc(matrixsize * sizeof(double **));
			coordinate = (double *) malloc(matrixsize * sizeof(double));
			
			for ( i = 0; i < SIZE; i++)
			{
				real[i] = (double **) malloc((int)(SIZE) * sizeof(double *));
				if (real[i]==NULL)
				{
					printf("NULL pointer at real[%d].\n", i);
					exit(0);
				}
				for ( j = 0; j < SIZE; j++)
				{
					real[i][j] = (double *) malloc((int)(SIZE) * sizeof(double));
					if (real[i][j]==NULL)
					{
						printf("NULL pointer at real[%d][%d].\n", i, j);
						exit(0);
					}
				}
			}

			for ( i = 0 ; i < matrixsize; i++)
			{
				distance[i] = (double **) malloc(matrixsize * sizeof(double *));
				if (distance[i]==NULL)
				{
					printf("NULL pointer at distance[%d].\n", i);
					exit(0);
				}
				for ( j = 0; j < matrixsize; j++)
				{
					distance[i][j] = (double *) malloc(matrixsize * sizeof(double));
					if (distance[i][j]==NULL)
					{
						printf("NULL pointer at distance[%d][%d].\n", i, j);
						exit(0);
					}
				}
				coordinate[i] = INT_L + i * distance_h;
			}

			for ( i = 0; i < matrixsize; i++)
			{
				for ( j = 0; j < matrixsize; j++)
				{
					for ( k = 0; k < matrixsize; k++)
					{
						fscanf(fpdist, "%lf", &distance[i][j][k]);
					}
				}
			}

			if (matrixsize < SIZE)
			{
				for ( i = 0; i < SIZE; i++)
				{
					zy = INT_L + H * i;
					indy = (int) ( (zy-INT_L)*matrixsize/(INT_U-INT_L) + ROUNDING );

					distancei = distance[indy];
					if (indz == matrixsize-1)
					{
						distancei1 = distance[indy];
					}
					else
					{
						distancei1 = distance[indy+1];
					}
					reali = real[i];
					

					for ( j = 0; j < SIZE; j++)
					{
						zx = INT_L + H * j;
						indx = (int) ( (zx-INT_L)*matrixsize/(INT_U-INT_L) + ROUNDING );

						distanceij = distancei[indx];
						distancei1j = distancei1[indx];
						if (indy == matrixsize-1)
						{
							distanceij1 = distancei[indx];
							distancei1j1 = distancei1[indx];
						}
						else
						{
							distanceij1 = distancei[indx+1];
							distancei1j1 = distancei1[indx+1];
						}
						realij = reali[j];

						for ( k = 0; k < SIZE; k++)
						{
							zz = INT_L + H * k;
							indz = (int) ( (zz-INT_L)*matrixsize/(INT_U-INT_L) + ROUNDING );

//							printf("ini = %d,  i = %d, j = %d, indx = %d, indy = %d, zx = %lf\n", ini, i, j, indx, indy, zx);

							templll = distanceij[indz];
							templhl = distanceij1[indz];
							temphll = distancei1j[indz];
							temphhl = distancei1j1[indz];
							if (indz == matrixsize-1)
							{
								templlh = distanceij[indz];
								templhh = distanceij1[indz];
								temphlh = distancei1j[indz];
								temphhh = distancei1j1[indz];
							}
							else
							{
								templlh = distanceij[indz+1];
								templhh = distanceij1[indz+1];
								temphlh = distancei1j[indz+1];
								temphhh = distancei1j1[indz+1];
							}
							

							templl = templll + (templlh - templll) * (zx - coordinate[indx])/ distance_h;
							templh = templhl + (templhh - templhl) * (zx - coordinate[indx])/ distance_h;
							temphl = temphll + (temphlh - temphll) * (zx - coordinate[indx])/ distance_h;
							temphh = temphhl + (temphhh - temphhl) * (zx - coordinate[indx])/ distance_h;

	   					        templ = templl + (templh - templl) * (zy - coordinate[indy]) / distance_h;
						        temph = temphl + (temphh - temphl) * (zy - coordinate[indy]) / distance_h;
						        realij[k] = templ + (temph - templ) * (zz - coordinate[indz])/ distance_h;
//							printf("here. indx = %d, indy = %d\n", indx, indy);
						}
					}
				}
			}
			ini++;
		}

//		printf("here. x = %lf, y = %lf, z = %lf\n", x, y, z);
		if ( matrixsize >= SIZE )
		{
			if (ini==1)
			{
				printf("here.\n");
				ini++;
			}
			indx = (int) ( (x-INT_L)*matrixsize/(INT_U-INT_L) + ROUNDING );
			indy = (int) ( (y-INT_L)*matrixsize/(INT_U-INT_L) + ROUNDING );
			indz = (int) ( (z-INT_L)*matrixsize/(INT_U-INT_L) + ROUNDING );

//		printf("matrixsize = %d, indx = %d, indy = %d, indz = %d.\n", matrixsize, indx, indy, indz);
			return -1.0 * distance[indz][indy][indx];
		}
		else
		{
			indx = (int) ( (x-INT_L)*N/(INT_U-INT_L) + ROUNDING );
			indy = (int) ( (y-INT_L)*N/(INT_U-INT_L) + ROUNDING );
			indz = (int) ( (z-INT_L)*N/(INT_U-INT_L) + ROUNDING );
			return -1.0 * real[indz][indy][indx];
		}
	}
}


double PolyWeight(double t, double tau)
{
	double term1, term2, term3, tau2, tau3, tau4;
	if ( (t <= tau) || (t >= 1.0) )
		return 0.0;
	else
	{
		tau2 = tau * tau;
//		tau3 = tau2 * tau;
//		tau4 = tau2 * tau2;
//		return 2./(1.-tau) - 4.*fabs(t-0.5*(1.+tau))/((1.-tau)*(1.-tau));
//		return 6. * (t-1.)*(t-tau)/((tau-1.)*(tau-1.)*(tau-1.));
//		return -30.*(t-1.)*(t-1.)*(t-tau)*(t-tau)/pow(tau-1,5.0);

//		tau2 = tau * tau;
		if (0)	//	SECOND MOMENT WITH SMOOTH ENDS	//
		{
			tau3 = tau2 * tau;
			tau4 = tau2 * tau2;

			term1 = -6. * (3. + 8.*tau + 3.* tau2)*t*t;
			term2 = 4. * (5. + 16.*tau + 16. * tau2 + 5. * tau3)*t;
			term3 = -1.* (5. + 20.*tau + 34.*tau2 + 20.*tau3 + 5.*tau4);

			return 210. * (t-1.)*(t-1.)*(t-tau)*(t-tau)*(term1 + term2 + term3)/pow(tau-1.,9.0);
		}
		else	//	FIRST MOMENT WITH SMOOTH ENDS	//
		{
			term1 = 7. * (1+tau)*t;
			term2 = -2. * (2. + 3. * tau + 2. * tau2);

			return 60. * (t-1.)*(t-1.)*(t-tau)*(t-tau)*(term1 + term2)/pow(tau-1., 7.0);
		}
//		SECOND MOMENT CONTINUOUS END	//
//		term1 = 420. * ( 1. + 3.*tau + tau2 ) * t * t;
//		term2 = -60. * ( 8. + 27.*tau + 27.*tau2 + 8.*tau3 ) * t;
//		term3 = 60. * ( 2. + 8.*tau + 15.*tau2 + 8.*tau3 + 2.*tau4 );
//		return ( t - 1. ) * ( t - tau ) * ( term1 + term2 + term3 ) / pow(tau-1., 7.0);
	}
}


double SineWeight(double t, double tau)
{
	double term1, term2, term3, tau2, PI2, PI3;
	double denom;

	if ( (t <= tau) || (t >= 1.0) )
		return 0.0;
	if (MOMENT_ORDER == 2)	//	SECOND MOMENT	//
	{
//		tau2 = tau * tau;
//		PI2 = PI * PI;
//		PI3 = PI2 * PI;
//
//		term1 = 2.*(20. - 40.*tau + 9.*PI2*tau + 20.*tau2) * (1. - cos(2.*PI*(t-tau)/(1.-tau))) / (3.*PI2 - 31.);
//		term2 = 9.*(PI2 + PI2 * tau) * (tau-1.)/32. * (cos(PI*(t-tau)/(1.-tau)) - cos(3.*PI*(t-tau)/(1.-tau)));
//		term3 = -9. * (3.*PI + PI3 - 6.*PI*tau + 4.*PI3*tau + 3.*PI*tau2 + PI3*tau2 ) * \
//			(3.*sin(PI*(t-tau)/(1.-tau)) - sin(3.*PI*(t-tau)/(1.-tau))) / ( 16. * (3. * PI2 - 31.) );
//
//		return (term1 + term2 + term3) / ( (tau-1.)*(tau-1.)*(tau-1.) );
///////////////////////////	TESTING EXPONENTIAL KERNEL APR 2016
		return exp( 0.5/( (t-0.1)*(t-1.) ) ) * ( CONSTA2*t*t + CONSTB2*t + CONSTC2 );
	}
	else	//	FIRST MOMENT	//
	{
//		denom = tau - 1.;
//		term1 = PI * pow( cos( PI * (t-1.) /(2.*denom) ), 3.0);
//		term2 = 7. + tau + 3. * (5. + 3. * tau) * cos(PI * (tau-t)/denom);
//		term3 = sin( PI * (t-1.)/(2. * denom))/(denom*denom);
//		return term1 * term2 * term3;
		return exp( 0.5/( (t-0.1)*(t-1.) ) ) * ( CONSTA1*t + CONSTB1 );
	}
}










