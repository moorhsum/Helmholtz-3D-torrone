
#include "tools.h"


double inner(int size, double vectorA[], double vectorB[])
{
	int i;
	double sum = 0.0;
	for (i = 0; i < size; i++)
	{
		sum += vectorA[i] * vectorB[i];
		if (sum!=sum)
		{
			printf("i = %d, ai = %lf, bi = %lf.\n", i, vectorA[i], vectorB[i]);
			exit(0);
		}
	}
	return sum;
}

double omp_inner(int size, double vectorA[], double vectorB[], int chunk)
{
	int i;
	double sum = 0.0;
	#pragma omp parallel for private(i) schedule(static, chunk) reduction(+:sum)
	for (i = 0; i < size; i++)
	{
		sum += vectorA[i] * vectorB[i];
	}
	return sum;
}

double norm(int size, double vectorA[])
{	return sqrt(fabs(inner(size, vectorA, vectorA)));	}

double BesselJ0(double x)
{
	int i, j;
	static const double polycoef[7] = J0COEFFICIENT;
	static const double fcoef[6] = F0COEFFICIENT;
	static const double thetacoef[6] = THETA0COEFFICIENT;

	double sum = 0.0;

	if (x < 3.0)
	{
		if (x <= 1.0)
		{
			for (i = 0; i < 7; i++)
			{
				sum += (pow(x/3.0, 2.*i) * polycoef[i]);
			}
			return sum;
		}
		else
		{
//			double sum = 0.0;
			for (i = 0; i < 7; i++)
			{
				sum += (pow(x/3.0, 2.*i) * polycoef[i]);
			}
			return sum;
//		for (int k = 0; k < BESSELJ_APROX; k++)
//		{
//			sum += pow(-1.0, k) * pow( pow(x/2.0, k)/factorial((double)k) , 2.0);
//		}
//		return sum;
		}
	}
	else
	{
		double f = 0.0;
		double theta = 0.0;
		if (x < 5.0)
		{
			for (j = 0; j < 6; j++)
			{
				f += (pow(3.0/x, 2*j) * fcoef[j]);
				theta += (pow(3.0/x, 2*j+1) * thetacoef[j]);
			}
			theta += (x - 0.25*PI);
			return f * cos(theta)/sqrt(x);
		}
		else if (x < 12.0)
		{
			for (j = 0; j < 6; j++)
			{
				f += (pow(3.0/x, 2*j) * fcoef[j]);
				theta += (pow(3.0/x, 2*j+1) * thetacoef[j]);
			}
			theta += (x - 0.25*PI);
			return f * cos(theta)/sqrt(x);
			
		}
		else
		{
			for ( j = 0; j < 6; j++)
			{
				f += (pow(3.0/x, 2*j) * fcoef[j]);
				theta += (pow(3.0/x, 2*j+1) * thetacoef[j]);
			}
			theta += (x - 0.25*PI);
			return f * cos(theta)/sqrt(x);
			
		}
	}

}

double BesselY0(double x)
{
	int i, j;
	static const double polycoef[7] = Y0COEFFICIENT;
	static const double fcoef[6] = F0COEFFICIENT;
	static const double thetacoef[6] = THETA0COEFFICIENT;

	double sum = 0.0;
//	double sum2 = EULER;
//	double sum3 = 0.0;
//	double tempJ0;
//	double result = 0.0;
//	int k = 0;
	
	if (x < 3.0)
	{
		if (x <= 0.8)
		{
			for ( i = 0; i < 7; i++)
			{
				sum += (pow(x/3.0, 2*i) * polycoef[i]);
			}
			return sum + (2.0/PI) * log(x/2.0) * BesselJ0(x);
		}
		else
		{
			for ( i = 0; i < 7; i++)
			{
				sum += (pow(x/3.0, 2*i) * polycoef[i]);
			}
			return sum + (2.0/PI) * log(x/2.0) * BesselJ0(x);
		}
	}
	else
	{
		double f = 0.0;
		double theta = 0.0;
		if (x < 5.0)
		{
			for ( j = 0; j < 6; j++)
			{
				f += (pow(3.0/x, 2*j) * fcoef[j]);
				theta += (pow(3.0/x, 2*j+1) * thetacoef[j]);
			}
			theta += (x - 0.25*PI);
			return f * sin(theta)/sqrt(x);
		}
		else if (x < 12.0)
		{
			for ( j = 0; j < 6; j++)
			{
				f += (pow(3.0/x, 2*j) * fcoef[j]);
				theta += (pow(3.0/x, 2*j+1) * thetacoef[j]);
			}
			theta += (x - 0.25*PI);
			return f * sin(theta)/sqrt(x);
		}
		else
		{
			for ( j = 0; j < 6; j++)
			{
				f += (pow(3.0/x, 2*j) * fcoef[j]);
				theta += (pow(3.0/x, 2*j+1) * thetacoef[j]);
			}
			theta += (x - 0.25*PI);
			return f * sin(theta)/sqrt(x);
		}
	}

//	tempJ0 = (log( fabs(x)/2.0 ) + sum2) * BesselJ0(x);

//	for (int k = 1; k < BESSELY_APROX; k++)
//	{
//		sum3 += 1.0/k;
//		sum += pow(-1.0, k+1) * pow( pow(x/2.0, k)/factorial(k), 2.0) * sum3;
//		sum2 += 1.0/k;
//		sum += pow(-1.0, k) * pow( pow(x/2.0, k) / factorial((double) k), 2.0) * (log(fabs(x)/2.0)-sum2);
//		printf("k = %d, sum = %lf, sum2 = %lf.\n", k, sum, sum2);
//		printf("log(x/2.0) = %lf", log(fabs(x)/2.0));
//		printf("preresult = %lf, without multiply = %lf.\n", pow(pow(x/2.0,k)/factorial(k),2.0), sum + tempJ0);
//		printf("k = %d, harmonic = %lf, Current result = %lf.\n", k, sum3, (sum+tempJ0)*(2.0/PI) );
//	}
//	return sum;
//	return (sum + tempJ0)*(2.0/PI);
}

double BesselJ1(double x)
{
	int i, j;
	static const double polycoef[7] = J1COEFFICIENT;
	static const double fcoef[6] = F1COEFFICIENT;
	static const double thetacoef[6] = THETA1COEFFICIENT;

	double sum = 0.0;
	if (x <= 3.0)
	{
		for ( i = 0; i < 7; i++) 
			sum += (pow(x/3.0, 2*i) * polycoef[i]);
		return sum * x;
	}
	else
	{
		double f1 = 0.0;
		double theta1 = 0.0;

		for ( j = 0; j < 6; j++)
		{
			f1 += ((pow(3.0/x, 2*j) * fcoef[j]));
			theta1 += ((pow(3.0/x, 2*j + 1) * thetacoef[j]));
		}
		theta1 += (x - 0.75*PI);
		return f1 * cos(theta1)/sqrt(x);
	}
}

double BesselY1(double x)
{
	int i, j;
	static const double polycoef[7] = Y1COEFFICIENT;
	static const double fcoef[6] = F1COEFFICIENT;
	static const double thetacoef[6] = THETA1COEFFICIENT;
	double sum = 0.0;

	if (x < 3.0)
	{
		sum = (2.0/PI)*(log(x/2.0)*BesselJ1(x)-1.0/x);
		for ( i = 0; i < 7; i++)
			sum += (pow(x/3.0, 2*i+1) * polycoef[i]);
		return sum;
	}
	else
	{
		double f1 = 0.0;
		double theta1 = 0.0;
		for ( j = 0; j < 6; j++)
		{
			f1 += (pow(3.0/x, 2*j) * fcoef[j]);
			theta1 += (pow(3.0/x, 2*j+1) * thetacoef[j]);
		}
		theta1 += (x - 0.75*PI);
		return f1 * sin(theta1)/sqrt(x);
	}
}

double besselj(int order, double x)
{
	switch (order)
	{
		case 0:
			return sin(x)/x;
		case 1:
			return sin(x)/(x*x) - cos(x)/x;
		case 2:
			return (3.0/(x*x)-1.0)*sin(x)/x - 3.0*cos(x)/(x*x);
		case 3:
			return (15.0/(x*x*x) - 6.0/x)*sin(x)/x - (15.0/(x*x) - 1.0)*cos(x)/x;
	}
	printf("Undefined behaviour in besselj when order = %d.\n", order);
	exit(0);
}

double bessely(int order, double x)
{
	switch (order)
	{
		case 0:
			return -1.0*cos(x)/x;
		case 1:
			return -1.0*cos(x)/(x*x) - sin(x)/x;
		case 2:
			return (-3.0/(x*x)+1.0)*cos(x)/x - 3.0*sin(x)/(x*x);
		case 3:
			return (-15.0/(x*x*x) + 6.0/x)*cos(x)/x - (15.0/(x*x) - 1.0)*sin(x)/x;
	}
	printf("Undefined behaviour in besselj when order = %d.\n", order);
	exit(0);
}

double factorial(double n)
{
	if (n <= 1.0)
		return 1.0;
	else
		return n * factorial(n-1.0);
}

double max(double a, double b)
{
	if ( a >= b )
		return a;
	else
		return b;
}

double min(double a, double b)
{
	if ( a <= b )
		return a;
	else
		return b;
}

double max3(double a, double b, double c)
{
	double p;
	p = (a>b)?a:b;
	return (p>c)?p:c;
}

double min3(double a, double b, double c)
{
	double p;
	p = (a<b)?a:b;
	return (p<c)?p:c;
}

double minmod(double x, double y)
{
	if (x * y < 0.0)
		return 0.0;
	else
		return (fabs(x)<fabs(y))?x:y;
}

double old_minmod(double x, double y)
{
	return (fabs(x)<fabs(y))?x:y;
}

double min_abs(double a, double b, double c)
{
	double p;
	p = (fabs(a)<fabs(b))?a:b;
	return (fabs(p)<fabs(c))?p:c;
}

double erf(double x)
{
    /* erf(z) = 2/sqrt(pi) * Integral(0..x) exp( -t^2) dt
    erf(0.01) = 0.0112834772 erf(3.7) = 0.9999998325
    Abramowitz/Stegun: p299, |erf(z)-erf| <= 1.5*10^(-7)
    */
    double y = 1.0 / ( 1.0 + 0.3275911 * x);

    return 1.0 - (((((
    + 1.061405429 * y
    - 1.453152027) * y
    + 1.421413741) * y
    - 0.284496736) * y
    + 0.254829592) * y)
    * exp (-x * x);

}


int cal_adj(double A[][3], double Adj[][3])
{
	int i, j;
//	int tmpi, tmpj;
	for (i = 0; i < 3; i++)
		for ( j = 0; j < 3; j++)	
		{
	Adj[j][i] = A[(i+1)%3][(j+1)%3] * A[(i+2)%3][(j+2)%3] - A[(i+2)%3][(j+1)%3] * A[(i+1)%3][(j+2)%3];
		}
	return 0;
}



