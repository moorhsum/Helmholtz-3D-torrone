
#include "finite_difference3D.h"


double Godunov(double signage, double a, double b, double c, double d, double e, double f)
{
	double aplus, aminus, bplus, bminus, cplus, cminus, dplus, dminus, eplus, eminus, fplus, fminus;
	aplus = max(a, 0.0) * max(a, 0.0);
	aminus = min(a, 0.0)* min(a, 0.0);
	bplus = max(b, 0.0) * max(b, 0.0);
	bminus = min(b, 0.0)* min(b, 0.0);
	cplus = max(c, 0.0) * max(c, 0.0);
	cminus = min(c, 0.0)* min(c, 0.0);
	dplus = max(d, 0.0) * max(d, 0.0);
	dminus = min(d, 0.0)* min(d, 0.0);
	eplus = max(e, 0.0) * max(e, 0.0);
	eminus = min(e, 0.0)* min(e, 0.0);
	fplus = max(f, 0.0) * max(f, 0.0);
	fminus = min(f, 0.0)* min(f, 0.0);

	if (signage > 0)
		return sqrt( max(aminus, bplus) + max(cminus, dplus) + max(eminus, fplus) );
	else if (signage < 0)
		return sqrt( max(aplus, bminus) + max(cplus, dminus) + max(eplus, fminus) );
	else
		return 0.0;
}

int FD_DX(double ***phi, int i, int j, int k, const int N, double H, double *Dxp, double *Dxn, int mode)
{
	double **phii = phi[i];
	if (mode >= 2)
		mode = 2;
	if (mode == 1)
	{
		if (j == 0)
		{
			*Dxp = (phii[j+1][k] - phii[j][k])/H;
			*Dxn = -1.0 * (*Dxp);
		}
		else if (j == N)
		{
			*Dxn = (phii[j][k] - phii[j-1][k])/H;
			*Dxp = -1.0 * (*Dxn);
		}
		else
		{
			*Dxp = (phii[j+1][k] - phii[j][k])/H;
			*Dxn = (phii[j][k] - phii[j-1][k])/H;
		}
	}
	else if (mode == 2)
	{
		double DDxp, DDxc, DDxn;
		if (j == 0)
		{
			*Dxp = (phii[j+1][k] - phii[j][k])/H;
			*Dxn = -1.0 * (*Dxp);
		}
		else if ((j == 1) || (j == N-1))
		{
			*Dxp = (phii[j+1][k] - phii[j][k])/H;
			*Dxn = (phii[j][k] - phii[j-1][k])/H;
		}
		else if (j == N)
		{
			*Dxn = (phii[j][k] - phii[j-1][k])/H;
			*Dxp = -1.0 * (*Dxn);
		}
		else
		{
			DDxp = XDER2(phi, i, j+1, k, N, H);
			DDxc = XDER2(phi, i, j, k, N, H);
			DDxn = XDER2(phi, i, j-1, k, N, H);
			*Dxp = (phii[j+1][k] - phii[j][k])/H - 0.5 * H * minmod(DDxp, DDxc);
			*Dxn = (phii[j][k] - phii[j-1][k])/H + 0.5 * H * minmod(DDxn, DDxc);
		}
		return 0;
	}
	printf("Weird in FD_DX. mode = %d\n", mode);
	exit(0);
//	if (mode )

}

int FD_DY(double ***phi, int i, int j, int k, const int N, double H, double *Dyp, double *Dyn, int mode)
{
	if (mode >= 2)
		mode = 2;
	if (mode == 1)
	{
		if (i == 0)
		{
			*Dyp = (phi[i+1][j][k] - phi[i][j][k])/H;
			*Dyn = -1.0 * (*Dyp);
		}
		else if (i == N)
		{
			*Dyn = (phi[i][j][k] - phi[i-1][j][k])/H;
			*Dyp = -1.0 * (*Dyn);
		}
		else
		{
			*Dyp = (phi[i+1][j][k] - phi[i][j][k])/H;
			*Dyn = (phi[i][j][k] - phi[i-1][j][k])/H;
		}
	}
	else if (mode == 2)
	{
		double DDyp, DDyn, DDyc;
		if (i == 0)
		{
			*Dyp = (phi[i+1][j][k] - phi[i][j][k])/H;
			*Dyn = -1.0 * (*Dyp);
		}
		else if ( (i == 1)|| (i == N-1) )
		{
			*Dyp = (phi[i+1][j][k] - phi[i][j][k])/H;
			*Dyn = (phi[i][j][k] - phi[i-1][j][k])/H;
		}
		else if (i == N)
		{
			*Dyn = (phi[i][j][k] - phi[i-1][j][k])/H;
			*Dyp = -1.0 * (*Dyn);
		}
		else
		{
			DDyp = YDER2(phi, i+1, j, k, N, H);
			DDyc = YDER2(phi, i, j, k, N, H);
			DDyn = YDER2(phi, i-1, j, k, N, H);
			*Dyp = (phi[i+1][j][k] - phi[i][j][k])/H - 0.5 * H * minmod(DDyp, DDyc);
			*Dyn = (phi[i][j][k] - phi[i-1][j][k])/H + 0.5 * H * minmod(DDyn, DDyc);
		}
		return 0;
	}
	printf("Weird in FD_DY. mode = %d\n", mode);
	exit(0);
}


int FD_DZ(double ***phi, int i, int j, int k, const int N, double H, double *Dzp, double *Dzn, int mode)
{
	if (mode >= 2)
		mode = 2;
	if (mode == 1)
	{
		if (k == 0)
		{
			*Dzp = (phi[i][j][k+1] - phi[i][j][k])/H;
			*Dzn = -1.0 * (*Dzp);
		}
		else if (k == N)
		{
			*Dzn = (phi[i][j][k] - phi[i][j][k-1])/H;
			*Dzp = -1.0 * (*Dzn);
		}
		else
		{
			*Dzp = (phi[i][j][k+1] - phi[i][j][k])/H;
			*Dzn = (phi[i][j][k] - phi[i][j][k-1])/H;
		}
	}
	else if (mode == 2)
	{
		double DDzp, DDzn, DDzc;
		if (k == 0)
		{
			*Dzp = (phi[i][j][k+1] - phi[i][j][k])/H;
			*Dzn = -1.0 * (*Dzp);
		}
		else if ( (k == 1)|| (k == N-1) )
		{
			*Dzp = (phi[i][j][k+1] - phi[i][j][k])/H;
			*Dzn = (phi[i][j][k] - phi[i][j][k-1])/H;
		}
		else if (k == N)
		{
			*Dzn = (phi[i][j][k] - phi[i][j][k-1])/H;
			*Dzp = -1.0 * (*Dzn);
		}
		else
		{
			DDzp = ZDER2(phi, i, j, k+1, N, H);
			DDzc = ZDER2(phi, i, j, k, N, H);
			DDzn = ZDER2(phi, i, j, k-1, N, H);
			*Dzp = (phi[i][j][k+1] - phi[i][j][k])/H - 0.5 * H * minmod(DDzp, DDzc);
			*Dzn = (phi[i][j][k] - phi[i][j][k-1])/H + 0.5 * H * minmod(DDzn, DDzc);
		}
		return 0;
	}
	printf("Weird in FD_DX. mode = %d\n", mode);
	exit(0);
}


double XDER_C(double ***level, int i, int j, int k, int N, double H)
{
	if ((j == 0)||(j == N))
		return 0.0;
	else
		return (level[i][j+1][k] - level[i][j-1][k])/(2.0*H);
}

double YDER_C(double ***level, int i, int j, int k, int N, double H)
{
	if ((i == 0)|| (i == N))
		return 0.0;
	else
		return (level[i+1][j][k] - level[i-1][j][k])/(2.0*H);
}

double ZDER_C(double ***level, int i, int j, int k, int N, double H)
{
	if ((k == 0)|| (k == N))
		return 0.0;
	else
		return (level[i][j][k+1] - level[i][j][k-1])/(2.0*H);
}

double XDER2(double ***level, int i, int j, int k, int N, double H)
{
	static double Dxpp, Dxp, Dxn, Dxnn;
	static double DDxp, DDxc, DDxn;
	double **temp;

	temp = level[i];

	if (j < 0)
		return (level[i][abs(j)+1][k] - 2.0 * level[i][abs(j)][k] + level[i][abs(j)-1][k])/(H*H);
	else if (j == 0)
		return 2.0*(level[i][j+1][k] - level[i][j][k])/(H*H);
	else if (j == N)
		return 2.0*(level[i][j-1][k] - level[i][j][k])/(H*H);
	else if (j > N)
		return (level[i][2*N-j+1][k] - 2.0 * level[i][2*N-j][k] + level[i][2*N-j-1][k])/(H*H);
	else if (j < 3)
	{
		DDxp = (2.0*level[i][j][k] - 5.0*level[i][j+1][k]+4.0*level[i][j+2][k]-level[i][j+3][k])/(H*H);
		DDxc = (level[i][j+1][k] - 2.0*level[i][j][k] + level[i][j-1][k])/(H*H);
		return (fabs(DDxp)<fabs(DDxc))?DDxp:DDxc;
	}
	else if (j > N - 3)
	{
		DDxn = (2.0*level[i][j][k] - 5.0*level[i][j-1][k]+4.0*level[i][j-2][k]-level[i][j-3][k])/(H*H);
		DDxc = (level[i][j+1][k] - 2.0*level[i][j][k] + level[i][j-1][k])/(H*H);
		return (fabs(DDxn)<fabs(DDxc))?DDxn:DDxc;
	}
	else
	{
		Dxpp = level[i][j+2][k] - level[i][j+1][k];
		Dxp = level[i][j+1][k] - level[i][j][k];
		Dxn = level[i][j][k] - level[i][j-1][k];
		Dxnn = level[i][j-1][k] - level[i][j-2][k];

		DDxp = (2.0*level[i][j][k] - 5.0*level[i][j+1][k] + 4.0*level[i][j+2][k] - level[i][j+3][k])/(H*H);
		DDxc = (level[i][j+1][k] - 2.0*level[i][j][k] + level[i][j-1][k])/(H*H);
		DDxn = (2.0*level[i][j][k] - 5.0*level[i][j-1][k] + 4.0*level[i][j-2][k] - level[i][j-3][k])/(H*H);

//		if ((DDxp!=DDxp)||(DDxc!=DDxc)||(DDxn!=DDxn))
//		{
//			printf("XDER2 at point (%d, %d, %d), (%lf,%lf,%lf).\n", j, i, k, INT_L+double(j)*H, \
//			INT_L+double(i)*H, INT_L+double(k)*H);
//			printf("distance = %lf, DDxp %lf, DDxc %lf, DDxn %lf.\n", temp[j][k].dist, DDxp, DDxc, DDxn)
//			exit(0);
//		}
		if ((Dxp * Dxnn <= 0.0 ) && (Dxpp * Dxn <= 0.0))
//			return old_minmod(DDxp, DDxn);
			return min_abs(DDxp, DDxc, DDxn);
		else if (Dxp * Dxnn <= 0.0)
			return old_minmod(DDxp, DDxc);
		else if (Dxpp * Dxn <= 0.0)
			return old_minmod(DDxn, DDxc);
		else
			return DDxc;
	}
}

double YDER2(double ***level, int i, int j, int k, int N, double H)
{
	static double DDyp, DDyc, DDyn;
	static double Dypp, Dyp, Dyn, Dynn;

	if (i < 0)	//	Periodic condition, mirror effect	//
		return (level[abs(i)+1][j][k] - 2.0 * level[abs(i)][j][k] + level[abs(i)-1][j][k])/(H*H);
	else if (i == 0)
		return 2.0*(level[i+1][j][k] - level[i][j][k])/(H*H);
	else if (i == N)
		return 2.0*(level[i-1][j][k] - level[i][j][k])/(H*H);
	else if (i > N)
		return (level[2*N-i+1][j][k] - 2.0 * level[2*N-i][j][k] + level[2*N-i-1][j][k])/(H*H);
	else if (i < 3)
	{
		DDyp = (2.0*level[i][j][k] - 5.0*level[i+1][j][k] + 4.0*level[i+2][j][k] - level[i+3][j][k])/(H*H);
		DDyc = (level[i+1][j][k] - 2.0*level[i][j][k] + level[i-1][j][k])/(H*H);
		return (fabs(DDyp)<fabs(DDyc))?DDyp:DDyc;
	}
	else if (i > N - 3)
	{
		DDyn = (2.0*level[i][j][k] - 5.0*level[i-1][j][k] + 4.0*level[i-2][j][k] - level[i-3][j][k])/(H*H);
		DDyc = (level[i+1][j][k] - 2.0*level[i][j][k] + level[i-1][j][k])/(H*H);
		return (fabs(DDyn)<fabs(DDyc))?DDyn:DDyc;
	}
	else
	{
		Dypp = level[i+2][j][k] - level[i+1][j][k];
		Dyp = level[i+1][j][k] - level[i][j][k];
		Dyn = level[i][j][k] - level[i-1][j][k];
		Dynn = level[i-1][j][k] - level[i-2][j][k];

		DDyp = (2.0*level[i][j][k] - 5.0*level[i+1][j][k] + 4.0*level[i+2][j][k] - level[i+3][j][k])/(H*H);
		DDyc = (level[i+1][j][k] - 2.0*level[i][j][k] + level[i-1][j][k])/(H*H);
		DDyn = (2.0*level[i][j][k] - 5.0*level[i-1][j][k] + 4.0*level[i-2][j][k] - level[i-3][j][k])/(H*H);

//		if ((DDyp!=DDyp)||(DDyc!=DDyc)||(DDyn!=DDyn))
//		{
//			printf("YDER2 at point (%d, %d, %d), (%lf,%lf,%lf).\n", j, i, k, INT_L+double(j)*H, \
//			INT_L+double(i)*H, INT_L+double(k)*H);
//			printf("distance = %lf, DDyp %lf, DDyc %lf, DDyn %lf.\n", ROOT[i][j][k].dist, DDyp, DDyc, DDyn);
//			exit(0);
//		}
		if ((Dyp * Dynn <= 0.0 ) && (Dypp * Dyn <= 0.0))
//			return old_minmod(DDyp, DDyn);
			return min_abs(DDyp, DDyn, DDyc);
		else if (Dyp * Dynn <= 0.0)
			return old_minmod(DDyp, DDyc);
		else if (Dypp * Dyn <= 0.0)
			return old_minmod(DDyn, DDyc);
		else
			return DDyc;
	}
}


double ZDER2(double ***level, int i, int j, int k, int N, double H)
{
	static double DDzp, DDzc, DDzn;
	static double Dzpp, Dzp, Dzn, Dznn;
	double *temp;
	
	temp = level[i][j];

	if (k < 0)	//	Periodic condition, mirror effect	//
		return (level[i][j][abs(k)+1] - 2.0 * level[i][j][abs(k)] + level[i][j][abs(k)-1])/(H*H);
	else if (k == 0)
		return 2.0*(level[i][j][k+1] - level[i][j][k])/(H*H);
	else if (k == N)
		return 2.0*(level[i][j][k-1] - level[i][j][k])/(H*H);
	else if (k > N)
		return (level[i][j][2*N-k+1] - 2.0 * level[i][j][2*N-k] + level[i][j][2*N-k-1])/(H*H);
	else if (k < 3)
	{
		DDzp = (2.0*temp[k] - 5.0*temp[k+1] + 4.0*temp[k+2] - temp[k+3])/(H*H);
		DDzc = (temp[k+1] - 2.0*temp[k] + temp[k-1])/(H*H);
		return (fabs(DDzp)<fabs(DDzc))?DDzp:DDzc;
	}
	else if (k > N - 3)
	{
		DDzn = (2.0*temp[k] - 5.0*temp[k-1] + 4.0*temp[k-2] - temp[k-3])/(H*H);
		DDzc = (temp[k+1] - 2.0*temp[k] + temp[k-1])/(H*H);
		return (fabs(DDzn)<fabs(DDzc))?DDzn:DDzc;
	}
	else
	{
		Dzpp = temp[k+2] - temp[k+1];
		Dzp = temp[k+1] - temp[k];
		Dzn = temp[k] - temp[k-1];
		Dznn = temp[k-1] - temp[k-2];

		DDzp = (2.0*temp[k] - 5.0*temp[k+1] + 4.0*temp[k+2] - temp[k+3])/(H*H);
		DDzc = (temp[k+1] - 2.0*temp[k] + temp[k-1])/(H*H);
		DDzn = (2.0*temp[k] - 5.0*temp[k-1] + 4.0*temp[k-2] - temp[k-3])/(H*H);

//		if ((DDzp!=DDzp)||(DDzc!=DDzc)||(DDzn!=DDzn))
//		{
//			printf("ZDER2 at point (%d, %d, %d), (%lf,%lf,%lf).\n", j, i, k, INT_L+double(j)*H, \
//			INT_L+double(i)*H, INT_L+double(k)*H);
//			printf("distance = %lf, DDzp %lf, DDzc %lf, DDzn %lf.\n", temp[k].dist, DDzp, DDzc, DDzn);
//			exit(0);
//		}
		if ( (Dzp * Dznn <= 0.0 ) && (Dzpp * Dzn <= 0.0))
//			return old_minmod(DDzp, DDzn);
			return min_abs(DDzp, DDzn, DDzc);
		else if (Dzp * Dznn <= 0.0)
			return old_minmod(DDzp, DDzc);
		else if (Dzpp * Dzn <= 0.0)
			return old_minmod(DDzn, DDzc);
		else
			return DDzc;
	}
}


double secondpartialxy(double ***level, int i, int j, int k, int N, double H)
{	
	int di, dj, counter, counti, countj;
	static int flag = 0;
	int sign, signi, signj, signij;

	static double Dxp, Dxpp, Dxn, Dxnn, Dxc;
	static double Dyp, Dypp, Dyn, Dynn, Dyc;
	static double Dypthis, Dyppthis, Dynthis, Dynnthis;
	static double DDxyp, DDxyc, DDxyn, DDxycc;
	static double DDxyelse;
	static double term1, term2, term3, result;
	static double final;

	static double dxy[3][3];
	static const double coef[3] = {1.0, -4.0, 3.0};
	static const double allcoef[3][3] = {{1.0, -4.0, 3.0}, {-4.0, 16.0, -12.0}, {3.0, -12.0, 9.0}};

	if ( (i == 0)||(i == N)||(j == 0)||(j == N))
		return 0.0;
	if ( (i < 2) || (i > N-2) || (j < 2) || (j > N-2))
	{
		return (level[i+1][j+1][k] - level[i-1][j+1][k] - \
			level[i+1][j-1][k] + level[i-1][j-1][k])/(4.0*H*H);
	}


	final = LARGENUM;

	for ( di = 0; di < 3; di++)
	{
		for ( dj = 0; dj < 3; dj++)
		{
			dxy[di][dj] = 0.0;
			if ( (di == 1) && (dj == 1))
			{
				dxy[di][dj] = ( level[i+1][j+1][k] - level[i-1][j+1][k] - \
						level[i+1][j-1][k] + level[i-1][j-1][k]);
			}
			else if  (di == 1)
			{
				sign = dj - 1;
				for ( counter = 2; counter >= 0; counter--)
				{
					dxy[di][dj] += (-1.0 * (double)(sign) * coef[2-counter] * \
						( level[i+1][j+sign*counter][k] - level[i-1][j+sign*counter][k] ) );

				}
//				if (flag < 1)
//				{
//					printf("di = %d, dj = %d, sign = %d, dxy = %lf.\n", di, dj, sign, dxy[di][dj]);
//
//				}
			}
			else if (dj == 1)
			{
				sign = di - 1;
				for ( counter = 2; counter >= 0; counter--)
				{
					dxy[di][dj] += (-1.0 * (double) (sign) * coef[2-counter] * \
						(level[i+sign*counter][j+1][k] - level[i+sign*counter][j-1][k]) );
				}
			}
			else
			{
				signi = di - 1;
				signj = dj - 1;
				signij = signi * signj;
				for ( counti = 2; counti >= 0; counti--)
				{ 
					for ( countj = 2; countj >= 0; countj--)
					{
						dxy[di][dj] += ( allcoef[2-counti][2-countj] * (double)(signij) * \
								level[i+signi*counti][j+signj*countj][k]);
					}
				}
			}
			if (fabs(dxy[di][dj]) < fabs(final))
			{
				final = dxy[di][dj];
			}
		}
	}

	if (flag < 0)
	{
		printf("The array:\n");
		for ( counti = -2; counti < 3; counti++)
		{
			for ( countj = -2; countj < 3; countj++)
			{
				printf("%lf\t", level[i+counti][j+countj][k]);
			}
			printf("\n");
		}
		printf("The candidate partials:\n");
		for ( counti = 0; counti < 3; counti++)
		{
			for ( countj = 0; countj < 3; countj++)
			{
				printf("%lf\t\t", dxy[counti][countj]/(4.0*H*H));
			}
			printf("\n");
		}
		printf("coef: %lf\t%lf\t%lf\n", coef[0], coef[1], coef[2]);
		printf("allcoef:\n");
		for ( counti =0; counti < 3; counti++)
		{
			for ( countj = 0; countj < 3; countj++)
			{
				printf("%lf\t\t", allcoef[counti][countj]);
			}
			printf("\n");
		}

		flag++;
	}


	return final/(4.0*H*H);
}


double secondpartialxz(double ***level, int i, int j, int k, int N, double H)
{
	static int flag = 0;
	int di, dj, counter, counti, countj;
	int sign, signi, signj, signij;

	static double Dxp, Dxpp, Dxn, Dxnn, Dxc;
	static double Dzp, Dzpp, Dzn, Dznn, Dzc;
	static double DDxzp, DDxzn, DDxzc, DDxzcc;
	static double DDxzelse;
	static double term1, term2, term3, result;
	static double final;

	static double dxy[3][3];
	static const double coef[3] = {1.0, -4.0, 3.0};
	static const double allcoef[3][3] = {{1.0, -4.0, 3.0}, {-4.0, 16.0, -12.0}, {3.0, -12.0, 9.0}};

	if ( (k == 0)||(k == N)||(j == 0)||(j == N))
		return 0.0;
	if ( (k < 2) || (k > N-2) || (j < 2) || (j > N-2))
	{
		return (level[i][j+1][k+1] - level[i][j-1][k+1] - \
			level[i][j+1][k-1] + level[i][j-1][k-1])/(4.0*H*H);
	}


	final = LARGENUM;

	for ( di = 0; di < 3; di++)
	{
		for ( dj = 0; dj < 3; dj++)
		{
			dxy[di][dj] = 0.0;
			if ( (di == 1) && (dj == 1))
			{
				dxy[di][dj] = ( level[i][j+1][k+1] - level[i][j-1][k+1] - \
						level[i][j+1][k-1] + level[i][j-1][k-1]);
			}
			else if  (di == 1)
			{
				sign = dj - 1;
				for ( counter = 2; counter >= 0; counter--)
				{
					dxy[di][dj] += (-1.0 * (double)(sign) * coef[2-counter] * \
						( level[i][j+1][k+sign*counter] - level[i][j-1][k+sign*counter] ) );
				}
			}
			else if (dj == 1)
			{
				sign = di - 1;
				for ( counter = 2; counter >= 0; counter--)
				{
					dxy[di][dj] += (-1.0 * (double) (sign) * coef[2-counter] * \
						(level[i][j+sign*counter][k+1] - level[i][j+sign*counter][k-1]) );
				}
			}
			else
			{
				signi = di - 1;
				signj = dj - 1;
				signij = signi * signj;
				for ( counti = 2; counti >= 0; counti--)
				{ 
					for ( countj = 2; countj >= 0; countj--)
					{
						dxy[di][dj] += ( allcoef[2-counti][2-countj] * (double)(signij) * \
								level[i][j+signi*counti][k+signj*countj]);
					}
				}
			}
			if (fabs(dxy[di][dj]) < fabs(final))
			{
				final = dxy[di][dj];
			}
		}
	}

	if (flag < 0)
	{
		printf("The array:\n");
		for ( counti = -2; counti < 3; counti++)
		{
			for ( countj = -2; countj < 3; countj++)
			{
				printf("%lf\t", level[i][j+counti][k+countj]);
			}
			printf("\n");
		}
		printf("The candidate partials:\n");
		for ( counti = 0; counti < 3; counti++)
		{
			for ( countj = 0; countj < 3; countj++)
			{
				printf("%lf\t\t", dxy[counti][countj]/(4.0*H*H));
			}
			printf("\n");
		}
		flag++;
	}

	return final/(4.0*H*H);

}


double secondpartialyz(double ***level, int i, int j, int k, int N, double H)
{	
	static int flag = 0;
	int di, dj, counter, counti, countj;
	int sign, signi, signj, signij;
	static double Dzp, Dzpp, Dzn, Dznn, Dzc;
	static double Dyp, Dypp, Dyn, Dynn, Dyc;
	static double DDyzp, DDyzn, DDyzc, DDyzcc;
	static double DDyzelse;
	static double term1, term2, term3, result;
	static double final;

	static double dxy[3][3];
	static const double coef[3] = {1.0, -4.0, 3.0};
	static const double allcoef[3][3] = {{1.0, -4.0, 3.0}, {-4.0, 16.0, -12.0}, {3.0, -12.0, 9.0}};

	if ( (i == 0)||(i == N)||(k == 0)||(k == N))
		return 0.0;
	if ( (i < 2) || (i > N-2) || (k < 2) || (k > N-2))
	{
		return (level[i+1][j][k+1] - level[i-1][j][k+1] - \
			level[i+1][j][k-1] + level[i-1][j][k-1])/(4.0*H*H);
	}


	final = LARGENUM;

	for ( di = 0; di < 3; di++)
	{
		for ( dj = 0; dj < 3; dj++)
		{
			dxy[di][dj] = 0.0;
			if ( (di == 1) && (dj == 1))
			{
				dxy[di][dj] = ( level[i+1][j][k+1] - level[i-1][j][k+1] - \
						level[i+1][j][k-1] + level[i-1][j][k-1]);
			}
			else if  (di == 1)
			{
				sign = dj - 1;
				for ( counter = 2; counter >= 0; counter--)
				{
					dxy[di][dj] += (-1.0 * (double)(sign) * coef[2-counter] * \
						( level[i+1][j][k+sign*counter] - level[i-1][j][k+sign*counter] ) );
				}
			}
			else if (dj == 1)
			{
				sign = di - 1;
				for ( counter = 2; counter >= 0; counter--)
				{
					dxy[di][dj] += (-1.0 * (double) (sign) * coef[2-counter] * \
						(level[i+sign*counter][j][k+1] - level[i+sign*counter][j][k-1]) );
				}
			}
			else
			{
				signi = di - 1;
				signj = dj - 1;
				signij = signi * signj;
				for ( counti = 2; counti >= 0; counti--)
				{ 
					for ( countj = 2; countj >= 0; countj--)
					{
						dxy[di][dj] += ( allcoef[2-counti][2-countj] * (double)(signij) * \
								level[i+signi*counti][j][k+signj*countj]);
					}
				}
			}
			if (fabs(dxy[di][dj]) < fabs(final))
			{
				final = dxy[di][dj];
			}
		}
	}

	if (flag < 0)
	{
		printf("The array:\n");
		for ( counti = -2; counti < 3; counti++)
		{
			for ( countj = -2; countj < 3; countj++)
			{
				printf("%lf\t", level[i+counti][j][k+countj]);
			}
			printf("\n");
		}
		printf("The candidate partials:\n");
		for ( counti = 0; counti < 3; counti++)
		{
			for ( countj = 0; countj < 3; countj++)
			{
				printf("%lf\t\t", dxy[counti][countj]/(4.0*H*H));
			}
			printf("\n");
		}
		flag++;
	}

	return final/(4.0*H*H);
}
