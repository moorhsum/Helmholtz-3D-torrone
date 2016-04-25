
#include "HelmholtzKernel.h"

int new_HelmholtzKernel(int size, double delta, double epsilon, double *zxlist, double *zylist, double *zzlist, \
	double *zxstarlist, double *zystarlist, double *zzstarlist, double *dist, double *gradx, double *grady, \
	double *gradz, double *curv, double *Jacobian, \
	double *PolyW, double *SineW, double *PolyResult, double *SineResult, \
	double **PolyKernelH, double **PolyKernelHimg, double **SineKernelH, double **SineKernelHimg)
{
	int counter;
	int i, j, k;
	double zx, zy, zz, zstarx, zstary, zstarz, thisdist;
	double thisgradx, thisgrady, thisgradz, stargradx, stargrady, stargradz;
	double xminusy1, xminusy2, xminusy3;
	double starxi, starxj, staryi, staryj, starzi, starzj;
	double partialx, partialy, partialz, partialximg, partialyimg, partialzimg, partial, partialimg;

	double regfactor, regfactorimg, tau, threshold, regterm1 ,regterm2, regterm3, regterm4;

	tau = 0.1 * H;
	threshold = 0.01*H*H;
	
	regterm1 = 1.0/(4.0*PI*tau*radius);		//	1/8pi*tau(2/r)	//
	regterm2 = (-5.0*tau)/(256.0*pow(radius,3.0)*PI);	//	-(1/pi)*(5/512)*(2/r^3)*tau	//
	regterm3 = (-50.0*tau)/(1536.0*pow(radius,3.0)*PI);//	-(1/pi)*(25/1536)*(2/r^3)*tau	//
	regterm4 = (pow(WAVE,2.0)*tau)/(24.0*PI*radius);	//	Helmholtz 3D, (k^2/48)*(k1+k2)*tau	//

	if (INTERIOR_MODE==1)
		regfactor = regterm1 + regterm2 + regterm3 + regterm4;
	else
		regfactor = -1.0*(regterm1 + regterm2 + regterm3 + regterm4);

	regfactorimg = 0.0;



	counter = 0;
	for ( i = 0; i < SIZE; i++)
	{
		zy = INT_L + (double)(i) * H;
	for ( j = 0; j < SIZE; j++)
	{	
		zx = INT_L + (double)(j) * H;
	for ( k = 0; k < SIZE; k++)
	{
		zz = INT_L + (double)(k) * H;
		thisdist = radius - sqrt( (zx-Cx)*(zx-Cx) + (zy-Cy)*(zy-Cy) + (zz-Cz)*(zz-Cz) );
//		thisdist = 	

		if (COMBO_ONESIDED == 0)		
		{
			if ( ( fabs(thisdist) > epsilon) || (fabs(thisdist) < delta) )
				continue;
		}
		else
		{
			if ( (thisdist > epsilon) | (thisdist < delta) )
				continue;
		}

		thisgradx = (Cx - zx) / sqrt( (zx-Cx)*(zx-Cx) + (zy-Cy)*(zy-Cy) + (zz-Cz)*(zz-Cz) );
		thisgrady = (Cy - zy) / sqrt( (zx-Cx)*(zx-Cx) + (zy-Cy)*(zy-Cy) + (zz-Cz)*(zz-Cz) );
		thisgradz = (Cz - zz) / sqrt( (zx-Cx)*(zx-Cx) + (zy-Cy)*(zy-Cy) + (zz-Cz)*(zz-Cz) );
			
		zstarx = zx - thisdist * thisgradx;
		zstary = zy - thisdist * thisgrady;
		zstarz = zz - thisdist * thisgradz;

	stargradx = (Cx - zstarx) / sqrt( (zstarx-Cx)*(zstarx-Cx) + (zstary-Cy)*(zstary-Cy) + (zstarz-Cz)*(zstarz-Cz) );
	stargrady = (Cy - zstary) / sqrt( (zstarx-Cx)*(zstarx-Cx) + (zstary-Cy)*(zstary-Cy) + (zstarz-Cz)*(zstarz-Cz) );
	stargradz = (Cz - zstarz) / sqrt( (zstarx-Cx)*(zstarx-Cx) + (zstary-Cy)*(zstary-Cy) + (zstarz-Cz)*(zstarz-Cz) );

			
		zxlist[counter] = zx;
		zylist[counter] = zy;
		zzlist[counter] = zz;
		zxstarlist[counter] = zstarx;
		zystarlist[counter] = zstary;
		zzstarlist[counter] = zstarz;
		dist[counter] = thisdist;
		gradx[counter] = stargradx;
		grady[counter] = stargrady;
		gradz[counter] = stargradz;
		curv[counter] = 1./(radius - thisdist);
			
			Jacobian[counter] = 1. + thisdist * curv[counter];
		//	Jacobian[counter] = 1.;
		//	Jacobian[counter] = cal_J(J_M, zx, zy, zz);

		if (POLY_TEST == 1)
		{
			if (COMBO_ONESIDED == 0)
				PolyW[counter] = PolyWeight(fabs(thisdist/epsilon), delta/epsilon)/(2.0*epsilon);
			else
				PolyW[counter] = PolyWeight(fabs(thisdist/epsilon), delta/epsilon)/(epsilon);
			
			PolyResult[counter] = HCUBED*Jacobian[counter] * PolyW[counter];
		}
		if (SINE_TEST == 1)
		{	
			if (COMBO_ONESIDED == 0)
				SineW[counter] = SineWeight(fabs(thisdist/epsilon), delta/epsilon)/(2.0*epsilon);
			else
				SineW[counter] = SineWeight(fabs(thisdist/epsilon), delta/epsilon)/(epsilon);
			SineResult[counter] = HCUBED*Jacobian[counter] * SineW[counter];
		}

		counter++;
	}
	}
	}

	if (counter != size)
	{
		printf("Size different for new Kernel.\n");
		exit(0);
	}

	for ( i = 0; i < size; i++)
	{
		starxi = zxstarlist[i];
		staryi = zystarlist[i];
		starzi = zzstarlist[i];
		for ( j = 0; j < size; j++)
		{
			starxj = zxstarlist[j];
			staryj = zystarlist[j];
			starzj = zzstarlist[j];

			partialx = cal_partial(1, starxi, staryi, starzi, starxj, staryj, starzj);// x coordinate real	//
			partialy = cal_partial(2, starxi, staryi, starzi, starxj, staryj, starzj);// y coordinate real	//
			partialz = cal_partial(3, starxi, staryi, starzi, starxj, staryj, starzj);// z coordinate real	//
			partialximg = cal_partial(4, starxi, staryi, starzi, starxj, staryj, starzj);// x coordinate imag
			partialyimg = cal_partial(5, starxi, staryi, starzi, starxj, staryj, starzj);// y coordinate imag
			partialzimg = cal_partial(6, starxi, staryi, starzi, starxj, staryj, starzj);// z coordinate imag

			partial =  -1. * (partialx * gradx[i] + partialy * grady[i] + partialz * gradz[i]);
			partialimg = -1. * (partialximg * gradx[i] + partialyimg * grady[i] + partialzimg * gradz[i]);

			if ( ( (starxi-starxj)*(starxi-starxj) + (staryi-staryj)*(staryi-staryj) + \
				(starzi-starzj)*(starzi-starzj) ) < threshold )
			{
				partial = regfactor;
				partialimg = regfactorimg;
			}
//			else
//			{
//				partial = partial;
//				partialimg = partialimg;
//			}
			
			if (POLY_TEST == 1)
			{
				PolyKernelH[i][j] = partial * PolyResult[j];
				PolyKernelHimg[i][j] = partialimg * PolyResult[j];
				if (i == j)
					PolyKernelH[i][j] -= 0.5;			
			}
			if (SINE_TEST == 1)
			{
				SineKernelH[i][j] = partial * SineResult[j];
				SineKernelHimg[i][j] = partialimg * SineResult[j];
				if (i == j)
					SineKernelH[i][j] -= 0.5;
			}
		}
	}
	return 0;
}



int new_HelmholtzKernelCombo(int size, double delta, double epsilon, double ***level, double ***leveldx, double ***leveldy,\
	double ***leveldz, double ***levelmcurv, double ***levelgcurv, double *zxlist, double *zylist, double *zzlist, \
	double *zxstarlist, double *zystarlist, double *zzstarlist, double *dist, double *gradx, double *grady, \
	double *gradz, double *Jacobian, double *regfactor, double *regfactorimg, \
	double *PolyW, double *SineW, double *PolyResult, double *SineResult, \
	double **PolyKernelH, double **PolyKernelHimg, double **SineKernelH, double **SineKernelHimg, double wavenum, \
	double eta)
{
	static int flag = 0;
	int counter, regcount, badgradnormcount;
	int i, j, k;
	double zx, zy, zz, zstarx, zstary, zstarz, thisdist;
	double thisgradx, thisgrady, thisgradz, stargradx, stargrady, stargradz, stargradnorm;

	double dxx, dyy, dzz, dxy, dxz, dyz, cur_meancurv, cur_gcurv;
	

	static const double fourpi = 4.0 * PI;
	static const double recipfourpi = 1./(4.*PI);
//	double eta;

	double tau, threshold, regterm1, regterm2, regterm3, regterm4;

//	tau = 0.01 * h;
//	threshold = 0.0001*h*h;

	tau = TAU_FACTOR * H;
	threshold = tau*tau;

//	regterm1 = 1.0/(4.0*PI*tau*radius);		//	1/8pi*tau(2/r)	//
//	regterm2 = (-5.0*tau)/(256.0*pow(radius,3.0)*PI);	//	-(1/pi)*(5/512)*(2/r^3)*tau	//
//	regterm3 = (-50.0*tau)/(1536.0*pow(radius,3.0)*PI);//	-(1/pi)*(25/1536)*(2/r^3)*tau	//
//	regterm4 = (pow(WAVE,2.0)*tau)/(24.0*PI*radius);	//	Helmholtz 3D, (k^2/48)*(k1+k2)*tau	//

//	if (INTERIOR_MODE==1)
//		regfactor = regterm1 + regterm2 + regterm3 + regterm4;
//	else
//		regfactor = -1.0*(regterm1 + regterm2 + regterm3 + regterm4);

//	regfactor = 0.0;
//	regfactorimg = 0.0;


////////////////////////////////////////////////////////////
//	eta = 0.5;
////////////////////////////////////////////////////////////

	counter = 0;
	regcount = 0;
	badgradnormcount = 0;
	for ( i = 0; i < SIZE; i++)
	{
		zy = INT_L + (double)(i) * H;
	for ( j = 0; j < SIZE; j++)
	{	
		zx = INT_L + (double)(j) * H;
	for ( k = 0; k < SIZE; k++)
	{
		zz = INT_L + (double)(k) * H;

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////	THIS NEEDS TO BE CHANGED WHEN GENERAL LEVEL SET KICKS IN	////////////////


//		thisdist = radius - sqrt( (zx-Cx)*(zx-Cx) + (zy-Cy)*(zy-Cy) + (zz-Cz)*(zz-Cz) );

		thisdist = level[i][j][k];

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	
		if (COMBO_ONESIDED == 0)
		{
			if ( ( fabs(thisdist) > epsilon) || (fabs(thisdist) < delta) )
				continue;
		}
		else
		{
			if (COMBO_INTERIOR == 1)
			{
				if ( ( thisdist > epsilon) || (thisdist < delta) )
					continue;
			}
			else
			{
				if ( ( thisdist < -1.*epsilon) || (thisdist > -1.*delta) )
					continue;
			}
		}

//		thisgradx = (Cx - zx) / sqrt( (zx-Cx)*(zx-Cx) + (zy-Cy)*(zy-Cy) + (zz-Cz)*(zz-Cz) );
//		thisgrady = (Cy - zy) / sqrt( (zx-Cx)*(zx-Cx) + (zy-Cy)*(zy-Cy) + (zz-Cz)*(zz-Cz) );
//		thisgradz = (Cz - zz) / sqrt( (zx-Cx)*(zx-Cx) + (zy-Cy)*(zy-Cy) + (zz-Cz)*(zz-Cz) );

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////	THIS NEEDS TO BE CHANGED WHEN GENERAL LEVEL SET KICKS IN	////////////////

		
//		thisgradx = (Cx - zx) / (radius - thisdist);
//		thisgrady = (Cy - zy) / (radius - thisdist);
//		thisgradz = (Cz - zz) / (radius - thisdist);
			
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

		zstarx = zx - thisdist * leveldx[i][j][k];
		zstary = zy - thisdist * leveldy[i][j][k];
		zstarz = zz - thisdist * leveldz[i][j][k];

//		stargradx = (Cx-zstarx)/sqrt( (zstarx-Cx)*(zstarx-Cx) + (zstary-Cy)*(zstary-Cy) + (zstarz-Cz)*(zstarz-Cz));
//		stargrady = (Cy-zstary)/sqrt( (zstarx-Cx)*(zstarx-Cx) + (zstary-Cy)*(zstary-Cy) + (zstarz-Cz)*(zstarz-Cz));
//		stargradz = (Cz-zstarz)/sqrt( (zstarx-Cx)*(zstarx-Cx) + (zstary-Cy)*(zstary-Cy) + (zstarz-Cz)*(zstarz-Cz));
		

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////	THIS NEEDS TO BE CHANGED WHEN GENERAL LEVEL SET KICKS IN	////////////////

		cal_interpgrad(N, H, zstarx, zstary, zstarz, leveldx, leveldy, leveldz, &stargradx, &stargrady, &stargradz);
		stargradnorm = sqrt(stargradx*stargradx + stargrady*stargrady + stargradz*stargradz);
		stargradx = stargradx/stargradnorm;
		stargrady = stargrady/stargradnorm;
		stargradz = stargradz/stargradnorm;
//		stargradnorm = sqrt(stargradx*stargradx + stargrady*stargrady + stargradz*stargradz);
//		stargradx = (Cx - zstarx)/radius;
//		stargrady = (Cy - zstary)/radius;
//		stargradz = (Cz - zstarz)/radius;

		if (fabs(stargradnorm - 1.0) > 0.1)
		{
			badgradnormcount++;
		}


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

		zxlist[counter] = zx;
		zylist[counter] = zy;
		zzlist[counter] = zz;
		zxstarlist[counter] = zstarx;
		zystarlist[counter] = zstary;
		zzstarlist[counter] = zstarz;
		dist[counter] = thisdist;
		gradx[counter] = stargradx;
		grady[counter] = stargrady;
		gradz[counter] = stargradz;

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////	THIS NEEDS TO BE CHANGED WHEN GENERAL LEVEL SET KICKS IN	////////////////

//		mcurv[counter] = levelmcurv[i][j][k];
//		gcurv[counter] = levelgcurv[i][j][k];
//		curv[counter] = 1./(radius - thisdist);

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			
//		Jacobian[counter] = 1. + thisdist * curv[counter];
//		Jacobian[counter] = 1.;
//		Jacobian[counter] = cal_J(J_M, zx, zy, zz);

//		dxx = secondpartialxx(zx, zy, zz);
//		dyy = secondpartialyy(zx, zy, zz);
//		dzz = secondpartialzz(zx, zy, zz);
//		dxy = secondpartialxy(zx, zy, zz);
//		dxz = secondpartialxz(zx, zy, zz);
//		dyz = secondpartialyz(zx, zy, zz);
//		cur_meancurv = -0.5 * (dxx + dyy + dzz);
//		cur_gcurv = dxx*dyy + dyy*dzz + dxx*dzz - dxy*dxy - dxz*dxz - dyz*dyz;
		Jacobian[counter] = 1 + 2. * thisdist * levelmcurv[i][j][k] + thisdist*thisdist * levelgcurv[i][j][k];

		///	REGULARIZATION		////////////////

//		cur_meancurv = cal_curv(level, zstarx, zstary, zstarz, N, H);
//		cur_gcurv = cal_GaussianCurv(level, zstarx, zstary, zstarz, N, H);
		cal_BothCurvs( zstarx, zstary, zstarz, H, levelmcurv, levelgcurv, &cur_meancurv, &cur_gcurv);

		regterm1 = cur_meancurv/(4.*PI*tau);		//	(k1+k2)/8pitau

		//	(5tau/512pi)) (k1+k2)(k1^2 + k1k2 + k2^2)	//
		regterm2 = -5.*tau*(2.*cur_meancurv)*(4.*cur_meancurv*cur_meancurv-cur_gcurv)/(512.*PI);	
		//	(25tau/1536pi) (k1+k2) * k1k2	//
		regterm3 = -25.*tau*(2.*cur_meancurv)*cur_gcurv/(1536.*PI);
		//	k^2 tau / (48pi) * (k1+k2)
		regterm4 = WAVE*WAVE * tau * 2.* cur_meancurv / (48.*PI);

		if (INTERIOR_MODE == 1)
			regfactor[counter] = regterm1 + regterm2 + regterm3 + regterm4;
		else
			regfactor[counter] = -1. * (regterm1 + regterm2 + regterm3 + regterm4);
		regfactorimg[counter] = 0.0;
	

		if (POLY_TEST == 1)
		{
			if (COMBO_ONESIDED == 0)
				PolyW[counter] = PolyWeight(fabs(thisdist/epsilon), delta/epsilon)/(2.0*epsilon);
			else
				PolyW[counter] = PolyWeight(fabs(thisdist/epsilon), delta/epsilon)/(epsilon);
				
			PolyResult[counter] = HCUBED*Jacobian[counter] * PolyW[counter];
		}
		if (SINE_TEST == 1)
		{
			if (COMBO_ONESIDED == 0)
				SineW[counter] = SineWeight(fabs(thisdist/epsilon), delta/epsilon)/(2.0*epsilon);
			else
				SineW[counter] = SineWeight(fabs(thisdist/epsilon), delta/epsilon)/(epsilon);
			
			SineResult[counter] = HCUBED*Jacobian[counter] * SineW[counter];
		}
		
		
//		if ( (gradx[counter] != gradx[counter]) || (grady[counter] != grady[counter]) || \
//			(gradz[counter] != gradz[counter]) || (fabs(stargradnorm - 1.0) > 0.3) )
//		{
//		printf("point(%lf, %lf, %lf), grad (%lf, %lf, %lf) project to (%lf, %lf, %lf) grad (%lf, %lf, %lf)\n", \
//		zx, zy, zz, leveldx[i][j][k], leveldy[i][j][k], leveldz[i][j][k], zstarx, zstary, zstarz, \
//		gradx[counter], grady[counter], gradz[counter]);
//		printf("dist %lf, mcurv %lf, gcurv %lf\n", level[i][j][k], levelmcurv[i][j][k], levelgcurv[i][j][k]);
//		flag++;
//		if (flag > 100)
//			exit(0);
//		}
		counter++;
	}
	}
	}

	printf("Total %d points have grad norm error bigger than 0.1\n", badgradnormcount);
	if (counter != size)
	{
		printf("Counter = %d, specified size = %d. Size different for new Kernel.\n", counter, size);
		printf("delta = %lf, epsilon = %lf.\n", delta, epsilon);
		exit(0);
	}

	regcount = 0;
	#pragma omp parallel for private(i,j) schedule(static, 16) reduction(+:regcount)
	for ( i = 0; i < size; i++)
	{
		double xminusy1, xminusy2, xminusy3, xminusynorm;
		double starxi, starxj, staryi, staryj, starzi, starzj;
		double partialx, partialy, partialz, partialximg, partialyimg, partialzimg, partial, partialimg;
		double Polysinglereal, Polysingleimg, Polydoublereal, Polydoubleimg, Polytempreal, Polytempimg;
		double Sinesinglereal, Sinesingleimg, Sinedoublereal, Sinedoubleimg, Sinetempreal, Sinetempimg;

		double x1, x2, x3, parameter, CosineValue, SineValue, commonterm1, commonterm2;
		double inner_pd1, inner_pd2, inner_pd3;
		double term1, term2, term3, term1R, term2R, term3R, term1I, term2I, term3I, term4R, term4I;

		double xmy1, xmy2, xmy3, xmynorm, CosValue, extSineValue, extCosValue, denom;
		double single_inner_pd, recipxminusynorm, recipxmynorm;	//	speed up techniques	//


		starxi = zxstarlist[i];
		staryi = zystarlist[i];
		starzi = zzstarlist[i];
		for ( j = 0; j < size; j++)
		{
			
			starxj = zxstarlist[j];
			staryj = zystarlist[j];
			starzj = zzstarlist[j];

			//	SINGLE LAYER PART, dG(x*,y*)/dnx W(y) J(y) * h^2	//

			xmy1 = starxi - starxj;
			xmy2 = staryi - staryj;
			xmy3 = starzi - starzj;
			xmynorm = sqrt( xmy1*xmy1 + xmy2*xmy2 + xmy3*xmy3 );


			if (xmynorm > tau)
			{
				recipxmynorm = 1./xmynorm;		//	1/|x-y|

				parameter = WAVE * xmynorm;		//	k|x-y|
				SineValue = sin(parameter);		//	sin(k|x-y|)
				CosValue = cos(parameter);		//	cos(k|x-y|)
				extSineValue = SineValue * recipxmynorm;	//	sin(k|x-y|)/|x-y|
				extCosValue = CosValue * recipxmynorm;		//	cos(k|x-y|)/|x-y|
				denom = recipfourpi*recipxmynorm*recipxmynorm;	//	1./4pi|x-y|^2
		
				commonterm1 = -1.* (WAVE * SineValue + extCosValue) * denom;
				commonterm2 = (WAVE * CosValue - extSineValue) * denom;
		

//				partialx = commonterm1 * (xmy1);
//				partialy = commonterm1 * (xmy2);
//				partialz = commonterm1 * (xmy3);
//				partialximg = commonterm2 * (xmy1);
//				partialyimg = commonterm2 * (xmy2);
//				partialzimg = commonterm2 * (xmy3);

				single_inner_pd = -1. * (xmy1 * gradx[i] + xmy2 * grady[i] + xmy3 * gradz[i]);

//			if (0)
//			{
//			partialx = cal_partial(1, starxi, staryi, starzi, starxj, staryj, starzj);// x coordinate real	//
//			partialy = cal_partial(2, starxi, staryi, starzi, starxj, staryj, starzj);// y coordinate real	//
//			partialz = cal_partial(3, starxi, staryi, starzi, starxj, staryj, starzj);// z coordinate real	//
//			partialximg = cal_partial(4, starxi, staryi, starzi, starxj, staryj, starzj);// x coordinate imag
//			partialyimg = cal_partial(5, starxi, staryi, starzi, starxj, staryj, starzj);// y coordinate imag
//			partialzimg = cal_partial(6, starxi, staryi, starzi, starxj, staryj, starzj);// z coordinate imag
			
//			}
			//	ik/4 * H1(k|x-y|)/|x-y| [(x-y) dot nx] = ik/4 * (J1 + iY1)/|x-y| [inner_pd]	//
			
				partial =  commonterm1 * single_inner_pd;
				partialimg = commonterm2 * single_inner_pd;
			}
			else
//			if ( ( (starxi-starxj)*(starxi-starxj) + (staryi-staryj)*(staryi-staryj) + \
//				(starzi-starzj)*(starzi-starzj)) < threshold )
			{
				regcount += 1;
				partial = regfactor[i];
				partialimg = 0.0;
			}
//			else
//			{
//				
//			}
//			else
//			{
//				partial = partial;
//				partialimg = partialimg;
//			}

			if (POLY_TEST == 1)
			{
				Polytempreal = partial * PolyResult[j];
				Polytempimg = partialimg * PolyResult[j];
				if (i == j)
					Polytempreal -= 0.5;
				//	-i * eta * dG(x,y*)/dnx		//
				Polysinglereal = eta * Polytempimg;
				Polysingleimg = -1. * eta * Polytempreal;
			}
			if (SINE_TEST == 1)
			{
				Sinetempreal = partial * SineResult[j];
				Sinetempimg = partialimg * SineResult[j];
				if (i == j)
					Sinetempreal -= 0.5;
				Sinesinglereal = eta * Sinetempimg;
				Sinesingleimg = -1. * eta * Sinetempreal;
			}

			//	DOUBLE LAYER PART, d^2G(x(d(y)), y*)/dnxdny W(y) J(y) h^2	//

			//	x = x(d(y)) = x* + |d(y)| nx	, only outside points	//
			x1 = starxi + fabs(dist[j]) * gradx[i];
			x2 = staryi + fabs(dist[j]) * grady[i];
			x3 = starzi + fabs(dist[j]) * gradz[i];
//			x1 = starxi - fabs(dist[j]) * gradx[i];
//			x2 = staryi - fabs(dist[j]) * grady[i];
//			x3 = starzi - fabs(dist[j]) * gradz[i];

			xminusy1 = x1 - starxj;
			xminusy2 = x2 - staryj;
			xminusy3 = x3 - starzj;
			xminusynorm = sqrt( xminusy1 * xminusy1 + xminusy2 * xminusy2 + xminusy3 * xminusy3);
			recipxminusynorm = 1./xminusynorm;	//	1/|x-y|
			parameter = wavenum * xminusynorm;


/////////////////////////////////	this part is/might be 3D specific	///////////////////////////
//
			inner_pd1 = gradx[i] * gradx[j] + grady[i] * grady[j] + gradz[i] * gradz[j];
			inner_pd2 = xminusy1 * gradx[i] + xminusy2 * grady[i] + xminusy3 * gradz[i];
			inner_pd3 = xminusy1 * gradx[j] + xminusy2 * grady[j] + xminusy3 * gradz[j];

			commonterm1 = inner_pd1 * recipxminusynorm;
			commonterm2 = inner_pd2 * recipxminusynorm * inner_pd3 * recipxminusynorm ;

			CosineValue = cos(parameter) * recipxminusynorm;
			SineValue = sin(parameter) * recipxminusynorm;

			//	term1 = e^ik|x-y| /|x-y|^2 * (ik - 1/|x-y|) * (nx dot ny)		//
			//	-1(cos(k|x-y|)/|x-y|^2 + k * sin(k|x-y|)/|x-y|) (nx dot ny)/|x-y|	//
			term1R = -1. * ( CosineValue * recipxminusynorm + wavenum * SineValue ) * commonterm1;
			//	(k * cos(k|x-y|)/|x-y| - sin(k|x-y|)/|x-y|^2) (nx dot ny)/|x-y|
			term1I = ( wavenum * CosineValue - SineValue * recipxminusynorm ) * commonterm1;
			//	term2 = -k^2 e^ik|x-y| /|x-y|^3	* (nx dot (x-y)) (ny dot (x-y))		//
			//	-k^2 cos(k|x-y|)/|x-y| * (inner_pd2 * inner_pd3)/ |x-y|^2		//
			term2R = -1. * wavenum * wavenum * CosineValue * commonterm2;
			//	-k^2 sin(k|x-y|)/|x-y| * (inner_pd2 * inner_pd3)/ |x-y|^2		//
			term2I = -1. * wavenum * wavenum * SineValue * commonterm2;
			//	term3 = -3ik * e^ik|x-y|/|x-y|^4 * (nx dot (x-y)) (ny dot (x-y))	//
			//	3ksin(k|x-y|)/|x-y| * 1/|x-y| * (inner_pd2 * inner_pd3)/|x-y|^2		//
			term3R = 3. * wavenum * SineValue * commonterm2 * recipxminusynorm;			
			//	-3k cos(k|x-y|)/|x-y| * 1/|x-y| * (inner_pd2 * inner_pd3)/|x-y|^2	//
			term3I = -3. * wavenum * CosineValue * commonterm2 * recipxminusynorm;
			//	term4 = 3e^ik|x-y|/|x-y|^5 * (inner_pd2 * inner_pd3)			//
			//	3cos(k|x-y|)/|x-y| * 1/|x-y|^2 * (inner_pd2 * inner_pd3)/|x-y|^2	//
			term4R = 3. * CosineValue * recipxminusynorm * commonterm2 * recipxminusynorm;
			//	3sin(k|x-y|)/|x-y| * 1/|x-y|^2 * (inner_pd2 * inner)pd3)/|x-y|^2	//
			term4I = 3. * SineValue * recipxminusynorm * commonterm2 * recipxminusynorm;

/////////////////////////////////////////////////////////////////////////////////////////////////////////////

			if (POLY_TEST == 1)
			{
				Polydoublereal = PolyResult[j] * (term1R + term2R + term3R + term4R) * recipfourpi;
				Polydoubleimg = PolyResult[j] * (term1I + term2I + term3I + term4I) * recipfourpi;
				PolyKernelH[i][j] = Polysinglereal + Polydoublereal;
				PolyKernelHimg[i][j] = Polysingleimg + Polydoubleimg;
			}
			if (SINE_TEST == 1)
			{
				Sinedoublereal = SineResult[j] * (term1R + term2R + term3R + term4R) * recipfourpi;
				Sinedoubleimg = SineResult[j] * (term1I + term2I + term3I + term4I) * recipfourpi;
				SineKernelH[i][j] = Sinesinglereal + Sinedoublereal;
				SineKernelHimg[i][j] = Sinesingleimg + Sinedoubleimg;

//				if ( (SineKernelH[i][j] != SineKernelH[i][j]) || (SineKernelHimg[i][j] != SineKernelHimg[i][j]))
//				{
//					printf("i, j = (%d, %d), Sinedouble (%lf, %lf), result %lf, partial (%lf, %lf)\n", \
//					i, j, Sinedoublereal, Sinedoubleimg, SineResult, partial, partialimg);
//					printf("commonterm1 %lf, 2 %lf, |x-y| %lf, Cosvalue %lf, Sinevalue %lf\n", \
//					commonterm1, commonterm2, xminusynorm, CosineValue, SineValue);
//					printf("gradi (%lf, %lf, %lf), gradj (%lf, %lf, %lf) dist %lf\n", \
//					gradx[i], grady[i], gradz[i], gradx[j], grady[j], gradz[j]);
//				}
			}

//			PolyKernelH[i][j] = Polysinglereal;
//			PolyKernelHimg[i][j] = Polysingleimg;
//			SineKernelH[i][j] = Sinesinglereal;
//			SineKernelHimg[i][j] = Sinesingleimg;

		}
	}
	printf("Regcount = %d, on average %d regs per point.\n", regcount, regcount/size);
	return regcount;
}


int DoubleLayerKernel(int size, double epsilon, double ***level, double ***leveldx, double ***leveldy, double ***leveldz, \
	double ***levelmcurv, double ***levelgcurv, double *zxlist, double *zylist, double *zzlist, \
	double *zxstarlist, double *zystarlist, double *zzstarlist, double *dist, double *gradx, double *grady, \
	double *gradz, double *delta, double *Jacobian, double *regfactor, \
	double *regfactorimg, double *result, double **Kernel, double **Kernelimg, double wavenum)
{
	int i, j, k, counti;
	double tau, threshold, term1, term2, term3, term4;
	double thisdist, thisnorm2;
	double tempzx, tempzy, tempzz, tstarx, tstary, tstarz;
	double tempdelta, tempJ;
	double starxi, staryi, starzi, starxj, staryj, starzj;
	double partialx, partialy, partialz, partialximg, partialyimg, partialzimg, partial, partialimg;

	double cur_meancurv, cur_gcurv;

	tau = H * TAU_FACTOR;
	threshold = tau * tau;

//	term1 = 1./(4.*PI*tau*radius);
//	term2 = (-5.*tau)/(256.*pow(radius,3.)*PI);
//	term3 = -50.*tau/(1536.*pow(radius,3.)*PI);
//	term4 = pow(WAVE,2.)*tau/(24.*PI*radius);

//	regfactor = -1.*(term1+term2+term3+term4);
//	regfactorimg = 0.;
	counti = 0;

	for (i = 0; i < SIZE; i++)
	{
		tempzy = INT_L + H*(double)(i);
	for (j = 0; j < SIZE ;j++)
	{
		tempzx = INT_L + H*(double)(j);
	for (k = 0; k < SIZE; k++)
	{
		tempzz = INT_L + H*(double)(k);

//		thisdist = distance(tempzx, tempzy, tempzz);
		thisdist = level[i][j][k];
		if (TRAD_ONESIDED == 0)
		{
			if (fabs(thisdist) > epsilon)
				continue;
		}
		else
		{
			if ((thisdist > epsilon) || (thisdist < 0.) )
				continue;
		}

		zxlist[counti] = tempzx;
		zylist[counti] = tempzy;
		zzlist[counti] = tempzz;
//		tstarx = starx(tempzx, tempzy, tempzz);
//		tstary = stary(tempzx, tempzy, tempzz);
//		tstarz = starz(tempzx, tempzy, tempzz);

		tstarx = tempzx - thisdist * leveldx[i][j][k];
		tstary = tempzy - thisdist * leveldy[i][j][k];
		tstarz = tempzz - thisdist * leveldz[i][j][k];

		zxstarlist[counti] = tstarx;
		zystarlist[counti] = tstary;
		zzstarlist[counti] = tstarz;
		dist[counti] = thisdist;

		if (TRAD_ONESIDED == 0)
			tempdelta = cal_delta(epsilon, thisdist);
		else
			tempdelta = 2. * cal_delta(epsilon, thisdist);

		delta[counti] = tempdelta;
//		tempJ = cal_J(J_M, tempzx, tempzy, tempzz);
		tempJ = 1. + 2. * thisdist * levelmcurv[i][j][k] + thisdist * thisdist * levelgcurv[i][j][k];
		Jacobian[counti] = tempJ;
//		gradx[counti] = cal_interpgrad(1, tstarx, tstary, tstarz);
//		grady[counti] = cal_interpgrad(2, tstarx, tstary, tstarz);
//		gradz[counti] = cal_interpgrad(3, tstarx, tstary, tstarz);
		cal_interpgrad(N,H,tstarx, tstary, tstarz, leveldx, leveldy, leveldz, \
			&gradx[counti], &grady[counti], &gradz[counti]);

		///	REGULARIZATION		////////////////

//		cur_meancurv = cal_curv(level, tstarx, tstary, tstarz, N, H);
//		cur_gcurv = cal_GaussianCurv(level, tstarx, tstary, tstarz, N, H);
		cal_BothCurvs( tstarx, tstary, tstarz, H, levelmcurv, levelgcurv, &cur_meancurv, &cur_gcurv);

		term1 = cur_meancurv/(4.*PI*tau);		//	(k1+k2)/8pitau

		//	(5tau/512pi)) (k1+k2)(k1^2 + k1k2 + k2^2)	//
		term2 = -5.*tau*(2.*cur_meancurv)*(4.*cur_meancurv*cur_meancurv-cur_gcurv)/(512.*PI);	
		//	(25tau/1536pi) (k1+k2) * k1k2	//
		term3 = -25.*tau*(2.*cur_meancurv)*cur_gcurv/(1536.*PI);
		//	k^2 tau / (48pi) * (k1+k2)
		term4 = WAVE*WAVE * tau * 2.* cur_meancurv / (48.*PI);

		regfactor[counti] = term1 + term2 + term3 + term4;
//		regfactor[counter] = -1. * (regterm1 + regterm2 + regterm3 + regterm4);
		regfactorimg[counti] = 0.0;

		result[counti] = HCUBED * tempdelta * tempJ;
		counti++;
	}	//	k
	}	//	j
	}	//	i

	if (size != counti)
	{
		printf("size= %d, counti = %d.\n", size, counti);
		exit(0);
	}

	for (i = 0; i < counti; i++)
	{
	for (j = 0; j < counti; j++)
	{
		starxi = zxstarlist[i];
		staryi = zystarlist[i];
		starzi = zzstarlist[i];
		starxj = zxstarlist[j];
		staryj = zystarlist[j];
		starzj = zzstarlist[j];

		thisnorm2 = (starxi-starxj)*(starxi-starxj)+(staryi-staryj)*(staryi-staryj)+(starzi-starzj)*(starzi-starzj);

		partialx = cal_partial(1, starxi, staryi, starzi, starxj, staryj, starzj);
		partialy = cal_partial(2, starxi, staryi, starzi, starxj, staryj, starzj);
		partialz = cal_partial(3, starxi, staryi, starzi, starxj, staryj, starzj);
		partialximg = cal_partial(4, starxi, staryi, starzi, starxj, staryj, starzj);
		partialyimg = cal_partial(5, starxi, staryi, starzi, starxj, staryj, starzj);
		partialzimg = cal_partial(6, starxi, staryi, starzi, starxj, staryj, starzj);

		partial = partialx * gradx[i] + partialy * grady[i] + partialz * gradz[i];
		partialimg = partialximg * gradx[i] + partialyimg * grady[i] + partialzimg * gradz[i];

		if (thisnorm2 < threshold)
		{
			partial = regfactor[i];
			partialimg = 0.0;
		}
		Kernel[i][j] = result[j] * partial;
		Kernelimg[i][j] = result[j] * partialimg;
		if (i == j)
			Kernel[i][j] += 0.5;
	}	//	j
	}	//	i
	return 0;
}
