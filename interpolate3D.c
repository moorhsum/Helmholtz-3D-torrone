
#include "interpolate3D.h"


int cal_interpgrad(int N, double H, double xj, double yj, double zj, double ***gradx, double ***grady, double ***gradz,\
		double *interpolatex, double *interpolatey, double *interpolatez)
{
//	static double interpolant;
//	static POINT tmplll, tmpllh, tmplhl, tmplhh, tmphll, tmphlh, tmphhl, tmphhh;

	double x, y, z, lowx, lowy, lowz;
//	double h = (INT_U-INT_L)/N;
	double temphh, temphl, templh, templl, temp, temp2;
	double templll, templlh, templhl, templhh, temphll, temphlh, temphhl, temphhh;
	double interpolatenorm;	

	int i, j, k;

	static int flag = 0;
	
	i = (int)((yj-INT_L)/H);
	j = (int)((xj-INT_L)/H);
	k = (int)((zj-INT_L)/H);

//	x = INT_L + H*double(i);
//	y = INT_L + H*double(j);
	lowx = xj - INT_L - (double)(j) * H;
	lowy = yj - INT_L - (double)(i) * H;
	lowz = zj - INT_L - (double)(k) * H;

	if ( (lowx > H) ||(lowx < 0.0) || (lowy > H) || (lowy < 0.0) || (lowz < 0.0 ) || (lowz > H) )
	{
		printf("At (%lf, %lf, %lf), low (%lf, %lf, %lf), index = (%d, %d, %d), H = %lf\n", \
		xj, yj, zj, lowx, lowy, lowz, j, i, k, H);
		exit(0);
	}


	//////////////////////////	ONLY INTERPOLATE DIRECTLY	////////////////////////////////////////

	templll = gradx[i][j][k];
	templlh = gradx[i][j][k+1];
	templhl = gradx[i+1][j][k];
	templhh = gradx[i+1][j][k+1];
	temphll = gradx[i][j+1][k];
	temphlh = gradx[i][j+1][k+1];
	temphhl = gradx[i+1][j+1][k];
	temphhh = gradx[i+1][j+1][k+1];
	*interpolatex = ( templll * (H-lowx)*(H-lowy)*(H-lowz) + templlh * (H-lowx) * (H-lowy) * lowz + \
			templhl * (H-lowx) * lowy * (H-lowz) + templhh * (H-lowx) * lowy * lowz + \
			temphll * lowx*(H-lowy)*(H-lowz) + temphlh * lowx * (H-lowy) * lowz + \
			temphhl * lowx * lowy * (H-lowz) + temphhh * lowx * lowy * lowz )/(H*H*H);

	templll = grady[i][j][k];
	templlh = grady[i][j][k+1];
	templhl = grady[i+1][j][k];
	templhh = grady[i+1][j][k+1];
	temphll = grady[i][j+1][k];
	temphlh = grady[i][j+1][k+1];
	temphhl = grady[i+1][j+1][k];
	temphhh = grady[i+1][j+1][k+1];
	*interpolatey = ( templll * (H-lowx)*(H-lowy)*(H-lowz) + templlh * (H-lowx) * (H-lowy) * lowz + \
			templhl * (H-lowx) * lowy * (H-lowz) + templhh * (H-lowx) * lowy * lowz + \
			temphll * lowx*(H-lowy)*(H-lowz) + temphlh * lowx * (H-lowy) * lowz + \
			temphhl * lowx * lowy * (H-lowz) + temphhh * lowx * lowy * lowz )/(H*H*H);

	templll = gradz[i][j][k];
	templlh = gradz[i][j][k+1];
	templhl = gradz[i+1][j][k];
	templhh = gradz[i+1][j][k+1];
	temphll = gradz[i][j+1][k];
	temphlh = gradz[i][j+1][k+1];
	temphhl = gradz[i+1][j+1][k];
	temphhh = gradz[i+1][j+1][k+1];

	*interpolatez = ( templll * (H-lowx)*(H-lowy)*(H-lowz) + templlh * (H-lowx) * (H-lowy) * lowz + \
			templhl * (H-lowx) * lowy * (H-lowz) + templhh * (H-lowx) * lowy * lowz + \
			temphll * lowx*(H-lowy)*(H-lowz) + temphlh * lowx * (H-lowy) * lowz + \
			temphhl * lowx * lowy * (H-lowz) + temphhh * lowx * lowy * lowz )/(H*H*H);

	if ( (*interpolatex != *interpolatex) || (*interpolatey != *interpolatey) || (*interpolatez != *interpolatez) )
	{
		printf("NaN in interplate gradient.\n");
		exit(0);
	}

	interpolatenorm = sqrt( (*interpolatex) * (*interpolatex) + (*interpolatey) * (*interpolatey) + \
				(*interpolatez) * (*interpolatez) );
	if ( fabs(interpolatenorm - 1.0) > 0.5)
	{
		printf("At (%lf, %lf, %lf), gradnorm = %lf, grad (%lf, %lf, %lf)\n", xj, yj, zj, interpolatenorm, \
		*interpolatex, *interpolatey, *interpolatez);
		*interpolatex = *interpolatex/interpolatenorm;
		*interpolatey = *interpolatey/interpolatenorm;
		*interpolatez = *interpolatez/interpolatenorm;
		return 1;
	}

	*interpolatex = *interpolatex/interpolatenorm;
	*interpolatey = *interpolatey/interpolatenorm;
	*interpolatez = *interpolatez/interpolatenorm;
	return 0;

	//////////////////////////////////////////////	INTERPOLATE GRADX	/////////////////////////////////////

	templll = gradx[i][j][k];
	templlh = gradx[i][j][k+1];
	templhl = gradx[i+1][j][k];
	templhh = gradx[i+1][j][k+1];
	temphll = gradx[i][j+1][k];
	temphlh = gradx[i][j+1][k+1];
	temphhl = gradx[i+1][j+1][k];
	temphhh = gradx[i+1][j+1][k+1];


	if ( ((templll>0)&&(templlh>0)&&(templhl>0)&&(templhh>0)&&(temphll>0)&&(temphlh>0)&&(temphhl>0)&&(temphhh>0) ) \
	 || ( (templll<0)&&(templlh<0)&&(templhl<0)&&(templhh<0)&&(temphll<0)&&(temphlh<0)&&(temphhl<0)&&(temphhh<0) ) )
	{
		*interpolatex = ( templll * (H-lowx)*(H-lowy)*(H-lowz) + templlh * (H-lowx) * (H-lowy) * lowz + \
				templhl * (H-lowx) * lowy * (H-lowz) + templhh * (H-lowx) * lowy * lowz + \
				temphll * lowx*(H-lowy)*(H-lowz) + temphlh * lowx * (H-lowy) * lowz + \
				temphhl * lowx * lowy * (H-lowz) + temphhh * lowx * lowy * lowz )/(H*H*H);
	}
	else if (lowx < H-lowx)
	{
		if (lowy < H-lowy)
		{
			if (lowz < H-lowz)
				*interpolatex = templll;
			else
				*interpolatex = templlh;
		}
		else
		{
			if (lowz < H-lowz)
				*interpolatex = templhl;
			else
				*interpolatex = templhh;
		}
	}
	else
	{
		if (lowy < H-lowy)
		{
			if (lowz < H-lowz)
				*interpolatex = temphll;
			else
				*interpolatex = temphlh;
		}
		else
		{
			if (lowz < H-lowz)
				*interpolatex = temphhl;
			else
				*interpolatex = temphhh;
		}
	}

	//////////////////////////////////////////////	INTERPOLATE GRADY	/////////////////////////////////////////

	templll = grady[i][j][k];
	templlh = grady[i][j][k+1];
	templhl = grady[i+1][j][k];
	templhh = grady[i+1][j][k+1];
	temphll = grady[i][j+1][k];
	temphlh = grady[i][j+1][k+1];
	temphhl = grady[i+1][j+1][k];
	temphhh = grady[i+1][j+1][k+1];


	if ( ((templll>0)&&(templlh>0)&&(templhl>0)&&(templhh>0)&&(temphll>0)&&(temphlh>0)&&(temphhl>0)&&(temphhh>0) ) \
	 || ( (templll<0)&&(templlh<0)&&(templhl<0)&&(templhh<0)&&(temphll<0)&&(temphlh<0)&&(temphhl<0)&&(temphhh<0) ) )
	{
		*interpolatey = ( templll * (H-lowx)*(H-lowy)*(H-lowz) + templlh * (H-lowx) * (H-lowy) * lowz + \
				templhl * (H-lowx) * lowy * (H-lowz) + templhh * (H-lowx) * lowy * lowz + \
				temphll * lowx*(H-lowy)*(H-lowz) + temphlh * lowx * (H-lowy) * lowz + \
				temphhl * lowx * lowy * (H-lowz) + temphhh * lowx * lowy * lowz )/(H*H*H);
	}
	else if (lowx < H-lowx)
	{
		if (lowy < H-lowy)
		{
			if (lowz < H-lowz)
				*interpolatey = templll;
			else
				*interpolatey = templlh;
		}
		else
		{
			if (lowz < H-lowz)
				*interpolatey = templhl;
			else
				*interpolatey = templhh;
		}
	}
	else
	{
		if (lowy < H-lowy)
		{
			if (lowz < H-lowz)
				*interpolatey = temphll;
			else
				*interpolatey = temphlh;
		}
		else
		{
			if (lowz < H-lowz)
				*interpolatey = temphhl;
			else
				*interpolatey = temphhh;
		}
	}


	//////////////////////////////////	INTERPOLATE GRADZ	////////////////////////////////////////////////


	templll = gradz[i][j][k];
	templlh = gradz[i][j][k+1];
	templhl = gradz[i+1][j][k];
	templhh = gradz[i+1][j][k+1];
	temphll = gradz[i][j+1][k];
	temphlh = gradz[i][j+1][k+1];
	temphhl = gradz[i+1][j+1][k];
	temphhh = gradz[i+1][j+1][k+1];


	if ( ((templll>0)&&(templlh>0)&&(templhl>0)&&(templhh>0)&&(temphll>0)&&(temphlh>0)&&(temphhl>0)&&(temphhh>0) ) \
	 || ( (templll<0)&&(templlh<0)&&(templhl<0)&&(templhh<0)&&(temphll<0)&&(temphlh<0)&&(temphhl<0)&&(temphhh<0) ) )
	{
		*interpolatez = ( templll * (H-lowx)*(H-lowy)*(H-lowz) + templlh * (H-lowx) * (H-lowy) * lowz + \
				templhl * (H-lowx) * lowy * (H-lowz) + templhh * (H-lowx) * lowy * lowz + \
				temphll * lowx*(H-lowy)*(H-lowz) + temphlh * lowx * (H-lowy) * lowz + \
				temphhl * lowx * lowy * (H-lowz) + temphhh * lowx * lowy * lowz )/(H*H*H);
	}
	else if (lowx < H-lowx)
	{
		if (lowy < H-lowy)
		{
			if (lowz < H-lowz)
				*interpolatez = templll;
			else
				*interpolatez = templlh;
		}
		else
		{
			if (lowz < H-lowz)
				*interpolatez = templhl;
			else
				*interpolatez = templhh;
		}
	}
	else
	{
		if (lowy < H-lowy)
		{
			if (lowz < H-lowz)
				*interpolatez = temphll;
			else
				*interpolatez = temphlh;
		}
		else
		{
			if (lowz < H-lowz)
				*interpolatez = temphhl;
			else
				*interpolatez = temphhh;
		}
	}





	if ((*interpolatex+1==*interpolatex)||(*interpolatex!=*interpolatex) || (fabs(*interpolatex) > 1.7) || \
	(*interpolatey+1 == *interpolatey)||(*interpolatey != *interpolatey) || (fabs(*interpolatey) > 1.7) || \
	(*interpolatez+1 == *interpolatez)||(*interpolatez != *interpolatez) || (fabs(*interpolatez) > 1.7))
	{
		printf("it's gradient (%lf, %lf, %lf), xj = %lf, yj = %lf, zj = %lf.\n", \
		*interpolatex, *interpolatey, *interpolatez, xj, yj, zj);
		printf("i = %d, j = %d, k = %d, point is (%lf, %lf, %lf).\n", i, j, k, INT_L + (double)(j)*H, \
			INT_L+(double)(i)*H, INT_L + (double)(k)*H);
		printf("lowx = %lf, lowy = %lf, lowz = %lf.\n", lowx, lowy, lowz);
//		printf("ZDERlll = %lf, ZDERlhl = %lf, ZDERhll = %lf, ZDERhhl = %lf, ZDERllh = %lf, ZDERlhh = %lf, ZDERhlh = %lf, \
//		ZDERhhh = %lf\n", ZDER_C(i, j, k, 1), ZDER_C(i+1, j, k, 1), ZDER_C(i, j+1, k, 1), ZDER_C(i+1, j+1, k, 1), \
//		ZDER_C(i, j, k+1, 1), ZDER_C(i+1, j, k+1, 1), ZDER_C(i, j+1, k+1, 1), ZDER_C(i+1, j+1, k+1, 1));
		printf("XDERlll = %lf, XDERlhl = %lf, XDERhll = %lf, XDERhhl = %lf, XDERllh = %lf, XDERlhh = %lf, \
		XDERhlh = %lf, XDERhhh = %lf\n", gradx[i][j][k], gradx[i][j+1][k], gradx[i+1][j][k], gradx[i+1][j+1][k], \
		gradx[i][j][k+1], gradx[i+1][j][k+1], gradx[i][j+1][k+1], gradx[i+1][j+1][k+1]);
		printf("YDERlll = %lf, YDERlhl = %lf, YDERhll = %lf, YDERhhl = %lf, YDERllh = %lf, YDERlhh = %lf, \
		YDERhlh = %lf, YDERhhh = %lf\n", grady[i][j][k], grady[i][j+1][k], grady[i+1][j][k], grady[i+1][j+1][k], \
		grady[i][j][k+1], grady[i+1][j][k+1], grady[i][j+1][k+1], grady[i+1][j+1][k+1]);
		printf("ZDERlll = %lf, ZDERlhl = %lf, ZDERhll = %lf, ZDERhhl = %lf, ZDERllh = %lf, ZDERlhh = %lf, \
		ZDERhlh = %lf, ZDERhhh = %lf\n", gradx[i][j][k], gradz[i][j+1][k], gradz[i+1][j][k], gradz[i+1][j+1][k], \
		gradz[i][j][k+1], gradz[i+1][j][k+1], gradz[i][j+1][k+1], gradz[i+1][j+1][k+1]);

		if ( (fabs(*interpolatex) > 1.7) || (fabs(*interpolatey) > 1.7) || (fabs(*interpolatez) > 1.7) )
			return *interpolatex;

		if ( !(fabs(*interpolatex) > 1.7) )
			exit(0);
	}
	return *interpolatex;
}



double cal_curv( double ***dist, double x, double y, double z, int N, double H)
{
	double dxx, dyy, dzz;
	double lowx, lowy, lowz;
	int i, j, k;
	j = (int)((x - INT_L)/H);
	i = (int)((y - INT_L)/H);
	k = (int)((z - INT_L)/H);

	lowx = x - INT_L - (double)(j) * H;
	lowy = y - INT_L - (double)(i) * H;
	lowz = z - INT_L - (double)(k) * H;


	dxx = XDER2(dist, i, j, k, N, H) * (H-lowx) * (H-lowy) * (H-lowz) + \
		XDER2(dist, i+1, j, k, N, H) * (H - lowx) * lowy * (H-lowz)+ \
		XDER2(dist, i, j+1, k, N, H) * lowx * (H-lowy) * (H-lowz) + \
		XDER2(dist, i+1, j+1, k, N, H) * lowx * lowy * (H-lowz) + \
		XDER2(dist, i, j, k+1, N, H) * (H - lowx) * (H-lowy) *lowz + \
		XDER2(dist, i+1, j, k+1, N, H) * (H - lowx) * lowy * lowz + \
		XDER2(dist, i, j+1, k+1, N, H) * lowx * (H-lowy) * lowz + \
		XDER2(dist, i+1, j+1, k+1, N, H) * lowx * lowy * lowz;

	dyy = YDER2(dist, i, j, k, N, H) * (H-lowx) * (H-lowy) * (H-lowz) + \
		YDER2(dist, i+1, j, k, N, H) * (H - lowx) * lowy * (H-lowz)+ \
		YDER2(dist, i, j+1, k, N, H) * lowx * (H-lowy) * (H-lowz) + \
		YDER2(dist, i+1, j+1, k, N, H) * lowx * lowy * (H-lowz) + \
		YDER2(dist, i, j, k+1, N, H) * (H - lowx) * (H-lowy) *lowz + \
		YDER2(dist, i+1, j, k+1, N, H) * (H - lowx) * lowy * lowz + \
		YDER2(dist, i, j+1, k+1, N, H) * lowx * (H-lowy) * lowz + \
		YDER2(dist, i+1, j+1, k+1, N, H) * lowx * lowy * lowz;

	dzz = ZDER2(dist, i, j, k, N, H) * (H-lowx) * (H-lowy) * (H-lowz) + \
		ZDER2(dist, i+1, j, k, N, H) * (H - lowx) * lowy * (H-lowz)+ \
		ZDER2(dist, i, j+1, k, N, H) * lowx * (H-lowy) * (H-lowz) + \
		ZDER2(dist, i+1, j+1, k, N, H) * lowx * lowy * (H-lowz) + \
		ZDER2(dist, i, j, k+1, N, H) * (H - lowx) * (H-lowy) *lowz + \
		ZDER2(dist, i+1, j, k+1, N, H) * (H - lowx) * lowy * lowz + \
		ZDER2(dist, i, j+1, k+1, N, H) * lowx * (H-lowy) * lowz + \
		ZDER2(dist, i+1, j+1, k+1, N, H) * lowx * lowy * lowz;



	if ((dxx+dyy+dzz+1 == dxx+dyy+dzz) || (dxx+dyy+dzz!=dxx+dyy+dzz))
	{
		printf("dxx = %lf, dyy = %lf , dzz = %lf, at (%lf, %lf, %lf).\n", dxx, dyy, dzz, x, y, z);
		printf("Interpolated with points starting at (%lf, %lf, %lf).\n", x-lowx, y-lowy, z-lowz);
		exit(0);
	}
	return -0.5 * (dxx+dyy+dzz)/(H*H*H);
}



double cal_GaussianCurv(double ***dist, double x, double y, double z, int N, double H)
{
	double dxx, dyy, dzz;
	double lowx, lowy, lowz;
	double tmpgausscurv;
	int i, j, k;

	j = (int)((x - INT_L)/H);
	i = (int)((y - INT_L)/H);
	k = (int)((z - INT_L)/H);

	lowx = x - INT_L - (double)(j) * H;
	lowy = y - INT_L - (double)(i) * H;
	lowz = z - INT_L - (double)(k) * H;

	



	tmpgausscurv = ( GaussianCurv(dist, i, j, k, N, H) * (H-lowx) * (H-lowy) * (H-lowz) + \
			 GaussianCurv(dist, i+1, j, k, N, H) * (H-lowx) * lowy * (H-lowz) + \
			 GaussianCurv(dist, i, j, k+1, N, H) * (H-lowx) * (H-lowy) * lowz + \
			 GaussianCurv(dist, i+1, j, k+1, N, H) * (H-lowx) * lowy * lowz + \
			 GaussianCurv(dist, i, j+1, k, N, H) * lowx * (H-lowy) * (H-lowz) + \
			 GaussianCurv(dist, i+1, j+1, k, N, H) * lowx * lowy * (H-lowz) + \
			 GaussianCurv(dist, i, j+1, k+1, N, H) * lowx * (H-lowy) * lowz + \
			 GaussianCurv(dist, i+1, j+1, k+1, N, H) * lowx * lowy * lowz )/(H*H*H);
	if ( (tmpgausscurv!=tmpgausscurv)||(tmpgausscurv+1==tmpgausscurv) )
	{
		printf("Gauss Curvature anomaly at (%lf, %lf, %lf).\nCurvature at each point is LLL:%lf, LHL:%lf, LLH:%lf, LHH:%lf, HLL:%lf, HHL:%lf, HLH:%lf, HHH:%lf.\n", x, y, z, \
			 GaussianCurv(dist, i, j, k, N, H) * (H-lowx) * (H-lowy) * (H-lowz), \
			 GaussianCurv(dist, i+1, j, k, N, H) * (H-lowx) * lowy * (H-lowz), \
			 GaussianCurv(dist, i, j, k+1, N, H) * (H-lowx) * (H-lowy) * lowz, \
			 GaussianCurv(dist, i+1, j, k+1, N, H) * (H-lowx) * lowy * lowz, \
			 GaussianCurv(dist, i, j+1, k, N, H) * lowx * (H-lowy) * (H-lowz), \
			 GaussianCurv(dist, i+1, j+1, k, N, H) * lowx * lowy * (H-lowz),  \
			 GaussianCurv(dist, i, j+1, k+1, N, H) * lowx * (H-lowy) * lowz, \
			 GaussianCurv(dist, i+1, j+1, k+1, N, H) * lowx * lowy * lowz);
		exit(0);
	}
	return tmpgausscurv;
}



double cal_dist(double ***dist, double x, double y, double z, int N, double H)
{
	int indx, indy, indz;
	double lowx, lowy, lowz, smallestxx, smallestyy, smallestzz, tempd2, tempd1, tempd11, tempd12, tempd21, tempd22;
	double tempd111, tempd112, tempd121, tempd122, tempd211, tempd212, tempd221, tempd222;

	indx = (int)((x - INT_L)/H);
	indy = (int)((y - INT_L)/H);
	indz = (int)((z - INT_L)/H);

//	printf("indx = %d, indy = %d\n", indx, indy);

//	smallestxx = LARGENUM;
//	smallestyy = LARGENUM;

	lowx = x - INT_L - indx * H;
	lowy = y - INT_L - indy * H;
	lowz = z - INT_L - indz * H;

	tempd111 = XDER2(dist, indy, indx, indz, N, H);
	tempd211 = XDER2(dist, indy+1, indx, indz, N, H);
	tempd121 = XDER2(dist, indy, indx+1, indz, N, H);
	tempd221 = XDER2(dist, indy+1, indx+1, indz, N, H);
	tempd112 = XDER2(dist, indy, indx, indz+1, N, H);
	tempd212 = XDER2(dist, indy+1, indx, indz+1, N, H);
	tempd122 = XDER2(dist, indy, indx+1, indz+1, N, H);
	tempd222 = XDER2(dist, indy+1, indx+1, indz+1, N, H);


	if (fabs(tempd111) < fabs(tempd112))
		tempd11 = tempd111;
	else
		tempd11 = tempd112;
	if (fabs(tempd121) < fabs(tempd122))
		tempd12 = tempd121;
	else
		tempd12 = tempd122;
	if (fabs(tempd211) < fabs(tempd212))
		tempd21 = tempd211;
	else
		tempd21 = tempd212;
	if (fabs(tempd221) < fabs(tempd222))
		tempd22 = tempd221;
	else
		tempd22 = tempd222;


	if (fabs(tempd11) < fabs(tempd12))
		tempd1 = tempd11;
	else
		tempd1 = tempd12;
	if (fabs(tempd21) < fabs(tempd22) )
		tempd2 = tempd21;
	else
		tempd2 = tempd22;
	if (fabs(tempd1) < fabs(tempd2))
		smallestxx = tempd1;
	else
		smallestxx = tempd2;

	tempd111 = YDER2(dist, indy, indx, indz, N, H);
	tempd211 = YDER2(dist, indy+1, indx, indz, N, H);
	tempd121 = YDER2(dist, indy, indx+1, indz, N, H);
	tempd221 = YDER2(dist, indy+1, indx+1, indz, N, H);
	tempd112 = YDER2(dist, indy, indx, indz+1, N, H);
	tempd212 = YDER2(dist, indy+1, indx, indz+1, N, H);
	tempd122 = YDER2(dist, indy, indx+1, indz+1, N, H);
	tempd222 = YDER2(dist, indy+1, indx+1, indz+1, N, H);


	if (fabs(tempd111) < fabs(tempd112))
		tempd11 = tempd111;
	else
		tempd11 = tempd112;
	if (fabs(tempd121) < fabs(tempd122))
		tempd12 = tempd121;
	else
		tempd12 = tempd122;
	if (fabs(tempd211) < fabs(tempd212))
		tempd21 = tempd211;
	else
		tempd21 = tempd212;
	if (fabs(tempd221) < fabs(tempd222))
		tempd22 = tempd221;
	else
		tempd22 = tempd222;
	if (fabs(tempd11) < fabs(tempd12))
		tempd1 = tempd11;
	else
		tempd1 = tempd12;
	if (fabs(tempd21) < fabs(tempd22) )
		tempd2 = tempd21;
	else
		tempd2 = tempd22;
	if (fabs(tempd1) < fabs(tempd2))
		smallestyy = tempd1;
	else
		smallestyy = tempd2;


	tempd111 = ZDER2(dist, indy, indx, indz, N, H);
	tempd211 = ZDER2(dist, indy+1, indx, indz, N, H);
	tempd121 = ZDER2(dist, indy, indx+1, indz, N, H);
	tempd221 = ZDER2(dist, indy+1, indx+1, indz, N, H);
	tempd112 = ZDER2(dist, indy, indx, indz+1, N, H);
	tempd212 = ZDER2(dist, indy+1, indx, indz+1, N, H);
	tempd122 = ZDER2(dist, indy, indx+1, indz+1, N, H);
	tempd222 = ZDER2(dist, indy+1, indx+1, indz+1, N, H);


	if (fabs(tempd111) < fabs(tempd112))
		tempd11 = tempd111;
	else
		tempd11 = tempd112;
	if (fabs(tempd121) < fabs(tempd122))
		tempd12 = tempd121;
	else
		tempd12 = tempd122;
	if (fabs(tempd211) < fabs(tempd212))
		tempd21 = tempd211;
	else
		tempd21 = tempd212;
	if (fabs(tempd221) < fabs(tempd222))
		tempd22 = tempd221;
	else
		tempd22 = tempd222;
	if (fabs(tempd11) < fabs(tempd12))
		tempd1 = tempd11;
	else
		tempd1 = tempd12;
	if (fabs(tempd21) < fabs(tempd22) )
		tempd2 = tempd21;
	else
		tempd2 = tempd22;
	if (fabs(tempd1) < fabs(tempd2))
		smallestzz = tempd1;
	else
		smallestzz = tempd2;


	return (dist[indy][indx][indz] * (H-lowx) * (H-lowy) * (H-lowz) + \
		dist[indy+1][indx][indz] * (H-lowx) * lowy * (H-lowz) + \
		dist[indy][indx+1][indz] * (lowx) * (H-lowy) *(H-lowz)+ \
		dist[indy+1][indx+1][indz] * lowx * lowy * (H-lowz) + \
		dist[indy][indx][indz+1] * (H-lowx) * (H-lowy) * lowz + \
		dist[indy+1][indx][indz+1] * (H-lowx) * lowy * lowz + \
		dist[indy][indx+1][indz+1] * (lowx) * (H-lowy) * lowz+ \
		dist[indy+1][indx+1][indz+1] * lowx * lowy * lowz )/(H*H*H) - \
		0.5 * ( smallestxx * (H-lowx) * lowx + smallestyy * (H-lowy) * lowy + smallestzz * (H-lowz) * lowz );

	return (dist[indy][indx][indz] * (H-lowx) * (H-lowy) * (H-lowz) + \
		dist[indy+1][indx][indz] * (H-lowx) * lowy * (H-lowz) + \
		dist[indy][indx+1][indz] * (lowx) * (H-lowy) *(H-lowz)+ \
		dist[indy+1][indx+1][indz] * lowx * lowy * (H-lowz) + \
		dist[indy][indx][indz+1] * (H-lowx) * (H-lowy) * lowz + \
		dist[indy+1][indx][indz+1] * (H-lowx) * lowy * lowz + \
		dist[indy][indx+1][indz+1] * (lowx) * (H-lowy) * lowz+ \
		dist[indy+1][indx+1][indz+1] * lowx * lowy * lowz )/(H*H*H);

}


int cal_BothCurvs(double x, double y, double z, double H, double ***mcurv, double ***gcurv, double *calmcurv, double *calgcurv)
{
	int i, j, k;
	double lowx, lowy, lowz;
	double hcubed;

	hcubed = H*H*H;

	j = (int)((x + 1e-19 - INT_L)/H);
	i = (int)((y + 1e-19 - INT_L)/H);
	k = (int)((z + 1e-19 - INT_L)/H);

	lowx = x - INT_L - (double)(j) * H;
	lowy = y - INT_L - (double)(i) * H;
	lowz = z - INT_L - (double)(k) * H;


	*calmcurv = ( mcurv[i][j][k] * (H - lowx) * (H - lowy) * (H - lowz) + \
			mcurv[i+1][j][k] * (H - lowx) * lowy * (H - lowz) + \
			mcurv[i][j+1][k] * lowx * (H - lowy) * (H - lowz) + \
			mcurv[i+1][j+1][k] * lowx * lowy * (H - lowz) + \
			mcurv[i][j][k+1] * (H - lowx) * (H - lowy) * lowz + \
			mcurv[i+1][j][k+1] * (H - lowx) * lowy * lowz + \
			mcurv[i][j+1][k+1] * lowx * (H - lowy) * lowz + \
			mcurv[i+1][j+1][k+1] * lowx * lowy * lowz ) / hcubed;

	*calgcurv = ( gcurv[i][j][k] * (H - lowx) * (H - lowy) * (H - lowz) + \
			gcurv[i+1][j][k] * (H - lowx) * lowy * (H - lowz) + \
			gcurv[i][j+1][k] * lowx * (H - lowy) * (H - lowz) + \
			gcurv[i+1][j+1][k] * lowx * lowy * (H - lowz) + \
			gcurv[i][j][k+1] * (H - lowx) * (H - lowy) * lowz + \
			gcurv[i+1][j][k+1] * (H - lowx) * lowy * lowz + \
			gcurv[i][j+1][k+1] * lowx * (H - lowy) * lowz + \
			gcurv[i+1][j+1][k+1] * lowx * lowy * lowz ) / hcubed;
	return 0;
	
}
