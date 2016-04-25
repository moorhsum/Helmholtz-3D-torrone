
#include "weno.h"


////////////////////////////////////////////////////////////////////////////////////////////
//	
//	This computes the derivative using WENO-3 or WENO-5 method.
//	
//	len is the total size of the array in each dimension (N), dx = H = (INT_U-INT_L)/N
//	ux_p stores the one-sided derivative from the positive side
//	ux_n stores the one-sided derivative from the negative side
//	uvalue is the function value
//	i, j denotes the position where we compute the derivative
//	order denotes the method type (WENO-3 or WENO-5)
//
////////////////////////////////////////////////////////////////////////////////////////////



/*
void WENO_DX(double **uvalue, const int i, const int j, int N, double H, double *ux_p, double *ux_n, const int order)
{
        static int weno_order = DEFAULT_WENO_ORDER;
//        static int m_shift = 3; //(weno_order==5?3:2);
        static int counter = 0;
	static int flag = 0;
	static double f_value1[5], f_value2[7];
//	POINT thispt;

	counter = 0;
	weno_order = order;

	if (weno_order == 5)
	{
		for (int xidx = j-3; xidx <= j+3; xidx++)
		{
			if (xidx < 0)
				f_value2[counter] = uvalue[i][-1 * xidx];
			else if (xidx > N)
				f_value2[counter] = uvalue[i][2*N-xidx];
			else
				f_value2[counter] = uvalue[i][xidx];

			counter++;
		}
	        weno(ux_p, ux_n, f_value2, H, weno_order);
	}
	else if (weno_order == 3)
	{
	        for(int xidx = j-2; xidx <= j+2; xidx++)
        	{
			if (xidx < 0)
				f_value1[counter] = uvalue[i][-1*xidx];
			else if (xidx > N)
				f_value1[counter] = uvalue[i][2*N-xidx];
			else
				f_value1[counter] = uvalue[i][xidx];

			counter++;
        	}
//		if (flag < 5)
//		{
//			printf("H = %lf\n", H);
//			flag++;
//		}
	        weno(ux_p, ux_n, f_value1, H, weno_order);
	}
}



void WENO_DY(double **uvalue, const int i, const int j, int N, double H, double *ux_p, double *ux_n, const int order)
{
        static int weno_order = DEFAULT_WENO_ORDER;
//        static int m_shift = 3; //(weno_order==5?3:2);
        static int counter = 0;
	static double f_value1[5], f_value2[7];

	counter = 0;
	weno_order = order;

	if (weno_order == 5)
	{
		for (int xidx = i-3; xidx <= i+3; xidx++)
		{
			if (xidx < 0)
				f_value2[counter] = uvalue[-1 * xidx][j];
			else if (xidx > N)
				f_value2[counter] = uvalue[2*N-xidx][j];
			else
				f_value2[counter] = uvalue[xidx][j];

			counter++;
		}
	        weno(ux_p, ux_n, f_value2, H, weno_order);
	}
	else if (weno_order == 3)
	{
	        for(int xidx = i-2; xidx <= i+2; xidx++)
        	{
			if (xidx < 0)
				f_value1[counter] = uvalue[-1 * xidx][j];
			else if (xidx > N)
				f_value1[counter] = uvalue[2*N-xidx][j];
			else
				f_value1[counter] = uvalue[xidx][j];

			counter++;
        	}
//		printf("Here, f_value = {%lf, %lf ,%lf, %lf, %lf}\n", f_value1[0], f_value1[1], f_value1[2], f_value1[3], f_value1[4]);
	        weno(ux_p, ux_n, f_value1, H, weno_order);
//		printf("ux_p = %lf, ux_n = %lf\n", *ux_p, *ux_n);
	}
}
*/

void WENO_DX(double ***uvalue, const int i, const int j, const int k, const int N, const double H, double *ux_p, double *ux_n, const int order)
{
        int weno_order = DEFAULT_WENO_ORDER;
//        static int m_shift = 3; //(weno_order==5?3:2);
        int counter = 0;
	double f_value1[5], f_value2[7];
	int xidx;
//	POINT thispt;

	counter = 0;
	weno_order = order;

	if (weno_order == 5)
	{
		for (xidx = j-3; xidx <= j+3; xidx++)
		{
			if (xidx < 0)
				f_value2[counter] = uvalue[i][-1 * xidx][k];
			else if (xidx > N)
				f_value2[counter] = uvalue[i][2*N-xidx][k];
			else
				f_value2[counter] = uvalue[i][xidx][k];

			counter++;
		}
	        weno(ux_p, ux_n, f_value2, H, weno_order);
	}
	else if (weno_order == 3)
	{
	        for( xidx = j-2; xidx <= j+2; xidx++)
        	{
			if (xidx < 0)
				f_value1[counter] = uvalue[i][-1*xidx][k];
			else if (xidx > N)
				f_value1[counter] = uvalue[i][2*N-xidx][k];
			else
				f_value1[counter] = uvalue[i][xidx][k];

			counter++;
        	}
	        weno(ux_p, ux_n, f_value1, H, weno_order);
	}
}



void WENO_DY(double ***uvalue, const int i, const int j, const int k, const int N, const double H, double *ux_p, double *ux_n, const int order)
{
        int weno_order = DEFAULT_WENO_ORDER;
//        static int m_shift = 3; //(weno_order==5?3:2);
        int counter = 0;
	double f_value1[5], f_value2[7];
	int xidx;

	counter = 0;
	weno_order = order;

	if (weno_order == 5)
	{
		for (xidx = i-3; xidx <= i+3; xidx++)
		{
			if (xidx < 0)
				f_value2[counter] = uvalue[-1 * xidx][j][k];
			else if (xidx > N)
				f_value2[counter] = uvalue[2*N-xidx][j][k];
			else
				f_value2[counter] = uvalue[xidx][j][k];

			counter++;
		}
	        weno(ux_p, ux_n, f_value2, H, weno_order);
	}
	else if (weno_order == 3)
	{
	        for( xidx = i-2; xidx <= i+2; xidx++)
        	{
			if (xidx < 0)
				f_value1[counter] = uvalue[-1 * xidx][j][k];
			else if (xidx > N)
				f_value1[counter] = uvalue[2*N-xidx][j][k];
			else
				f_value1[counter] = uvalue[xidx][j][k];

			counter++;
        	}
//		printf("Here, f_value = {%lf, %lf ,%lf, %lf, %lf}\n", f_value1[0], f_value1[1], f_value1[2], f_value1[3], f_value1[4]);
	        weno(ux_p, ux_n, f_value1, H, weno_order);
//		printf("ux_p = %lf, ux_n = %lf\n", *ux_p, *ux_n);
	}
}




void WENO_DZ(double ***uvalue, const int i, const int j, const int k, const int N, const double H, double *ux_p, double *ux_n, const int order)
{
        int weno_order = DEFAULT_WENO_ORDER;
//        static int m_shift = 3; //(weno_order==5?3:2);
        int counter = 0;
	int xidx;
	double f_value1[5], f_value2[7];

	counter = 0;
	weno_order = order;

	if (weno_order == 5)
	{
		for (xidx = k-3; xidx <= k+3; xidx++)
		{
			if (xidx < 0)
				f_value2[counter] = uvalue[i][j][-1 * xidx];
			else if (xidx > N)
				f_value2[counter] = uvalue[i][j][2*N-xidx];
			else
				f_value2[counter] = uvalue[i][j][xidx];

			counter++;
		}
	        weno(ux_p, ux_n, f_value2, H, weno_order);
	}
	else if (weno_order == 3)
	{
	        for( xidx = k-2; xidx <= k+2; xidx++)
        	{
			if (xidx < 0)
				f_value1[counter] = uvalue[i][j][-1 * xidx];
			else if (xidx > N)
				f_value1[counter] = uvalue[i][j][2*N-xidx];
			else
				f_value1[counter] = uvalue[i][j][xidx];

			counter++;
        	}
//		printf("Here, f_value = {%lf, %lf ,%lf, %lf, %lf}\n", f_value1[0], f_value1[1], f_value1[2], f_value1[3], f_value1[4]);
	        weno(ux_p, ux_n, f_value1, H, weno_order);
//		printf("ux_p = %lf, ux_n = %lf\n", *ux_p, *ux_n);
	}
}




void weno(double *dfp, double *dfn, double fi[], const double dx, const int order)
{
  double aa, bb, cc, dd, ee, ff, df, ih;
  double IS0, IS1, IS2, a0, a1, a2, w0, w1, w2;
//	const double eps = 1.e-6;
	const double eps = pow(10.0, -12);	//	1.e-12;
  int i;

  ih  = 1. / dx;


  if ( order == 3 ) {                // WENO-3
	  i=2;
    aa = fi[i-1] - fi[i-2];
    bb = fi[i]   - fi[i-1];
    cc = fi[i+1] - fi[i];
    dd = fi[i+2] - fi[i+1]; 

    df = bb + cc;

    a0 = bb - aa;
    a1 = cc - bb;
    a2 = dd - cc;

    IS0 = ( eps + a2*a2 ) / (eps + a1*a1);
    w0  = 1. / (1. + 2.*IS0*IS0);

//	printf("Is it here?\n");
    *dfp = .5 * ih * ( df - w0*(dd - 2.*cc + bb) );
//	printf("OR HERE?\n");

    IS0 = ( eps + a0*a0 ) / (eps + a1*a1);
    w0  = 1 / (1 + 2.*IS0*IS0);

    *dfn = .5 * ih * ( df - w0*(aa - 2.*bb + cc) );
  }

  else {                             // WENO-5
	  i=3;
    aa = (fi[i-2] - fi[i-3]) * ih;
    bb = (fi[i-1] - fi[i-2]) * ih;
    cc = (fi[i]   - fi[i-1]) * ih;
    dd = (fi[i+1] - fi[i]  ) * ih;
    ee = (fi[i+2] - fi[i+1]) * ih;
    ff = (fi[i+3] - fi[i+2]) * ih;

    df = (-bb + 7.* (cc + dd) - ee) / 12.;

    aa = bb - aa;
    bb = cc - bb;
    cc = dd - cc;
    dd = ee - dd;
    ee = ff - ee;

    //  dfp = df + fweno(ee, dd, cc, bb);
    //  dfn = df - fweno(aa, bb, cc, dd);

    IS0 = 13. * (ee-dd) * (ee-dd) + 3. * (ee-3.*dd) * (ee-3.*dd);
    IS1 = 13. * (dd-cc) * (dd-cc) + 3. * (dd +  cc) * (dd +  cc);
    IS2 = 13. * (cc-bb) * (cc-bb) + 3. * (3.*cc-bb) * (3.*cc-bb);
    
    a0 = 1. / ( (eps + IS0) * (eps + IS0) );
    a1 = 6. / ( (eps + IS1) * (eps + IS1) );
    a2 = 3. / ( (eps + IS2) * (eps + IS2) );

    w1 = 1. / ( a0 + a1 + a2 );
    w0 = a0 * w1;
    w2 = a2 * w1;

    *dfp = df + w0 * ( ee - 2.*dd + cc ) / 3.  + ( w2 - .5 ) * ( dd - 2.*cc + bb ) / 6.;

    IS0 = 13. * (aa-bb) * (aa-bb) + 3. * (aa-3.*bb) * (aa-3.*bb);
    IS1 = 13. * (bb-cc) * (bb-cc) + 3. * (bb +  cc) * (bb +  cc);
    IS2 = 13. * (cc-dd) * (cc-dd) + 3. * (3.*cc-dd) * (3.*cc-dd);

    a0 = 1. / ( (eps + IS0) * (eps + IS0) );
    a1 = 6. / ( (eps + IS1) * (eps + IS1) );
    a2 = 3. / ( (eps + IS2) * (eps + IS2) );

    w1 = 1. / ( a0 + a1 + a2 );
    w0 = a0 * w1;
    w2 = a2 * w1;

    *dfn = df - w0 * ( aa - 2.*bb + cc ) / 3.  - ( w2 - .5 ) * ( bb - 2.*cc + dd ) / 6.;
  }

}



