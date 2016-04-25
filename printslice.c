
#include "printslice.h"


int print_apruz(int size, double ***level, double ***leveldx, double ***leveldy, double ***leveldz, \
		double ***levelmcurv, double ***levelgcurv)
{
	int i, j, k;
	FILE *fpapru, *fpapruimg, *fpx, *fpximg, *fpinwave, *fpinwaveimg;
	double x[size], ximg[size];
	double zx, zy, tempwx, tempwy, tempwz, tstarx, tstary, tstarz, partialx, partialy, partialz;
	double partialximg, partialyimg, partialzimg;
	double wx[size], wy[size], wz[size], wstarx[size], wstary[size], wstarz[size];
	double J[size], delta[size], gradx[size], grady[size], gradz[size], result[size];

	double thisdist;

	double sum, sumimg, pfunction, temp, tempimg;
	double epsilon = FACTOR * H;

	double nonhomogeneous[2] = {0.0, 0.0};

	int counter = 0;

	if ((fpapru = fopen( "apru.txt", "w+"))==NULL)
	{
		printf("apru.txt error.\n");
		exit(0);
	}

	if ((fpx = fopen( "large3Dx.txt", "r"))==NULL)
	{
		printf("largex.txt error.\n");
		exit(0);
	}

	if (LAPLACE==0)
	{
		if ((fpapruimg = fopen( "apruimg.txt", "w+"))==NULL)
		{
			printf("apruimg.txt error.\n");
			exit(0);
		}
		if ((fpximg = fopen("large3Dximg.txt", "r"))==NULL)
		{
			printf("largeximg.txt error.\n");
			exit(0);
		}
		if ((fpinwave = fopen("inwave.txt", "w+"))==NULL)
		{
			printf("Open inwave.txt error.\n");
			exit(0);
		}
		if ((fpinwaveimg = fopen("inwaveimg.txt", "w+"))==NULL)
		{
			printf("Open inwaveimg.txt error.\n");
			exit(0);
		}
	}


	for ( i = 0; i < size; i++)
	{
		fscanf(fpx, "%lf", &x[i]);
		if (LAPLACE==0)
			fscanf(fpximg, "%lf", &ximg[i]);
		if (i < 0)
		{
			printf("here in apru.\n");
			printf("x[%d] = %lf, ", i, x[i]);
			printf("ximg[%d] = %lf\n", i, ximg[i]);
		}
	//	printf("%lf ", x[i]);	
	}

	for ( i = 0; i < SIZE; i++)
	{
		tempwz = INT_L + H*(double)(i);
		for ( j = 0; j < SIZE; j++)
		{
			tempwy = INT_L + H*(double)(j);
			for ( k = 0; k < SIZE; k++)
			{
				tempwx = INT_L + H*(double)(k);
//				printf("distance = %lf, epsilon = %lf.\n", distance(tempwx, tempwy, tempwz), epsilon);
				thisdist = level[i][j][k];				
				if ( fabs(thisdist) > epsilon )
					continue;

//				printf("zx = %lf, zy = %lf, zz = %lf.\n", tempwx, tempwy, SLICE);
//				tstarx = starx(tempwx, tempwy, tempwz);
//				printf("counter = %d.\n", counter);
//				tstary = stary(tempwx, tempwy, tempwz);
//				printf("counter = %d.\n", counter);
//				tstarz = starz(tempwx, tempwy, tempwz);
				tstarx = tempwx - thisdist * leveldx[i][j][k];
				tstary = tempwy - thisdist * leveldy[i][j][k];
				tstarz = tempwz - thisdist * leveldz[i][j][k];
				wx[counter] = tempwx;
				wy[counter] = tempwy;
				wz[counter] = tempwz;
				wstarx[counter] = tstarx;
				wstary[counter] = tstary;
				wstarz[counter] = tstarz;


				delta[counter] = cal_delta(epsilon, thisdist);
//				J[counter] = cal_J(J_M, tempwx, tempwy, tempwz);
				J[counter] = 1. + 2.*thisdist*levelmcurv[i][j][k] + thisdist*thisdist*levelgcurv[i][j][k];

//				gradx[counter] = cal_interpgrad(1, tstarx, tstary, tstarz);
//				grady[counter] = cal_interpgrad(2, tstarx, tstary, tstarz);
//				gradz[counter] = cal_interpgrad(3, tstarx, tstary, tstarz);
				cal_interpgrad(N, H, tstarx, tstary, tstarz, leveldx, leveldy, leveldz, \
					&gradx[counter], &grady[counter], &gradz[counter]);
				result[counter] = HCUBED * delta[counter] * J[counter];

				if (result[counter] != result[counter])
				{
					printf("counter = %d, delta = %lf.\n", counter, delta[counter]);
					exit(0);
				}

				counter++;
			}
		}
	}
	

//	printf("here.\n");
	if (counter!=size)
	{
		printf("Apru size wrong.\n");
		exit(0);
	}
	counter = 0;
	for ( i = 0; i < N2; i++)
	{
		zy = INT_L + H2 * (double)(i);
		for ( j = 0; j < N2; j++)
		{
			zx = INT_L + H2 * (double)(j);
			sum = 0.0;
			sumimg = 0.0;
//			counter = 0;
			if (LAPLACE == 0)
			{
				if ( ((zx > 0.0) && (fabs(zy) < radius)) || (distance(zx, zy, ZSLICE) > -0.06) )
				{
					fprintf(fpinwave, "%.8lf ", 0.0);
					fprintf(fpinwaveimg, "%.8lf ", 0.0);
				}
				else
				{
					fprintf(fpinwave, "%.8lf ", cos(WAVE * zx));
					fprintf(fpinwaveimg, "%.8lf ", sin(WAVE * zx));
				}

				if (FORCING == 1)
				{
					if (distance(zx,zy, ZSLICE) > -0.06)
					{
						fprintf(fpapru, "%.8lf ", sum);
						fprintf(fpapruimg, "%.8lf ", sum);
						continue;
					}
				}
			}

			for ( k = 0; k < size; k++)
			{
				tstarx = wstarx[k];
				tstary = wstary[k];
				tstarz = wstarz[k];

				partialx = cal_partial(1, zx, zy, ZSLICE, tstarx, tstary, tstarz);
				partialy = cal_partial(2, zx, zy, ZSLICE, tstarx, tstary, tstarz);
				partialz = cal_partial(3, zx, zy, ZSLICE, tstarx, tstary, tstarz);
				partialximg = cal_partial(4, zx, zy, ZSLICE, tstarx, tstary, tstarz);
				partialyimg = cal_partial(5, zx, zy, ZSLICE, tstarx, tstary, tstarz);
				partialzimg = cal_partial(6, zx, zy, ZSLICE, tstarx, tstary, tstarz);
			
			
				//	for interior problem	//	
				if (INTERIOR_MODE==1)
				{
					temp = (-1.0) * (partialx*gradx[k] + partialy*grady[k] + partialz*gradz[k]); // outer normal
					if (LAPLACE==0)
						tempimg = (-1.0) * (partialximg*gradx[k] + partialyimg*grady[k] + partialzimg*gradz[k]);
				}
				//	for exterior problem	//
				else
				{
					temp = partialx*gradx[k] + partialy*grady[k] + partialz*gradz[k];		// inner normal
					if (LAPLACE==0)
						tempimg = partialximg*gradx[k] + partialyimg*grady[k] + partialzimg*gradz[k];
				}
				if (LAPLACE==1)
				{
					sum += (temp * result[k] * x[k]);
				}
				else
				{
					sum += (result[k] * (temp * x[k] - tempimg * ximg[k]));
					sumimg += (result[k] * (tempimg * x[k] + temp * ximg[k]));
				}

				if ((temp!=temp)||(temp+1.0==temp)||(tempimg!=tempimg)||(tempimg+1.0==tempimg))
				{
					printf("zx = %lf, zy = %lf, zz = %lf, counter = %d, ", zx, zy, ZSLICE, counter);
					printf("wx = %lf, wy = %lf, wstarx = %lf, wstary = %lf.\n", wx[k], wy[k], wstarx[k], wstary[k]);
					printf("delta = %lf, pfunction = %lf, temp = %lf\n", delta[k], pfunction, temp);
				}
//				if (counter < 10)
//					printf("%lf", sum);
//				counter++;

			}

			if (NONHOM_MODE == 1)
			{
				cal_nonhom(zx, zy, ZSLICE, nonhomogeneous);
				sum += nonhomogeneous[0];
				sumimg += nonhomogeneous[1];
			}
			fprintf(fpapru, "%.8lf ", sum);
			if (LAPLACE==0)
				fprintf(fpapruimg, "%.8lf ", sumimg);

		}
		fprintf(fpapru, "\n");
		if (LAPLACE == 0)
		{
			fprintf(fpinwave, "\n");
			fprintf(fpinwaveimg, "\n");
			fprintf(fpapruimg, "\n");
		}
	}

	fclose(fpapru);
	fclose(fpx);
	if (LAPLACE==0)
	{
		fclose(fpinwave);
		fclose(fpximg);
		fclose(fpinwaveimg);
		fclose(fpapruimg);
	}
	return 0;
}


int print_apruy(int size, double ***level, double ***leveldx, double ***leveldy, double ***leveldz, \
		double ***levelmcurv, double ***levelgcurv)
{
	int i, j, k;
	FILE *fpapruy, *fpapruyimg, *fpx, *fpximg, *fpinwave, *fpinwaveimg;
	double x[size], ximg[size];
	double zx, zz, tempwx, tempwy, tempwz, tstarx, tstary, tstarz, partialx, partialy, partialz;
	double partialximg, partialyimg, partialzimg;
	double wx[size], wy[size], wz[size], wstarx[size], wstary[size], wstarz[size];
	double J[size], delta[size], gradx[size], grady[size], gradz[size], result[size];

	double sum, sumimg, pfunction, temp, tempimg;
	double epsilon = FACTOR * H;
	double thisdist;

	double nonhomogeneous[2] = {0.0, 0.0};

	int counter = 0;

	if ((fpapruy = fopen( "apruyslice.txt", "w+"))==NULL)
	{
		printf("apruyslice.txt error.\n");
		exit(0);
	}

	if ((fpx = fopen( "large3Dx.txt", "r"))==NULL)
	{
		printf("largex.txt error.\n");
		exit(0);
	}

	if (LAPLACE==0)
	{
		if ((fpapruyimg = fopen( "apruysliceimg.txt", "w+"))==NULL)
		{
			printf("apruysliceimg.txt error.\n");
			exit(0);
		}
		if ((fpximg = fopen("large3Dximg.txt", "r"))==NULL)
		{
			printf("largeximg.txt error.\n");
			exit(0);
		}
		if ((fpinwave = fopen("inwave.txt", "w+"))==NULL)
		{
			printf("Open inwave.txt error.\n");
			exit(0);
		}
		if ((fpinwaveimg = fopen("inwaveimg.txt", "w+"))==NULL)
		{
			printf("Open inwaveimg.txt error.\n");
			exit(0);
		}
	}


	for ( i = 0; i < size; i++)
	{
		fscanf(fpx, "%lf", &x[i]);
		if (LAPLACE==0)
			fscanf(fpximg, "%lf", &ximg[i]);
		if (i < 0)
		{
			printf("here in apru.\n");
			printf("x[%d] = %lf, ", i, x[i]);
			printf("ximg[%d] = %lf\n", i, ximg[i]);
		}
	//	printf("%lf ", x[i]);	
	}

	for ( i = 0; i < SIZE; i++)
	{
		tempwz = INT_L + H*(double)(i);
		for ( j = 0; j < SIZE; j++)
		{
			tempwy = INT_L + H*(double)(j);
			for ( k = 0; k < SIZE; k++)
			{
				tempwx = INT_L + H*(double)(k);
				thisdist = level[i][j][k];
//				printf("distance = %lf, epsilon = %lf.\n", distance(tempwx, tempwy, tempwz), epsilon);
				if ( fabs(thisdist) > epsilon )
					continue;

//				printf("zx = %lf, zy = %lf, zz = %lf.\n", tempwx, tempwy, SLICE);
//				tstarx = starx(tempwx, tempwy, tempwz);
////				printf("counter = %d.\n", counter);
//				tstary = stary(tempwx, tempwy, tempwz);
//				printf("counter = %d.\n", counter);
//				tstarz = starz(tempwx, tempwy, tempwz);
				tstarx = tempwx - thisdist * leveldx[i][j][k];
				tstary = tempwy - thisdist * leveldy[i][j][k];
				tstarz = tempwz - thisdist * leveldz[i][j][k];
				wx[counter] = tempwx;
				wy[counter] = tempwy;
				wz[counter] = tempwz;
				wstarx[counter] = tstarx;
				wstary[counter] = tstary;
				wstarz[counter] = tstarz;


				delta[counter] = cal_delta(epsilon, distance(tempwx, tempwy, tempwz));
				J[counter] = 1. + 2.*thisdist*levelmcurv[i][j][k] + thisdist*thisdist*levelgcurv[i][j][k];

//				gradx[counter] = cal_interpgrad(1, tstarx, tstary, tstarz);
//				grady[counter] = cal_interpgrad(2, tstarx, tstary, tstarz);
//				gradz[counter] = cal_interpgrad(3, tstarx, tstary, tstarz);
				cal_interpgrad(N, H, tstarx, tstary, tstarz, leveldx, leveldy, leveldz, \
					&gradx[counter], &grady[counter], &gradz[counter]);
				result[counter] = HCUBED * delta[counter] * J[counter];

				if (result[counter] != result[counter])
				{
					printf("counter = %d, delta = %lf.\n", counter, delta[counter]);
					exit(0);
				}

				counter++;
			}
		}
	}
	

//	printf("here.\n");
	if (counter!=size)
	{
		printf("Apru size wrong.\n");
		exit(0);
	}
	counter = 0;
	for ( i = 0; i < N2; i++)
	{
		zz = INT_L + H2 * (double)(i);
		for ( j = 0; j < N2; j++)
		{
			zx = INT_L + H2 * (double)(j);
			sum = 0.0;
			sumimg = 0.0;
//			counter = 0;
			if (LAPLACE == 0)
			{
				if ( ((zx > 0.0) && (fabs(YSLICE) < radius)) || (distance(zx, YSLICE, zz) > -0.06) )
				{
					fprintf(fpinwave, "%.8lf ", 0.0);
					fprintf(fpinwaveimg, "%.8lf ", 0.0);
				}
				else
				{
					fprintf(fpinwave, "%.8lf ", cos(WAVE * zx));
					fprintf(fpinwaveimg, "%.8lf ", sin(WAVE * zx));
				}

				if (FORCING == 1)
				{
					if (distance(zx, YSLICE, zz) > -0.06)
					{
						fprintf(fpapruy, "%.8lf ", sum);
						fprintf(fpapruyimg, "%.8lf ", sum);
						continue;
					}
				}
			}

			for ( k = 0; k < size; k++)
			{
				tstarx = wstarx[k];
				tstary = wstary[k];
				tstarz = wstarz[k];

				partialx = cal_partial(1, zx, YSLICE, zz, tstarx, tstary, tstarz);
				partialy = cal_partial(2, zx, YSLICE, zz, tstarx, tstary, tstarz);
				partialz = cal_partial(3, zx, YSLICE, zz, tstarx, tstary, tstarz);
				partialximg = cal_partial(4, zx, YSLICE, zz, tstarx, tstary, tstarz);
				partialyimg = cal_partial(5, zx, YSLICE, zz, tstarx, tstary, tstarz);
				partialzimg = cal_partial(6, zx, YSLICE, zz, tstarx, tstary, tstarz);
			
			
				//	for interior problem	//	
				if (INTERIOR_MODE==1)
				{
					temp = (-1.0) * (partialx*gradx[k] + partialy*grady[k] + partialz*gradz[k]); // outer normal
					if (LAPLACE==0)
						tempimg = (-1.0) * (partialximg*gradx[k] + partialyimg*grady[k] + partialzimg*gradz[k]);
				}
				//	for exterior problem	//
				else
				{
					temp = partialx*gradx[k] + partialy*grady[k] + partialz*gradz[k];		// inner normal
					if (LAPLACE==0)
						tempimg = partialximg*gradx[k] + partialyimg*grady[k] + partialzimg*gradz[k];
				}
				if (LAPLACE==1)
				{
					sum += (temp * result[k] * x[k]);
				}
				else
				{
					sum += (result[k] * (temp * x[k] - tempimg * ximg[k]));
					sumimg += (result[k] * (tempimg * x[k] + temp * ximg[k]));
				}

				if ((temp!=temp)||(temp+1.0==temp)||(tempimg!=tempimg)||(tempimg+1.0==tempimg))
				{
					printf("zx = %lf, zy = %lf, zz = %lf, counter = %d, ", zx, YSLICE, zz, counter);
					printf("wx = %lf, wy = %lf, wstarx = %lf, wstary = %lf.\n", wx[k], wy[k], wstarx[k], wstary[k]);
					printf("delta = %lf, pfunction = %lf, temp = %lf\n", delta[k], pfunction, temp);
				}
//				if (counter < 10)
//					printf("%lf", sum);
//				counter++;

			}

			if (NONHOM_MODE == 1)
			{
				cal_nonhom(zx, YSLICE, zz, nonhomogeneous);
				sum += nonhomogeneous[0];
				sumimg += nonhomogeneous[1];
			}
			fprintf(fpapruy, "%.8lf ", sum);
			if (LAPLACE==0)
				fprintf(fpapruyimg, "%.8lf ", sumimg);

		}
		fprintf(fpapruy, "\n");
		if (LAPLACE == 0)
		{
			fprintf(fpinwave, "\n");
			fprintf(fpinwaveimg, "\n");
			fprintf(fpapruyimg, "\n");
		}
	}

	fclose(fpapruy);
	fclose(fpx);
	if (LAPLACE==0)
	{
		fclose(fpinwave);
		fclose(fpximg);
		fclose(fpinwaveimg);
		fclose(fpapruyimg);
	}
	return 0;
}


int print_aprux(int size, double ***level, double ***leveldx, double ***leveldy, double ***leveldz, \
		double ***levelmcurv, double ***levelgcurv)
{
	int i, j, k;
	FILE *fpaprux, *fpapruximg, *fpx, *fpximg, *fpinwave, *fpinwaveimg;
	double x[size], ximg[size];
	double zz, zy, tempwx, tempwy, tempwz, tstarx, tstary, tstarz, partialx, partialy, partialz;
	double partialximg, partialyimg, partialzimg;
	double wx[size], wy[size], wz[size], wstarx[size], wstary[size], wstarz[size];
	double J[size], delta[size], gradx[size], grady[size], gradz[size], result[size];

	double sum, sumimg, pfunction, temp, tempimg;
	double epsilon = FACTOR * H;
	double thisdist;

	double nonhomogeneous[2] = {0.0, 0.0};

	int counter = 0;

	if ((fpaprux = fopen( "apruxslice.txt", "w+"))==NULL)
	{
		printf("apruxslice.txt error.\n");
		exit(0);
	}

	if ((fpx = fopen( "large3Dx.txt", "r"))==NULL)
	{
		printf("largex.txt error.\n");
		exit(0);
	}

	if (LAPLACE==0)
	{
		if ((fpapruximg = fopen( "apruxsliceimg.txt", "w+"))==NULL)
		{
			printf("apruxsliceimg.txt error.\n");
			exit(0);
		}
		if ((fpximg = fopen("large3Dximg.txt", "r"))==NULL)
		{
			printf("largeximg.txt error.\n");
			exit(0);
		}
		if ((fpinwave = fopen("inwave.txt", "w+"))==NULL)
		{
			printf("Open inwave.txt error.\n");
			exit(0);
		}
		if ((fpinwaveimg = fopen("inwaveimg.txt", "w+"))==NULL)
		{
			printf("Open inwaveimg.txt error.\n");
			exit(0);
		}
	}


	for ( i = 0; i < size; i++)
	{
		fscanf(fpx, "%lf", &x[i]);
		if (LAPLACE==0)
			fscanf(fpximg, "%lf", &ximg[i]);
		if (i < 0)
		{
			printf("here in apru.\n");
			printf("x[%d] = %lf, ", i, x[i]);
			printf("ximg[%d] = %lf\n", i, ximg[i]);
		}
	//	printf("%lf ", x[i]);	
	}

	for ( i = 0; i < SIZE; i++)
	{
		tempwz = INT_L + H*(double)(i);
		for ( j = 0; j < SIZE; j++)
		{
			tempwy = INT_L + H*(double)(j);
			for ( k = 0; k < SIZE; k++)
			{
				tempwx = INT_L + H*(double)(k);
				thisdist = level[i][j][k];
//				printf("distance = %lf, epsilon = %lf.\n", distance(tempwx, tempwy, tempwz), epsilon);
				if ( fabs(thisdist) > epsilon )
					continue;

//				printf("zx = %lf, zy = %lf, zz = %lf.\n", tempwx, tempwy, SLICE);
//				tstarx = starx(tempwx, tempwy, tempwz);
////				printf("counter = %d.\n", counter);
//				tstary = stary(tempwx, tempwy, tempwz);
//				printf("counter = %d.\n", counter);
//				tstarz = starz(tempwx, tempwy, tempwz);
				tstarx = tempwx - thisdist * leveldx[i][j][k];				
				tstary = tempwy - thisdist * leveldy[i][j][k];				
				tstarz = tempwz - thisdist * leveldz[i][j][k];				
				wx[counter] = tempwx;
				wy[counter] = tempwy;
				wz[counter] = tempwz;
				wstarx[counter] = tstarx;
				wstary[counter] = tstary;
				wstarz[counter] = tstarz;


				delta[counter] = cal_delta(epsilon, distance(tempwx, tempwy, tempwz));
//				J[counter] = cal_J(J_M, tempwx, tempwy, tempwz);

				J[counter] = 1. * 2.*thisdist*levelmcurv[i][j][k] + thisdist*thisdist*levelgcurv[i][j][k];
//				gradx[counter] = cal_interpgrad(1, tstarx, tstary, tstarz);
//				grady[counter] = cal_interpgrad(2, tstarx, tstary, tstarz);
//				gradz[counter] = cal_interpgrad(3, tstarx, tstary, tstarz);
				cal_interpgrad(N, H, tstarx, tstary, tstarz, leveldx, leveldy, leveldz, \
					&gradx[counter], &grady[counter], &gradz[counter]);
				result[counter] = HCUBED * delta[counter] * J[counter];

				if (result[counter] != result[counter])
				{
					printf("counter = %d, delta = %lf.\n", counter, delta[counter]);
					exit(0);
				}

				counter++;
			}
		}
	}
	

//	printf("here.\n");
	if (counter!=size)
	{
		printf("Apru size wrong.\n");
		exit(0);
	}
	counter = 0;
	for ( i = 0; i < N2; i++)
	{
		zz = INT_L + H2 * (double)(i);
		for ( j = 0; j < N2; j++)
		{
			zy = INT_L + H2 * (double)(j);
			sum = 0.0;
			sumimg = 0.0;
//			counter = 0;
			if (LAPLACE == 0)
			{
				if ( ((zz > 0.0) && (fabs(zy) < radius)) || (distance(XSLICE, zy, zz) > -0.06) )
				{
					fprintf(fpinwave, "%.8lf ", 0.0);
					fprintf(fpinwaveimg, "%.8lf ", 0.0);
				}
				else
				{
					fprintf(fpinwave, "%.8lf ", cos(WAVE * XSLICE));
					fprintf(fpinwaveimg, "%.8lf ", sin(WAVE * XSLICE));
				}

				if (FORCING == 1)
				{
					if (distance(XSLICE,zy,zz) > -0.06)
					{
						fprintf(fpaprux, "%.8lf ", sum);
						fprintf(fpapruximg, "%.8lf ", sum);
						continue;
					}
				}
			}

			for ( k = 0; k < size; k++)
			{
				tstarx = wstarx[k];
				tstary = wstary[k];
				tstarz = wstarz[k];

				partialx = cal_partial(1, XSLICE, zy, zz, tstarx, tstary, tstarz);
				partialy = cal_partial(2, XSLICE, zy, zz, tstarx, tstary, tstarz);
				partialz = cal_partial(3, XSLICE, zy, zz, tstarx, tstary, tstarz);
				partialximg = cal_partial(4, XSLICE, zy, zz, tstarx, tstary, tstarz);
				partialyimg = cal_partial(5, XSLICE, zy, zz, tstarx, tstary, tstarz);
				partialzimg = cal_partial(6, XSLICE, zy, zz, tstarx, tstary, tstarz);
			
			
				//	for interior problem	//	
				if (INTERIOR_MODE==1)
				{
					temp = (-1.0) * (partialx*gradx[k] + partialy*grady[k] + partialz*gradz[k]); // outer normal
					if (LAPLACE==0)
						tempimg = (-1.0) * (partialximg*gradx[k] + partialyimg*grady[k] + partialzimg*gradz[k]);
				}
				//	for exterior problem	//
				else
				{
					temp = partialx*gradx[k] + partialy*grady[k] + partialz*gradz[k];		// inner normal
					if (LAPLACE==0)
						tempimg = partialximg*gradx[k] + partialyimg*grady[k] + partialzimg*gradz[k];
				}
				if (LAPLACE==1)
				{
					sum += (temp * result[k] * x[k]);
				}
				else
				{
					sum += (result[k] * (temp * x[k] - tempimg * ximg[k]));
					sumimg += (result[k] * (tempimg * x[k] + temp * ximg[k]));
				}

				if ((temp!=temp)||(temp+1.0==temp)||(tempimg!=tempimg)||(tempimg+1.0==tempimg))
				{
					printf("zx = %lf, zy = %lf, zz = %lf, counter = %d, ", XSLICE, zy, zz, counter);
					printf("wx = %lf, wy = %lf, wstarx = %lf, wstary = %lf.\n", wx[k], wy[k], wstarx[k], wstary[k]);
					printf("delta = %lf, pfunction = %lf, temp = %lf\n", delta[k], pfunction, temp);
				}
//				if (counter < 10)
//					printf("%lf", sum);
//				counter++;

			}

			if (NONHOM_MODE == 1)
			{
				cal_nonhom(XSLICE, zy, zz, nonhomogeneous);
				sum += nonhomogeneous[0];
				sumimg += nonhomogeneous[1];
			}
			fprintf(fpaprux, "%.8lf ", sum);
			if (LAPLACE==0)
				fprintf(fpapruximg, "%.8lf ", sumimg);

		}
		fprintf(fpaprux, "\n");
		if (LAPLACE == 0)
		{
			fprintf(fpinwave, "\n");
			fprintf(fpinwaveimg, "\n");
			fprintf(fpapruximg, "\n");
		}
	}

	fclose(fpaprux);
	fclose(fpx);
	if (LAPLACE==0)
	{
		fclose(fpinwave);
		fclose(fpximg);
		fclose(fpinwaveimg);
		fclose(fpapruximg);
	}
	return 0;
}



