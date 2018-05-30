/*!
*    Copyright 2017 by ChengFeng Liu, Jin Zhu<zhuj37mail2.sysu.edu.cn>
*     This program is free software: you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation, either version 3 of the License, or
* (at your option) any later version.

* This program is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
* GNU General Public License for more details.

* You should have received a copy of the GNU General Public License
* along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifdef _OPENMP
# include "omp.h"
#endif

#include "math.h"
#include "string.h" 
#include "stdlib.h" 
#include "stdio.h"
#include "utilities.h"
#include "bcor.h"
#include "BI.h"


void Ball_Correlation(double *bcor_stat, int *n, int *p, double *x, double **Dy, int **yidx)
{
	int i, j, k, pi, src, lastpos, *yrank, *isource, *icount, *xy_index, *xy_temp, **xyidx;
	double pxy, px, py, lastval, *xx_cpy, *yy_cpy;
	double bcov_weight0 = 0.0, bcov_weight_prob = 0.0, bcov_weight_hhg = 0.0, bcov_fixed_ball = 0.0;
	double bcov_weight0_x = 0.0, bcov_weight_prob_x = 0.0;
	double bcov_weight0_y = 0.0, bcov_weight_prob_y = 0.0;
	double hhg_ball_num = 0.0, minor_ball_prop = 2 / (*n);

	double **Dx, *x_cpy;
	int **xidx;
	Dx = alloc_matrix(*n, *n);
	xidx = alloc_int_matrix(*n, *n);
	x_cpy = (double *)malloc(*n * sizeof(double));

	for (i = 0; i<*n; i++)
	{
		for (j = 0; j<*n; j++)
		{
			xidx[i][j] = j;
		}
	}
	Euclidean_distance(x, Dx, *n, *p);
	for (i = 0; i<(*n); i++)
	{
		memcpy(x_cpy, Dx[i], *n * sizeof(double));
		quicksort(x_cpy, xidx[i], 0, *n - 1);
	}
	free(x_cpy);

	yrank = (int *)malloc(*n * sizeof(int));
	isource = (int *)malloc(*n * sizeof(int));
	icount = (int *)malloc(*n * sizeof(int));
	xy_index = (int *)malloc(*n * sizeof(int));
	xy_temp = (int *)malloc(*n * sizeof(int));
	xyidx = alloc_int_matrix(*n, *n);
	xx_cpy = (double *)malloc(*n * sizeof(double));
	yy_cpy = (double *)malloc(*n * sizeof(double));

	for (i = 0; i<*n; i++)
		for (j = 0; j<*n; j++)
			xyidx[i][j] = j;

	for (i = 0; i<(*n); i++)
	{
		memcpy(xx_cpy, Dx[i], *n * sizeof(double));
		for (j = 0; j<*n; j++)
			yy_cpy[j] = Dy[i][j];
		quicksort2(xx_cpy, yy_cpy, xyidx[i], 0, *n - 1);
	}
	free(xx_cpy);
	free(yy_cpy);

	for (i = 0; i<(*n); i++)
	{
		pi = i;
		lastval = 0;
		lastpos = -1;
		for (j = *n - 1, k = *n - 1; j >= 1; --j, --k)
		{
			k -= (yidx[pi][k] == pi);
			if (lastpos == -1 || Dy[pi][yidx[pi][k]] != lastval)
			{
				lastval = Dy[pi][yidx[pi][k]];
				lastpos = j;
			}
			src = yidx[pi][k];
			src -= (src > i);
			yrank[src] = lastpos;
		}

		for (j = 0, k = 0; j < *n - 1; ++j, ++k)
		{
			k += (xyidx[i][k] == i);
			src = xyidx[i][k]; // NOTE: k may be different than from the line above
			src -= (src > i);
			xy_index[j] = yrank[src];
			isource[j] = j;
			icount[j] = 0;
			xy_temp[j] = xy_index[j];
		}
		Inversions(xy_temp, isource, icount, *n - 1, *n);
		lastval = 0;
		lastpos = -1;

		for (j = *n - 2, k = *n - 1; j >= 0; --j, --k)
		{
			k -= (xidx[i][k] == i);
			if (lastpos == -1 || Dx[i][xidx[i][k]] != lastval) {
				lastval = Dx[i][xidx[i][k]];
				lastpos = j;
			}

			pxy = lastpos - icount[j] + 2;
			px = lastpos + 2;
			py = xy_index[j] + 1;
			px /= (*n);
			py /= (*n);
			pxy /= (*n);
			// compute BCov(X, Y)
			bcov_fixed_ball = pow(pxy - px*py, 2);
			bcov_weight0 += bcov_fixed_ball;
			bcov_weight_prob += bcov_fixed_ball / (px*py);
			if (px > minor_ball_prop && py > minor_ball_prop && px != 1 && py != 1)
			{
				bcov_weight_hhg += bcov_fixed_ball / ((px - minor_ball_prop)*(1.0 - px + minor_ball_prop)*(py - minor_ball_prop)*(1.0 - py + minor_ball_prop));
				hhg_ball_num += 1;
			}
			// compute BCov(X, X)
			bcov_weight_prob_x += (1.0 - px)*(1.0 - px);
			bcov_weight0_x += px*px*(1.0 - px)*(1.0 - px);

			// compute BCov(Y, Y)
			bcov_weight_prob_y += (1.0 - py)*(1.0 - py);
			bcov_weight0_y += py*py*(1.0 - py)*(1.0 - py);

		}
		pxy = 0;
		px = 0;
		py = 0;
		for (j = 0; j<*n; j++)
		{
			if (Dx[i][xidx[i][j]] == 0)
			{
				px += 1;
				if (Dy[pi][xidx[i][j]] == 0)
				{
					pxy += 1;
					py += 1;
				}
			}
			else if (Dy[pi][xidx[i][j]] == 0)
				py += 1;
		}
		px /= (*n);
		py /= (*n);
		pxy /= (*n);
		bcov_fixed_ball = pow(pxy - px*py, 2);
		bcov_weight0 += bcov_fixed_ball;
		bcov_weight_prob += bcov_fixed_ball / (px*py);
		if (px > minor_ball_prop && py > minor_ball_prop && px != 1 && py != 1)
		{
			bcov_weight_hhg += bcov_fixed_ball / ((px - minor_ball_prop)*(1.0 - px + minor_ball_prop)*(py - minor_ball_prop)*(1.0 - py + minor_ball_prop));
			hhg_ball_num += 1;
		}
		// compute BCov(X, X)
		bcov_weight_prob_x += (1.0 - px)*(1.0 - px);
		bcov_weight0_x += px*px*(1.0 - px)*(1.0 - px);

		// compute BCov(Y, Y)
		bcov_weight_prob_y += (1.0 - py)*(1.0 - py);
		bcov_weight0_y += py*py*(1.0 - py)*(1.0 - py);
	}
	bcor_stat[0] = bcov_weight0 / (sqrt(bcov_weight0_x) * sqrt(bcov_weight0_y));
	bcor_stat[1] = bcov_weight_prob / (sqrt(bcov_weight_prob_x) * sqrt(bcov_weight_prob_y));
	bcor_stat[2] = bcov_weight_hhg / (hhg_ball_num);

	// free memory
	free_matrix(Dx, *n, *n);
	free_int_matrix(xidx, *n, *n);
	free_int_matrix(xyidx, *n, *n);
	free(isource);
	free(icount);
	free(xy_index);
	free(yrank);
	free(xy_temp);
	return;
}


void U_Ball_Correlation(double *bcor_stat, int *n, double *x, int *yrank, int **lowyidx, int **higyidx)
{
	int **Rank, **lowxidx, **higxidx, *xidx, *xrank;
	xidx = (int *)malloc(*n * sizeof(int));
	xrank = (int *)malloc(*n * sizeof(int));
	Rank = alloc_int_matrix((*n) + 1, (*n) + 1);
	lowxidx = alloc_int_matrix(*n, *n);
	higxidx = alloc_int_matrix(*n, *n);

	for (int k = 0; k<*n; k++)
	{
		xidx[k] = k;
	}
	quicksort(x, xidx, 0, *n - 1);
	ranksort(n, xrank, x, xidx);
	createidx(n, xidx, x, lowxidx, higxidx);
	initRank_bcor(*n, Rank, xrank, yrank);
	
	free(xrank);
	free(xidx);

	int i, j, pi, pj;
	double px, py, pxy;
	double bcov_weight0 = 0.0, bcov_weight_prob = 0.0, bcov_weight_hhg = 0.0, bcov_fixed_ball = 0.0;
	double bcov_weight0_x = 0.0, bcov_weight_prob_x = 0.0;
	double bcov_weight0_y = 0.0, bcov_weight_prob_y = 0.0;
	double minor_ball_prop = 2 / (*n), hhg_ball_num = 0.0;
	for (i = 0; i < *n; i++) {
		for (j = 0; j < *n; j++) {
			pi = i;
			pj = j;
			px = higxidx[i][j] - lowxidx[i][j] + 1;
			py = higyidx[pi][pj] - lowyidx[pi][pj] + 1;
			pxy = Rank[higxidx[i][j]][higyidx[pi][pj]] + Rank[lowxidx[i][j] - 1][lowyidx[pi][pj] - 1];
			pxy -= (Rank[higxidx[i][j]][lowyidx[pi][pj] - 1] + Rank[lowxidx[i][j] - 1][higyidx[pi][pj]]);

			pxy /= (*n);
			px /= (*n);
			py /= (*n);
			bcov_fixed_ball = pow(pxy - px*py, 2);
			bcov_weight0 += bcov_fixed_ball;
			bcov_weight_prob += bcov_fixed_ball / (px*py);
			if (px > minor_ball_prop && py > minor_ball_prop && px != 1 && py != 1)
			{
				bcov_weight_hhg += bcov_fixed_ball / ((px - minor_ball_prop)*(1.0 - px + minor_ball_prop)*(py - minor_ball_prop)*(1.0 - py + minor_ball_prop));
				hhg_ball_num += 1;
			}

			bcov_weight_prob_x += (1.0 - px)*(1.0 - px);
			bcov_weight0_x += px*px*(1.0 - px)*(1.0 - px);

			bcov_weight_prob_y += (1.0 - py)*(1.0 - py);
			bcov_weight0_y += py*py*(1.0 - py)*(1.0 - py);
		}
	}
	bcor_stat[0] = bcov_weight0 / (sqrt(bcov_weight0_x) * sqrt(bcov_weight0_y));
	bcor_stat[1] = bcov_weight_prob / (sqrt(bcov_weight_prob_x) * sqrt(bcov_weight_prob_y));
	bcor_stat[2] = bcov_weight_hhg / (hhg_ball_num);
	
	free_int_matrix(Rank, (*n) + 1, (*n) + 1);
	free_int_matrix(lowxidx, *n, *n);
	free_int_matrix(higxidx, *n, *n);
	return;
}


void _bcor_test(double *bcorsis_stat, double *y, double *x, int *x_number, int *f_number, int *n, int *p, int *nthread)
{
	// Start the multi-threaded support only when feature number is large
	if (*f_number > 99)
	{
#ifdef _OPENMP
		omp_set_num_threads(*nthread);
#endif
	}
	else {
		int single_thread = 1;
		*nthread = single_thread;
#ifdef _OPENMP
		omp_set_num_threads(*nthread);
#endif
	}

	int i, j, **yidx;
	double **Dy, *y_cpy;

	Dy = alloc_matrix(*n, *n);
	yidx = alloc_int_matrix(*n, *n);
	y_cpy = (double *)malloc(*n * sizeof(double));

	if (*p == 0)
	{
		vector2matrix(y, Dy, *n, *n, 1);
	}
	else {
		Euclidean_distance(y, Dy, *n, *p);
	}

	for (i = 0; i<*n; i++)
	{
		for (j = 0; j<*n; j++)
		{
			yidx[i][j] = j;
		}
	}

	// sort distance matrix y:
	for (i = 0; i<(*n); i++)
	{
		memcpy(y_cpy, Dy[i], *n * sizeof(double));
		quicksort(y_cpy, yidx[i], 0, *n - 1);
	}
	free(y_cpy);

#pragma omp parallel
	{
		int f_thread, i_thread, k_thread, s_thread, stop_index_thread, x_size_thread;
		double *x_thread;
		double bcorsis_stat_tmp[3];

		// main loop for calculate bcor statistic
#pragma omp for
			for (f_thread = 0; f_thread < *f_number; f_thread++)
			{
				// extract value to x_thread
				k_thread = 0;
				i_thread = (f_thread) * (*n);
				x_size_thread = x_number[f_thread] * (*n);
				stop_index_thread = (f_thread) * (*n) + x_size_thread;
				x_thread = (double *)malloc(x_size_thread * sizeof(double));

				while (i_thread < stop_index_thread) 
				{
					x_thread[k_thread] = x[i_thread];
					k_thread++;
					i_thread++;
				}

				Ball_Correlation(bcorsis_stat_tmp, n, &x_number[f_thread], x_thread, Dy, yidx);
				for (s_thread = 0; s_thread < 3; s_thread++)
				{
					bcorsis_stat[3 * f_thread + s_thread] = bcorsis_stat_tmp[s_thread];
				}
				free(x_thread);
			}
	}
	return;
}


void _fast_bcor_test(double *bcorsis_stat, double *y, double *x, int *f_number, int *n, int *nthread)
{
	// Start the multi-threaded support only when feature number is large
	if (*f_number > 99)
	{
#ifdef _OPENMP
		omp_set_num_threads(*nthread);
#endif
	}
	else {
		int single_thread = 1;
		*nthread = single_thread;
#ifdef _OPENMP
		omp_set_num_threads(*nthread);
#endif
	}

	int i, *yidx, *yrank, **lowyidx, **higyidx;

	yidx = (int *)malloc(*n * sizeof(int));
	yrank = (int *)malloc(*n * sizeof(int));
	lowyidx = alloc_int_matrix(*n, *n);
	higyidx = alloc_int_matrix(*n, *n);

	for (i = 0; i<*n; i++)
	{
		yidx[i] = i;
	}
	quicksort(y, yidx, 0, *n - 1);
	ranksort(n, yrank, y, yidx);
	createidx(n, yidx, y, lowyidx, higyidx);
	free(yidx);

#pragma omp parallel
	{
		int f_thread, i_thread, s_thread, start_index;
		double *x_thread;
		double bcorsis_stat_tmp[3];
		
		x_thread = (double *)malloc((*n) * sizeof(double));

		// main loop for calculate bcor statistic
#pragma omp for
			for (f_thread = 0; f_thread < *f_number; f_thread++)
			{
				// extract value to x_thread
				start_index = f_thread * (*n);
				for (i_thread = 0; i_thread < (*n); i_thread++)
				{
					x_thread[i_thread] = x[start_index + i_thread];
				}

				U_Ball_Correlation(bcorsis_stat_tmp, n, x_thread, yrank, lowyidx, higyidx);
				for (s_thread = 0; s_thread < 3; s_thread++)
				{
					bcorsis_stat[3 * f_thread + s_thread] = bcorsis_stat_tmp[s_thread];
				}
			}
		free(x_thread);
	}

	//int f_thread, i_thread, s_thread, start_index;
	//double *x_thread;
	//double bcorsis_stat_tmp[3];

	//// main loop for calculate bcor statistic
	//for (f_thread = 0; f_thread < *f_number; f_thread++)
	//{
	//	// extract value to x_thread
	//	x_thread = (double *)malloc((*n) * sizeof(double));
	//	start_index = f_thread * (*n);
	//	for (i_thread = 0; i_thread < (*n); i_thread++)
	//	{
	//		x_thread[i_thread] = x[start_index + i_thread];
	//	}

	//	U_Ball_Correlation(bcorsis_stat_tmp, n, x_thread, yrank, lowyidx, higyidx);
	//	for (s_thread = 0; s_thread < 3; s_thread++)
	//	{
	//		bcorsis_stat[3 * f_thread + s_thread] = bcorsis_stat_tmp[s_thread];
	//	}
	//	free(x_thread);
	//}

	free(yrank);
	free_int_matrix(lowyidx, *n, *n);
	free_int_matrix(higyidx, *n, *n);
	return;
}


void _bcor_stat(double *bcorsis_stat, double *y, double *x, int *n)
{
	double **Dx, **Dy, *x_cpy, *y_cpy;
	int **xidx, **yidx;
	Dx = alloc_matrix(*n, *n);
	Dy = alloc_matrix(*n, *n);
	xidx = alloc_int_matrix(*n, *n);
	yidx = alloc_int_matrix(*n, *n);

	x_cpy = (double *)malloc((*n) * sizeof(double));
	y_cpy = (double *)malloc((*n) * sizeof(double));
	vector2matrix(x, Dx, *n, *n, 1);
	vector2matrix(y, Dy, *n, *n, 1);

	int s, t;
	for (s = 0; s < *n; s++)
	{
		for (t = 0; t<*n; t++)
		{
			xidx[s][t] = t;
			yidx[s][t] = t;
		}
	}

	for (s = 0; s<(*n); s++)
	{
		// copy site to x_cpy and y_cpy
		memcpy(x_cpy, Dx[s], *n * sizeof(double));
		memcpy(y_cpy, Dy[s], *n * sizeof(double));
		quicksort(x_cpy, xidx[s], 0, *n - 1);
		quicksort(y_cpy, yidx[s], 0, *n - 1);
	}
	free(x_cpy);
	free(y_cpy);

	int i, j, k, pi, src, lastpos, *yrank, *isource, *icount, *xy_index, *xy_temp, **xyidx;
	double pxy, px, py, lastval, *xx_cpy, *yy_cpy;
	double bcov_weight0 = 0.0, bcov_weight_prob = 0.0, bcov_weight_hhg = 0.0, bcov_fixed_ball = 0.0;
	double bcov_weight0_x = 0.0, bcov_weight_prob_x = 0.0;
	double bcov_weight0_y = 0.0, bcov_weight_prob_y = 0.0;
	double hhg_ball_num = 0.0, minor_ball_prop = 2 / (*n);

	yrank = (int *)malloc(*n * sizeof(int));
	isource = (int *)malloc(*n * sizeof(int));
	icount = (int *)malloc(*n * sizeof(int));
	xy_index = (int *)malloc(*n * sizeof(int));
	xy_temp = (int *)malloc(*n * sizeof(int));
	xyidx = alloc_int_matrix(*n, *n);
	xx_cpy = (double *)malloc(*n * sizeof(double));
	yy_cpy = (double *)malloc(*n * sizeof(double));

	for (i = 0; i<*n; i++)
		for (j = 0; j<*n; j++)
			xyidx[i][j] = j;

	for (i = 0; i<(*n); i++)
	{
		memcpy(xx_cpy, Dx[i], *n * sizeof(double));
		for (j = 0; j<*n; j++)
			yy_cpy[j] = Dy[i][j];
		quicksort2(xx_cpy, yy_cpy, xyidx[i], 0, *n - 1);
	}
	free(xx_cpy);
	free(yy_cpy);

	for (i = 0; i<(*n); i++)
	{
		pi = i;
		lastval = 0;
		lastpos = -1;
		for (j = *n - 1, k = *n - 1; j >= 1; --j, --k)
		{
			k -= (yidx[pi][k] == pi);
			if (lastpos == -1 || Dy[pi][yidx[pi][k]] != lastval)
			{
				lastval = Dy[pi][yidx[pi][k]];
				lastpos = j;
			}
			src = yidx[pi][k];
			src -= (src > i);
			yrank[src] = lastpos;
		}

		for (j = 0, k = 0; j < *n - 1; ++j, ++k)
		{
			k += (xyidx[i][k] == i);
			src = xyidx[i][k]; // NOTE: k may be different than from the line above
			src -= (src > i);
			xy_index[j] = yrank[src];
			isource[j] = j;
			icount[j] = 0;
			xy_temp[j] = xy_index[j];
		}
		Inversions(xy_temp, isource, icount, *n - 1, *n);
		lastval = 0;
		lastpos = -1;

		for (j = *n - 2, k = *n - 1; j >= 0; --j, --k)
		{
			k -= (xidx[i][k] == i);
			if (lastpos == -1 || Dx[i][xidx[i][k]] != lastval) {
				lastval = Dx[i][xidx[i][k]];
				lastpos = j;
			}

			pxy = lastpos - icount[j] + 2;
			px = lastpos + 2;
			py = xy_index[j] + 1;
			px /= (*n);
			py /= (*n);
			pxy /= (*n);
			// compute BCov(X, Y)
			bcov_fixed_ball = pow(pxy - px*py, 2);
			bcov_weight0 += bcov_fixed_ball;
			bcov_weight_prob += bcov_fixed_ball / (px*py);
			if (px > minor_ball_prop && py > minor_ball_prop && px != 1 && py != 1)
			{
				bcov_weight_hhg += bcov_fixed_ball / ((px - minor_ball_prop)*(1.0 - px + minor_ball_prop)*(py - minor_ball_prop)*(1.0 - py + minor_ball_prop));
				hhg_ball_num += 1;
			}
			// compute BCov(X, X)
			bcov_weight_prob_x += (1.0 - px)*(1.0 - px);
			bcov_weight0_x += px*px*(1.0 - px)*(1.0 - px);

			// compute BCov(Y, Y)
			bcov_weight_prob_y += (1.0 - py)*(1.0 - py);
			bcov_weight0_y += py*py*(1.0 - py)*(1.0 - py);
		}
		pxy = 0;
		px = 0;
		py = 0;
		for (j = 0; j<*n; j++)
		{
			if (Dx[i][xidx[i][j]] == 0)
			{
				px += 1;
				if (Dy[pi][xidx[i][j]] == 0)
				{
					pxy += 1;
					py += 1;
				}
			}
			else if (Dy[pi][xidx[i][j]] == 0)
				py += 1;
		}
		px /= (*n);
		py /= (*n);
		pxy /= (*n);
		// compute BCov(X, Y)
		bcov_fixed_ball = pow(pxy - px*py, 2);
		bcov_weight0 += bcov_fixed_ball;
		bcov_weight_prob += bcov_fixed_ball / (px*py);
		if (px > minor_ball_prop && py > minor_ball_prop && px != 1 && py != 1)
		{
			bcov_weight_hhg += bcov_fixed_ball / ((px - minor_ball_prop)*(1.0 - px + minor_ball_prop)*(py - minor_ball_prop)*(1.0 - py + minor_ball_prop));
			hhg_ball_num += 1;
		}
		// compute BCov(X, X)
		bcov_weight_prob_x += (1.0 - px)*(1.0 - px);
		bcov_weight0_x += px*px*(1.0 - px)*(1.0 - px);

		// compute BCov(Y, Y)
		bcov_weight_prob_y += (1.0 - py)*(1.0 - py);
		bcov_weight0_y += py*py*(1.0 - py)*(1.0 - py);
	}
	bcorsis_stat[0] = bcov_weight0 / (sqrt(bcov_weight0_x) * sqrt(bcov_weight0_y));
	bcorsis_stat[1] = bcov_weight_prob / (sqrt(bcov_weight_prob_x) * sqrt(bcov_weight_prob_y));
	bcorsis_stat[2] = bcov_weight_hhg / (hhg_ball_num);

	// free memory
	free_int_matrix(xidx, *n, *n);
	free_int_matrix(yidx, *n, *n);
	free_matrix(Dx, *n, *n);
	free_matrix(Dy, *n, *n);
	free(isource);
	free(icount);
	free(xy_index);
	free(yrank);
	free(xy_temp);
	free_int_matrix(xyidx, *n, *n);
	return;
}


// R API function:
/*
bcorsis_stat: bcor statistics or p-value for screening
y: response
x: covariate
x_number: if x_number = [1, 2, 2, 1, 1], then x[:, 1], x[:, 2:3], x[:, 4:5], x[:, 6], x[:, 7] will be used to compute ball correlation
f_number: the number of covariate
size: sample size of each group
n: total sample size
p: dimensionlity of response variable
k: group number
dst_y: whether y should be recompute as distance
dst_x: whether x should be recompute as distance
R: permutation replication
nthread: control the number threads used to compute statistics
*/
void bcor_test(double *bcorsis_stat, double *y, double *x, int *x_number, int *f_number, int *size, int *n, int *p, int *k, int *dst_y, int *dst_x, int *nthread)
{
	// 
	if (*dst_y == 1) {
		if (*dst_x == 0) {
			_bcor_test(bcorsis_stat, y, x, x_number, f_number, n, p, nthread);
		} else {
			_bcor_stat(bcorsis_stat, y, x, n);
		}
	}
	else {
		// If y is univariate, and y indicator for K classes. We can simplify the computation for ball correlation by 
		// bcor_two_sample() and bcor_k_sample() function;
		// x may be univariate, multivariate or their combination
		if (*k > 1) {
			// TODO:
			if (*k == 2) {
				//bcor_two_sample(bcorsis_stat, x, x_number, size, n, nthread);
			}
			// TODO:
			else {
				//bcor_k_sample(bcorsis_stat, x, x_number, size, n, k, nthread);
			}
		}
		else {
			_fast_bcor_test(bcorsis_stat, y, x, f_number, n, nthread);
		}
	}
	return;
}



//double bcov_value(double *x, double *y, int *n, int *weight, int *thread)
//{
//	/*  computes RCT(x,y)  */
//	int    i, j, **xidx, **yidx, *i_perm, *i_perm_inv;
//	double **Dx, **Dy, *x_cpy, *y_cpy, RCTV0;
//	double bcov[3];
//
//	Dx = alloc_matrix(*n, *n);
//	Dy = alloc_matrix(*n, *n);
//	xidx = alloc_int_matrix(*n, *n);
//	yidx = alloc_int_matrix(*n, *n);
//	i_perm = (int *)malloc(*n * sizeof(int));
//	i_perm_inv = (int *)malloc(*n * sizeof(int));
//	x_cpy = (double *)malloc(*n * sizeof(double));
//	y_cpy = (double *)malloc(*n * sizeof(double));
//
//	vector2matrix(x, Dx, *n, *n, 1);
//	vector2matrix(y, Dy, *n, *n, 1);
//
//	for (i = 0; i<*n; i++)
//	{
//		for (j = 0; j<*n; j++)
//		{
//			xidx[i][j] = j;
//			yidx[i][j] = j;
//		}
//		i_perm[i] = i;
//		i_perm_inv[i] = i;
//	}
//
//	for (i = 0; i<(*n); i++)
//	{
//		// copy site to x_cpy and y_cpy
//		memcpy(x_cpy, Dx[i], *n * sizeof(double));
//		memcpy(y_cpy, Dy[i], *n * sizeof(double));
//		quicksort(x_cpy, xidx[i], 0, *n - 1);
//		quicksort(y_cpy, yidx[i], 0, *n - 1);
//	}
//	free(x_cpy);
//	free(y_cpy);
//
//	Ball_Information_wrapper(bcov, n, Dx, Dy, xidx, yidx, i_perm, i_perm_inv, thread);
//
//	if (*weight == 0)
//	{
//		RCTV0 = bcov[0];
//	}
//	else {
//		RCTV0 = bcov[1];
//	}
//
//	free_matrix(Dx, *n, *n);
//	free_matrix(Dy, *n, *n);
//	free_int_matrix(xidx, *n, *n);
//	free_int_matrix(yidx, *n, *n);
//	free(i_perm);
//	free(i_perm_inv);
//	return(RCTV0);
//}
//
//
//double bcor_value(double *x, double *y, int *n, int *weight, int *dst, int *thread)
//{
//	double bcov_stat_xy, bcov_stat_xx, bcov_stat_yy;
//	double bcor_stat;
//	if ((*dst)) {
//		bcov_stat_xy = bcov_value(x, y, n, weight, thread);
//		bcov_stat_xx = bcov_value(x, x, n, weight, thread);
//		bcov_stat_yy = bcov_value(x, x, n, weight, thread);
//	}
//	else {
//		bcov_stat_xy = ubcov_value(x, y, n, weight, thread);
//		bcov_stat_xx = ubcov_value(x, x, n, weight, thread);
//		bcov_stat_yy = ubcov_value(x, y, n, weight, thread);
//	}
//	bcor_stat = bcov_stat_xy / sqrt(bcov_stat_xx*bcov_stat_yy);
//	return(bcor_stat);
//}
//
//
//void UBCOR(double *bcor, double *permuted_bcor, double *x, double *y, int *n, int *R, int *weight, int *thread)
//{
//	int i;
//	double bcov_stat_xx, bcov_stat_yy;
//	UBI(bcor, permuted_bcor, x, y, n, R, thread);
//	bcov_stat_xx = ubcov_value(x, x, n, weight, thread);
//	bcov_stat_yy = ubcov_value(y, y, n, weight, thread);
//	(*bcor) = (*bcor) / bcov_stat_xx / bcov_stat_yy;
//	for (i = 0; i < (*R); i++) {
//		permuted_bcor[i] = permuted_bcor[i] / bcov_stat_xx / bcov_stat_yy;
//	}
//}
//
//
//void BCOR(double *bcor, double *permuted_bcor, double *x, double *y, int *n, int *R, int *weight, int *thread)
//{
//	int i;
//	double bcov_stat_xx, bcov_stat_yy;
//	BI(bcor, permuted_bcor, x, y, n, R, thread);
//	bcov_stat_xx = bcov_value(x, x, n, weight, thread);
//	bcov_stat_yy = bcov_value(y, y, n, weight, thread);
//	(*bcor) = (*bcor) / bcov_stat_xx / bcov_stat_yy;
//	for (i = 0; i < (*R); i++) {
//		permuted_bcor[i] = permuted_bcor[i] / bcov_stat_xx / bcov_stat_yy;
//	}
//}
//
//
//void bcov_stat(double *bcov, double *x, double *y, int *n, int *weight, int *dst, int *type, int *thread)
//{
//	double bcov_stat;
//	if ((*type) == 1) {
//		if ((*dst)) {
//			bcov_stat = bcov_value(x, y, n, weight, thread);
//		}
//		else {
//			bcov_stat = ubcov_value(x, y, n, weight, thread);
//		}
//		*bcov = bcov_stat;
//	}
//	else {
//		*bcov = bcor_value(x, y, n, weight, dst, thread);
//	}
//	return;
//}