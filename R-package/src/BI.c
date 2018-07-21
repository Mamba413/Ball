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
#include "stdlib.h" 
#include "string.h"
#include "stdio.h"
#include "utilities.h"
#include "BI.h"


void Ball_Information(double *bcov_stat, int *n, double **Dx, double **Dy, int **xidx, int **yidx, int *i_perm, int *i_perm_inv)
{
	int i, j, k, pi, src, lastpos, *yrank, *isource, *icount, *xy_index, *xy_temp, **xyidx;
	double pxy, px, py, lastval, *xx_cpy, *yy_cpy;
	double bcov_weight0 = 0.0, bcov_weight_prob = 0.0, bcov_weight_hhg = 0.0, bcov_fixed_ball = 0.0;
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
			yy_cpy[j] = Dy[i_perm[i]][i_perm[j]];
		quicksort2(xx_cpy, yy_cpy, xyidx[i], 0, *n - 1);
	}
	free(xx_cpy);
	free(yy_cpy);

	for (i = 0; i<(*n); i++)
	{
		pi = i_perm[i];
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
			src = i_perm_inv[yidx[pi][k]];
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
			bcov_fixed_ball = pow(pxy - px*py, 2);
			bcov_weight0 += bcov_fixed_ball;
			bcov_weight_prob += bcov_fixed_ball / (px*py);
			if (px > minor_ball_prop && py > minor_ball_prop && px != 1 && py != 1)
			{
				bcov_weight_hhg += bcov_fixed_ball / ((px - minor_ball_prop)*(1.0 - px + minor_ball_prop)*(py - minor_ball_prop)*(1.0 - py + minor_ball_prop));
				hhg_ball_num += 1.0;
			}
		}
		pxy = 0;
		px = 0;
		py = 0;
		for (j = 0; j<*n; j++)
		{
			if (Dx[i][xidx[i][j]] == 0)
			{
				px += 1;
				if (Dy[pi][i_perm[xidx[i][j]]] == 0)
				{
					pxy += 1;
					py += 1;
				}
			}
			else if (Dy[pi][i_perm[xidx[i][j]]] == 0)
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
			hhg_ball_num += 1.0;
		}
	}
	bcov_stat[0] = bcov_weight0 / (1.0*(*n)*(*n));
	bcov_stat[1] = bcov_weight_prob / (1.0*(*n)*(*n));
	bcov_stat[2] = bcov_weight_hhg / (hhg_ball_num);

	// free memory
	free(isource);
	free(icount);
	free(xy_index);
	free(yrank);
	free(xy_temp);
	free_int_matrix(xyidx, *n, *n);
	return;
}


void Ball_Information_parallel(double *bcov_stat, int *n, double **Dx, double **Dy, int **xidx, int **yidx, int *i_perm, int *i_perm_inv, int *nthread)
{
	double bcov_weight0 = 0.0, bcov_weight_prob = 0.0, bcov_weight_hhg = 0.0;
	double minor_ball_prop = 2 / (*n);
	int i, jjjj, **xyidx;
	double *xx_cpy, *yy_cpy;
	xyidx = alloc_int_matrix(*n, *n);
	for (i = 0; i<*n; i++)
		for (jjjj = 0; jjjj<*n; jjjj++)
			xyidx[i][jjjj] = jjjj;

	xx_cpy = (double *)malloc(*n * sizeof(double));
	yy_cpy = (double *)malloc(*n * sizeof(double));

	for (i = 0; i<(*n); i++) {
		memcpy(xx_cpy, Dx[i], *n * sizeof(double));
		for (jjjj = 0; jjjj<*n; jjjj++)
			yy_cpy[jjjj] = Dy[i_perm[i]][i_perm[jjjj]];
		quicksort2(xx_cpy, yy_cpy, xyidx[i], 0, *n - 1);
	}
	free(xx_cpy);
	free(yy_cpy);

#ifdef _OPENMP
	omp_set_num_threads(*nthread);
#endif
#pragma omp parallel
	{
		int j, k, pi, src, lastpos, *yrank, *isource, *icount, *xy_index, *xy_temp;
		double pxy, px, py, lastval;
		double bcov_weight0_thread = 0.0, bcov_weight_prob_thread = 0.0, bcov_weight_hhg_thread = 0.0, bcov_fixed_ball = 0.0;
		yrank = (int *)malloc(*n * sizeof(int));
		isource = (int *)malloc(*n * sizeof(int));
		icount = (int *)malloc(*n * sizeof(int));
		xy_index = (int *)malloc(*n * sizeof(int));
		xy_temp = (int *)malloc(*n * sizeof(int));

#pragma omp for
		for (i = 0; i<(*n); i++)
		{
			pi = i_perm[i];
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
				src = i_perm_inv[yidx[pi][k]];
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
				bcov_fixed_ball = pow(pxy - px*py, 2);
				bcov_weight0_thread += bcov_fixed_ball;
				bcov_weight_prob_thread += bcov_fixed_ball / (px*py);
				if (px > minor_ball_prop && py > minor_ball_prop && px != 1 && py != 1)
				{
					bcov_weight_hhg_thread += bcov_fixed_ball / ((px - minor_ball_prop)*(1.0 - px + minor_ball_prop)*(py - minor_ball_prop)*(1.0 - py + minor_ball_prop));
				}
			}

			pxy = 0;
			px = 0;
			py = 0;

			for (j = 0; j<*n; j++)
			{
				if (Dx[i][xidx[i][j]] == 0)
				{
					px += 1;
					if (Dy[pi][i_perm[xidx[i][j]]] == 0)
					{
						pxy += 1;
						py += 1;
					}
				}
				else if (Dy[pi][i_perm[xidx[i][j]]] == 0)
					py += 1;
			}

			px /= (*n);
			py /= (*n);
			pxy /= (*n);
			bcov_fixed_ball = (pxy - px*py) * (pxy - px*py);
			bcov_weight0_thread += bcov_fixed_ball;
			bcov_weight_prob_thread += bcov_fixed_ball / (px*py);
			if (px > minor_ball_prop && py > minor_ball_prop && px != 1 && py != 1)
			{
				bcov_weight_hhg_thread += bcov_fixed_ball / ((px - minor_ball_prop)*(1.0 - px + minor_ball_prop)*(py - minor_ball_prop)*(1.0 - py + minor_ball_prop));
			}
		}
#pragma omp critical
		{
			bcov_weight0 += bcov_weight0_thread;
			bcov_weight_prob += bcov_weight_prob_thread;
			bcov_weight_hhg += bcov_weight_hhg_thread;
		}
		bcov_stat[0] = bcov_weight0 / (1.0*(*n)*(*n));
		bcov_stat[1] = bcov_weight_prob / (1.0*(*n)*(*n));
		bcov_stat[2] = bcov_weight_hhg / (1.0*(*n)*(*n));

		free(isource);
		free(icount);
		free(xy_index);
		free(yrank);
		free(xy_temp);
	}
	free_int_matrix(xyidx, *n, *n);
	return;
}


void Ball_Information_wrapper(double *bcov_stat, int *n, double **Dx, double **Dy, int **xidx, int **yidx, int *i_perm, int *i_perm_inv, int *nthread)
{
	if (*nthread == 1) {
		Ball_Information(bcov_stat, n, Dx, Dy, xidx, yidx, i_perm, i_perm_inv);
	}
	else {
		Ball_Information_parallel(bcov_stat, n, Dx, Dy, xidx, yidx, i_perm, i_perm_inv, nthread);
	}
	return;
}


void BI(double *bcov, double *pvalue, double *x, double *y, int *n, int *R, int *thread)
{
	/*  computes RCT(x,y)  */
	int    i, j, **xidx, **yidx, *i_perm, *i_perm_inv;
	double **Dx, **Dy, *x_cpy, *y_cpy;

	Dx = alloc_matrix(*n, *n);
	Dy = alloc_matrix(*n, *n);
	xidx = alloc_int_matrix(*n, *n);
	yidx = alloc_int_matrix(*n, *n);
	i_perm = (int *)malloc(*n * sizeof(int));
	i_perm_inv = (int *)malloc(*n * sizeof(int));
	x_cpy = (double *)malloc(*n * sizeof(double));
	y_cpy = (double *)malloc(*n * sizeof(double));

	// convert distance vector to distance matrix:
	vector2matrix(x, Dx, *n, *n, 1);
	vector2matrix(y, Dy, *n, *n, 1);

	for (i = 0; i<*n; i++)
	{
		for (j = 0; j<*n; j++)
		{
			xidx[i][j] = j;
			yidx[i][j] = j;
		}
		i_perm[i] = i;
		i_perm_inv[i] = i;
	}

	// sort each row of Dx and Dy
	// computational complexity: O(n^2 * logn)
	// xidx, yidx is index of each row after sorted
	for (i = 0; i<(*n); i++)
	{
		// copy site to x_cpy and y_cpy
		memcpy(x_cpy, Dx[i], *n * sizeof(double));
		memcpy(y_cpy, Dy[i], *n * sizeof(double));
		quicksort(x_cpy, xidx[i], 0, *n - 1);
		quicksort(y_cpy, yidx[i], 0, *n - 1);
	}
	free(x_cpy);
	free(y_cpy);

	Ball_Information_wrapper(bcov, n, Dx, Dy, xidx, yidx, i_perm, i_perm_inv, thread);
	// printf("RCTV0 = %f\n", bcov[0]);
	if (*R > 0)
	{
		double bcov_tmp[3], *permuted_bcov_weight0, *permuted_bcov_weight_prob, *permuted_bcov_weight_hhg;
		permuted_bcov_weight0 = (double *)malloc(*R * sizeof(double));
		permuted_bcov_weight_prob = (double *)malloc(*R * sizeof(double));
		permuted_bcov_weight_hhg = (double *)malloc(*R * sizeof(double));
    
    // int iii;
		for (i = 0; i<*R; i++)
		{
			// stop permutation if user stop it manually:
			if (pending_interrupt()) {
				print_stop_message();
				break;
			}
			resample(i_perm, i_perm_inv, n);
			/*
			for (iii = 0; iii < *n; iii++)
			{
			  printf("%d ", i_perm[iii]);
			}
			printf("\n");
			*/
			Ball_Information_wrapper(bcov_tmp, n, Dx, Dy, xidx, yidx, i_perm, i_perm_inv, thread);
			permuted_bcov_weight0[i] = bcov_tmp[0]; permuted_bcov_weight_prob[i] = bcov_tmp[1]; permuted_bcov_weight_hhg[i] = bcov_tmp[2];
	    // printf("i = %d, RCTV1 = %f\n", i, bcov_tmp[0]);
		}
		pvalue[0] = compute_pvalue(bcov[0], permuted_bcov_weight0, i);
		pvalue[1] = compute_pvalue(bcov[1], permuted_bcov_weight_prob, i);
		pvalue[2] = compute_pvalue(bcov[2], permuted_bcov_weight_hhg, i);

		free(permuted_bcov_weight0); free(permuted_bcov_weight_prob); free(permuted_bcov_weight_hhg);
	}

	free_matrix(Dx, *n, *n);
	free_matrix(Dy, *n, *n);
	free_int_matrix(xidx, *n, *n);
	free_int_matrix(yidx, *n, *n);
	free(i_perm);
	free(i_perm_inv);
	return;
}


void BI_parallel(double *bcov, double *pvalue, double *x, double *y, int *n, int *R, int *thread)
{
	int    i, j, **xidx, **yidx, *i_perm, *i_perm_inv;
	double **Dx, **Dy, *x_cpy, *y_cpy;

	Dx = alloc_matrix(*n, *n);
	Dy = alloc_matrix(*n, *n);
	xidx = alloc_int_matrix(*n, *n);
	yidx = alloc_int_matrix(*n, *n);
	i_perm = (int *)malloc(*n * sizeof(int));
	i_perm_inv = (int *)malloc(*n * sizeof(int));
	x_cpy = (double *)malloc(*n * sizeof(double));
	y_cpy = (double *)malloc(*n * sizeof(double));

	vector2matrix(x, Dx, *n, *n, 1);
	vector2matrix(y, Dy, *n, *n, 1);

	for (i = 0; i<*n; i++)
	{
		for (j = 0; j<*n; j++)
		{
			xidx[i][j] = j;
			yidx[i][j] = j;
		}
		i_perm[i] = i;
		i_perm_inv[i] = i;
	}

	for (i = 0; i<(*n); i++)
	{
		// copy site to x_cpy and y_cpy
		memcpy(x_cpy, Dx[i], *n * sizeof(double));
		memcpy(y_cpy, Dy[i], *n * sizeof(double));
		quicksort(x_cpy, xidx[i], 0, *n - 1);
		quicksort(y_cpy, yidx[i], 0, *n - 1);
	}
	free(x_cpy);
	free(y_cpy);

	Ball_Information(bcov, n, Dx, Dy, xidx, yidx, i_perm, i_perm_inv);
	// printf("RCTV0 = %f\n", bcov[0]);
	/*
	free_int_matrix(xidx, *n, *n);
	free_int_matrix(yidx, *n, *n);
	 */
	free(i_perm);
	free(i_perm_inv);

	if (*R > 0) {
		double *permuted_bcov_weight0, *permuted_bcov_weight_prob, *permuted_bcov_weight_hhg;
		permuted_bcov_weight0 = (double *)malloc(*R * sizeof(double));
		permuted_bcov_weight_prob = (double *)malloc(*R * sizeof(double));
		permuted_bcov_weight_hhg = (double *)malloc(*R * sizeof(double));

#ifdef _OPENMP
		omp_set_num_threads(*thread);
#endif
#pragma omp parallel
		{
			// int kk, **xidx_thread, **yidx_thread;
			// xidx_thread = alloc_int_matrix(*n, *n);
			// yidx_thread = alloc_int_matrix(*n, *n);
			int *i_perm_thread, *i_perm_inv_thread, k, i_thread;
			double bcov_tmp[3];
			i_perm_thread = (int *)malloc(*n * sizeof(int));
			i_perm_inv_thread = (int *)malloc(*n * sizeof(int));
#pragma omp critical
			{
				for (k = 0; k < *n; k++)
				{
				  /*
					for (kk = 0; kk < *n; kk++)
					{
						xidx_thread[k][kk] = kk;
						yidx_thread[k][kk] = kk;
					}
				   */
					i_perm_thread[k] = k;
					i_perm_inv_thread[k] = k;
				}
			}

#pragma omp for
			for (i_thread = 0; i_thread < (*R); i_thread++)
			{
#pragma omp critical
				{
					resample(i_perm_thread, i_perm_inv_thread, n);
				}
				Ball_Information(bcov_tmp, n, Dx, Dy, xidx, yidx, i_perm_thread, i_perm_inv_thread);
				permuted_bcov_weight0[i_thread] = bcov_tmp[0]; 
				permuted_bcov_weight_prob[i_thread] = bcov_tmp[1]; 
				permuted_bcov_weight_hhg[i_thread] = bcov_tmp[2];
				// printf("RCTV1 = %f\n", bcov_tmp[0]);
			}
			free(i_perm_thread);
			free(i_perm_inv_thread);
			/*
			free_int_matrix(xidx_thread, *n, *n);
			free_int_matrix(yidx_thread, *n, *n);
			 */
		}
		pvalue[0] = compute_pvalue(bcov[0], permuted_bcov_weight0, *R);
		pvalue[1] = compute_pvalue(bcov[1], permuted_bcov_weight_prob, *R);
		pvalue[2] = compute_pvalue(bcov[2], permuted_bcov_weight_hhg, *R);
		free(permuted_bcov_weight0); free(permuted_bcov_weight_prob); free(permuted_bcov_weight_hhg);
	}

	free_matrix(Dx, *n, *n);
	free_matrix(Dy, *n, *n);
	return;
}


void U_Ball_Information(double *bcov_stat, int *n, int **Rank, int **lowxidx, int **higxidx, int **lowyidx, int **higyidx, int *i_perm)
{
	int i, j, pi, pj;
	double px, py, pxy;
	double bcov_weight0 = 0.0, bcov_weight_prob = 0.0, bcov_weight_hhg = 0.0, bcov_fixed_ball = 0.0;
	double hhg_ball_num = 0.0, minor_ball_prop = 2 / (*n);
	for (i = 0; i < *n; i++) {
		for (j = 0; j < *n; j++) {
			pi = i_perm[i];
			pj = i_perm[j];
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
		}
	}
	bcov_stat[0] = bcov_weight0 / (1.0*(*n)*(*n));
	bcov_stat[1] = bcov_weight_prob / (1.0*(*n)*(*n));
	bcov_stat[2] = bcov_weight_hhg / hhg_ball_num;

	return;
}


void U_Ball_Information_parallel(double *bcov_stat, int *n, int **Rank, int **lowxidx, int **higxidx, int **lowyidx, int **higyidx, int *i_perm, int *nthread)
{
	double bcov_weight0 = 0.0, bcov_weight_prob = 0.0, bcov_weight_hhg = 0.0;
	double minor_ball_prop = 2 / (*n);
#ifdef _OPENMP
	omp_set_num_threads(*nthread);
#endif
#pragma omp parallel
	{
		int i, j, pi, pj;
		double px, py, pxy;
		double bcov_weight0_thread = 0.0, bcov_weight_prob_thread = 0.0, bcov_weight_hhg_thread = 0.0, bcov_fixed_ball = 0.0;

#pragma omp for
		for (i = 0; i < *n; i++) {
			for (j = 0; j < *n; j++) {
				pi = i_perm[i];
				pj = i_perm[j];
				px = higxidx[i][j] - lowxidx[i][j] + 1;
				py = higyidx[pi][pj] - lowyidx[pi][pj] + 1;
				pxy = Rank[higxidx[i][j]][higyidx[pi][pj]] + Rank[lowxidx[i][j] - 1][lowyidx[pi][pj] - 1];
				pxy -= (Rank[higxidx[i][j]][lowyidx[pi][pj] - 1] + Rank[lowxidx[i][j] - 1][higyidx[pi][pj]]);

				pxy /= (*n);
				px /= (*n);
				py /= (*n);

				bcov_fixed_ball = (pxy - px*py) * (pxy - px*py);
				bcov_weight0_thread += bcov_fixed_ball;
				bcov_weight_prob_thread += bcov_fixed_ball / (px*py);
				if (px > minor_ball_prop && py > minor_ball_prop && px != 1 && py != 1)
				{
					bcov_weight_hhg_thread += bcov_fixed_ball / ((px - minor_ball_prop)*(1.0 - px + minor_ball_prop)*(py - minor_ball_prop)*(1.0 - py + minor_ball_prop));
				}
			}
		}

#pragma omp critical
		{
			bcov_weight0 += bcov_weight0_thread;
			bcov_weight_prob += bcov_weight_prob_thread;
			bcov_weight_hhg += bcov_weight_hhg_thread;
		}
	}
	bcov_stat[0] = bcov_weight0;
	bcov_stat[1] = bcov_weight_prob;
	bcov_stat[2] = bcov_weight_hhg;
	return;
}


void U_Ball_Information_wrapper(double *bcov_stat, int *n, int **Rank, int **lowxidx, int **higxidx, int **lowyidx, int **higyidx, int *i_perm, int *nthread)
{
	if ((*nthread) == 1)
	{
		U_Ball_Information(bcov_stat, n, Rank, lowxidx, higxidx, lowyidx, higyidx, i_perm);
	}
	else {
		U_Ball_Information_parallel(bcov_stat, n, Rank, lowxidx, higxidx, lowyidx, higyidx, i_perm, nthread);
	}
	return;
}


void UBI(double *bcov, double *pvalue, double *x, double *y, int *n, int *R, int *thread)
{
	int i, j, *xidx, *yidx, *xrank, *yrank, *i_perm, **Rank, **lowxidx, **higxidx, **lowyidx, **higyidx;

	xidx = (int *)malloc(*n * sizeof(int));
	yidx = (int *)malloc(*n * sizeof(int));
	xrank = (int *)malloc(*n * sizeof(int));
	yrank = (int *)malloc(*n * sizeof(int));
	i_perm = (int *)malloc(*n * sizeof(int));
	Rank = alloc_int_matrix((*n) + 1, (*n) + 1);
	lowxidx = alloc_int_matrix(*n, *n);
	higxidx = alloc_int_matrix(*n, *n);
	lowyidx = alloc_int_matrix(*n, *n);
	higyidx = alloc_int_matrix(*n, *n);

	for (i = 0; i<*n; i++)
	{
		xidx[i] = i;
		yidx[i] = i;
		i_perm[i] = i;
	}

	// first step: sort the x, y and obtain index of x, y (after the x, y have been sorted)
	quicksort(x, xidx, 0, *n - 1);
	quicksort(y, yidx, 0, *n - 1);
	ranksort(n, xrank, x, xidx);
	ranksort(n, yrank, y, yidx);
	createidx(n, xidx, x, lowxidx, higxidx);
	createidx(n, yidx, y, lowyidx, higyidx);


	initRank(*n, Rank, xrank, yrank, i_perm);
	U_Ball_Information_wrapper(bcov, n, Rank, lowxidx, higxidx, lowyidx, higyidx, i_perm, thread);
	if (*R > 0)
	{
		double bcov_tmp[3], *permuted_bcov_weight0, *permuted_bcov_weight_prob, *permuted_bcov_weight_hhg;
		permuted_bcov_weight0 = (double *)malloc(*R * sizeof(double));
		permuted_bcov_weight_prob = (double *)malloc(*R * sizeof(double));
		permuted_bcov_weight_hhg = (double *)malloc(*R * sizeof(double));

		for (j = 0; j<*R; j++)
		{
			// stop permutation if user stop it manually:
			if (pending_interrupt()) {
				print_stop_message();
				break;
			}
			resample2(i_perm, n);
			initRank(*n, Rank, xrank, yrank, i_perm);
			U_Ball_Information_wrapper(bcov_tmp, n, Rank, lowxidx, higxidx, lowyidx, higyidx, i_perm, thread);
			permuted_bcov_weight0[j] = bcov_tmp[0]; permuted_bcov_weight_prob[j] = bcov_tmp[1]; permuted_bcov_weight_hhg[j] = bcov_tmp[2];
		}
		pvalue[0] = compute_pvalue(bcov[0], permuted_bcov_weight0, j);
		pvalue[1] = compute_pvalue(bcov[1], permuted_bcov_weight_prob, j);
		pvalue[2] = compute_pvalue(bcov[2], permuted_bcov_weight_hhg, j);

		free(permuted_bcov_weight0); free(permuted_bcov_weight_prob); free(permuted_bcov_weight_hhg);
	}

	free_int_matrix(Rank, (*n) + 1, (*n) + 1);
	free_int_matrix(lowxidx, *n, *n);
	free_int_matrix(higxidx, *n, *n);
	free_int_matrix(lowyidx, *n, *n);
	free_int_matrix(higyidx, *n, *n);
	free(xidx);
	free(yidx);
	free(xrank);
	free(yrank);
	free(i_perm);
	return;
}


void UBI_parallel(double *bcov, double *pvalue, double *x, double *y, int *n, int *R, int *thread)
{
	int i, *xidx, *yidx, *xrank, *yrank, *i_perm, **Rank, **lowxidx, **higxidx, **lowyidx, **higyidx;

	xidx = (int *)malloc(*n * sizeof(int));
	yidx = (int *)malloc(*n * sizeof(int));
	xrank = (int *)malloc(*n * sizeof(int));
	yrank = (int *)malloc(*n * sizeof(int));
	i_perm = (int *)malloc(*n * sizeof(int));
	Rank = alloc_int_matrix((*n) + 1, (*n) + 1);
	lowxidx = alloc_int_matrix(*n, *n);
	higxidx = alloc_int_matrix(*n, *n);
	lowyidx = alloc_int_matrix(*n, *n);
	higyidx = alloc_int_matrix(*n, *n);

	for (i = 0; i<*n; i++)
	{
		xidx[i] = i;
		yidx[i] = i;
		i_perm[i] = i;
	}

	// first step: sort the x
	quicksort(x, xidx, 0, *n - 1);
	quicksort(y, yidx, 0, *n - 1);
	ranksort(n, xrank, x, xidx);
	ranksort(n, yrank, y, yidx);
	createidx(n, xidx, x, lowxidx, higxidx);
	createidx(n, yidx, y, lowyidx, higyidx);

	initRank(*n, Rank, xrank, yrank, i_perm);
	U_Ball_Information(bcov, n, Rank, lowxidx, higxidx, lowyidx, higyidx, i_perm);

	free_int_matrix(Rank, (*n) + 1, (*n) + 1);
	free(i_perm);
	free(xidx);
	free(yidx);

	if (*R > 0)
	{
		double *permuted_bcov_weight0, *permuted_bcov_weight_prob, *permuted_bcov_weight_hhg;
		permuted_bcov_weight0 = (double *)malloc(*R * sizeof(double));
		permuted_bcov_weight_prob = (double *)malloc(*R * sizeof(double));
		permuted_bcov_weight_hhg = (double *)malloc(*R * sizeof(double));

#ifdef _OPENMP
		omp_set_num_threads(*thread);
#endif
#pragma omp parallel
		{
			int **Rank_thread, *i_perm_thread;
			int k, j_thread;
			double bcov_tmp[3];
			i_perm_thread = (int *)malloc(*n * sizeof(int));
			Rank_thread = alloc_int_matrix((*n) + 1, (*n) + 1);
#pragma omp critical
			{
				for (k = 0; k < *n; k++)
				{
					i_perm_thread[k] = k;
				}
			}
#pragma omp for
			for (j_thread = 0; j_thread < (*R); j_thread++) {
#pragma omp critical
				{
					// TODO: I don't know why this command can only be runned in a single thread. But it makes R console temporarily available and seems to be correct.
					resample2(i_perm_thread, n);
				}
				initRank(*n, Rank_thread, xrank, yrank, i_perm_thread);
				U_Ball_Information(bcov_tmp, n, Rank_thread, lowxidx, higxidx, lowyidx, higyidx, i_perm_thread);
				permuted_bcov_weight0[j_thread] = bcov_tmp[0];
				permuted_bcov_weight_prob[j_thread] = bcov_tmp[1];
				permuted_bcov_weight_hhg[j_thread] = bcov_tmp[2];
			}
			free(i_perm_thread);
			free_int_matrix(Rank_thread, (*n) + 1, (*n) + 1);
		}

		pvalue[0] = compute_pvalue(bcov[0], permuted_bcov_weight0, *R);
		pvalue[1] = compute_pvalue(bcov[1], permuted_bcov_weight_prob, *R);
		pvalue[2] = compute_pvalue(bcov[2], permuted_bcov_weight_hhg, *R);

		free(permuted_bcov_weight0); free(permuted_bcov_weight_prob); free(permuted_bcov_weight_hhg);
	}


	free_int_matrix(lowxidx, *n, *n);
	free_int_matrix(higxidx, *n, *n);
	free_int_matrix(lowyidx, *n, *n);
	free_int_matrix(higyidx, *n, *n);
	free(xrank);
	free(yrank);
	
	return;
}


/////////////////////////////////////////////////////////////////
//////////// two functions below wrap key function //////////////
// bcov_stat return ball covariance statistic 
// bcov_test execute ball covariance statistic based test
/////////////////////////////////////////////////////////////////

void bcov_test(double *bcov, double *pvalue, double *x, double *y, int *n, int *R, int *dst, int *thread)
{
	//parallel method
	// if parallel_type == 1, we parallel the computation through statistics.
	// if parallel_type == 2, we parallel the computation through permutation.
	int single_thread = 1;
	int parallel_type = 2;
	if ((*n) >= 500)
	{
		parallel_type = 1;
	}
	if ((*R) <= 200)
	{
		*thread = single_thread;
	}
	if ((*dst)) {
		if (parallel_type == 2 && *thread > 1)
		{
		  // printf("BI parallel\n");
			BI_parallel(bcov, pvalue, x, y, n, R, thread);
		}
		else {
		  // printf("BI\n");
		  *thread = single_thread;
			BI(bcov, pvalue, x, y, n, R, thread);
		}
	}
	else {
		if (parallel_type == 2 && *thread > 1)
		{
			UBI_parallel(bcov, pvalue, x, y, n, R, thread);
		}
		else
		{
		  *thread = single_thread;
			UBI(bcov, pvalue, x, y, n, R, thread);
		}
	}
	return;
}
