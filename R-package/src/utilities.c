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
#include "math.h"
#include "utilities.h"
#include "stdio.h"
#include "stdlib.h"
#include "utilize_R.h"


 /*
 The rank computation (initRank, computeRank) is refer to: http://www.jmlr.org/papers/volume17/14-441/14-441.pdf [section 3.1.2]
 */
void computeRank(int n, int **Rank)
{
	int i, j;
	for (i = 1; i < n; i++)
		for (j = 1; j < n; j++)
			Rank[i][j] += (Rank[i][j - 1] + Rank[i - 1][j] - Rank[i - 1][j - 1]);
}

void initRank(int n, int **Rank, int *xrank, int *yrank, int *i_perm)
{
	int i, j;
	for (i = 0; i < n + 1; i++)
		for (j = 0; j < n + 1; j++)
			Rank[i][j] = 0;
	for (i = 0; i < n; i++)
		Rank[xrank[i] + 1][yrank[i_perm[i]] + 1] += 1;
	computeRank(n + 1, Rank);
}


//void initRank_surv(int n, int **Rank, int *xrank, int *yrank, int *i_perm)
//{
//	int i, j, old_x_rank = -1, old_t_rank = -1;
//	for (i = 0; i < n + 1; i++)
//		for (j = 0; j < n + 1; j++)
//			Rank[i][j] = 0;
//	//
//	for (i = 0; i < n; i++)
//	{
//		if (old_x_rank != (xrank[i] + 1) & old_t_rank != (yrank[i_perm[i]] + 1))
//		{
//			Rank[xrank[i] + 1][yrank[i_perm[i]] + 1] += 1;
//			old_x_rank = xrank[i] + 1;
//			old_t_rank = yrank[i_perm[i]] + 1;
//		}
//		else {
//			Rank[old_x_rank][old_t_rank] += 1;
//		}
//	}
//	printf("\n");
//	for (i = n - 1; i >= 0; i--)
//	{
//		for (int j = 0; j < (n); j++)
//			printf("%d ", Rank[i][j]);
//		printf("\n");
//	}
//	printf("\n");
//	computeRank(n + 1, Rank);
//}


/*
Input: n=6, zrank = []; z = [1, 2, 3, 4, 5, 5], zidx = [3, 1, 5, 2, 6, 4];
Output: zrank = [2, 4, 1, 5, 3, 5];
*/
void ranksort(int *n, int *zrank, double *z, int *zidx)
{
	int i, lastpos = 0;
	double lastval = -1.0;

	for (i = *n - 1; i >= 0; i--) {
		if (lastval != z[i])
			lastpos = i;
		lastval = z[i];
		zrank[zidx[i]] = lastpos;
	}
}

void quicksort(double *a, int *idx, int l, int u)
{
  int i, m, idx_temp;
  double a_temp;
  if (l >= u)
    return;
  m = l;
  for (i=l+1; i<=u; i++)
  {
    if (a[i] < a[l])
    {
      ++m;
      idx_temp = idx[m];
      idx[m] = idx[i];
      idx[i] = idx_temp;
      
      a_temp = a[m];
      a[m] = a[i];
      a[i] = a_temp;
    }
  }
  idx_temp = idx[l];
  idx[l] = idx[m];
  idx[m] = idx_temp;
  
  a_temp = a[l];
  a[l] = a[m];
  a[m] = a_temp;
  
  quicksort(a, idx, l, m-1);
  quicksort(a, idx, m+1, u);
}

void quicksort2(double *a, double *b, int *idx, int l, int u)
{
  int i, m, idx_temp;
  double a_temp;
  if (l >= u)
    return;
  m = l;
  for (i=l+1; i<=u; i++)
  {
    if (a[i] < a[l] || (a[i] == a[l] && b[i] > b[l]))
    {
      ++m;
      idx_temp = idx[m];
      idx[m] = idx[i];
      idx[i] = idx_temp;
      
      a_temp = a[m];
      a[m] = a[i];
      a[i] = a_temp;
      
      a_temp = b[m];
      b[m] = b[i];
      b[i] = a_temp;
    }
  }
  idx_temp = idx[l];
  idx[l] = idx[m];
  idx[m] = idx_temp;
  
  a_temp = a[l];
  a[l] = a[m];
  a[m] = a_temp;
  
  a_temp = b[l];
  b[l] = b[m];
  b[m] = a_temp;
  
  quicksort2(a, b, idx, l, m-1);
  quicksort2(a, b, idx, m+1, u);
}

double **alloc_matrix(int r, int c)
{
  /* allocate a matrix with r rows and c columns */
  int i;
  double **matrix;
  matrix = (double **) calloc(r, sizeof(double *));
  for (i = 0; i < r; i++)
    matrix[i] = (double *) calloc(c, sizeof(double));
  return matrix;
}

int **alloc_int_matrix(int r, int c)
{
  /* allocate a matrix with r rows and c columns */
  int i;
  int **matrix;
  matrix = (int **) calloc(r, sizeof(int *));
  for (i = 0; i < r; i++)
    matrix[i] = (int *) calloc(c, sizeof(int));
  return matrix;
}

void free_matrix(double **matrix, int r, int c)
{
  /* free a matrix with r rows and c columns */
  int i;
  for (i = 0; i < r; i++){
    free(matrix[i]);	
  }
  free(matrix);
}

void free_int_matrix(int **matrix, int r, int c)
{
  /* free a matrix with r rows and c columns */
  int i;
  for (i = 0; i < r; i++){
    free(matrix[i]);	
  }
  free(matrix);
}

void vector2matrix(double *x, double **y, int N, int d, int isroworder)
{
  /* copy a d-variate sample into a matrix, N samples in rows */
  int i, k;
  if (isroworder == 1) {
    for (k=0; k<d; k++)
      for (i=0; i<N; i++)
        y[i][k] = (*(x+i*d+k));
  }
  else {
    for (k=0; k<N; k++)
      for (i=0; i<d; i++)
        y[i][k] = (*(x+k*N+i));
  }
  return;
}

void Euclidean_distance(double *x, double **Dx, int n, int d)
{
  /*
   interpret x as an n by d matrix, in row order (n vectors in R^d)
   compute the Euclidean distance matrix Dx
   */
  int i, j, k, p, q;
  double dsum, dif;
  for (i=1; i<n; i++) {
    Dx[i][i] = 0.0;
    p = i*d;
    for (j=0; j<i; j++) {
      dsum = 0.0;
      q = j*d;
      for (k=0; k<d; k++) {
        dif = *(x+p+k) - *(x+q+k);
        dsum += dif*dif;
      }
      Dx[i][j] = Dx[j][i] = sqrt(dsum);
    }
  }
}

void distance(double *x, double *Dx, int *n, int *d)
{
  /*
   interpret x as an n by d matrix, in row order (n vectors in R^d)
   compute the Euclidean distance matrix Dx
   */
  int i, j, k, p, q;
  double dsum, dif;
  for (i=1; i<n[0]; i++) {
    p = i*d[0];
    for (j=0; j<i; j++) {
      dsum = 0.0;
      q = j*d[0];
      for (k=0; k<d[0]; k++) {
        dif = *(x+p+k) - *(x+q+k);
        dsum += dif*dif;
      }
      Dx[i*n[0]+j] = Dx[j*n[0]+i] = sqrt(dsum);
    }
  }
}


void resample(int *i_perm, int *i_perm_inv, int *n)
{
  int i, j, temp;
  for (i = *n - 1; i > 0; --i) {
    j = random_index2(i);
    temp = i_perm[j];
    i_perm[j] = i_perm[i];
    i_perm[i] = temp;
  }
  for (i = 0; i < *n; ++i) {
    i_perm_inv[i_perm[i]] = i;
  }
}


void resample2(int *i_perm, int *n)
{
  int i, j, temp;
  for (i = *n - 1; i > 0; --i) {
    // j = rand() % (i + 1);
    j = random_index2(i);
    temp = i_perm[j];
    i_perm[j] = i_perm[i];
    i_perm[i] = temp;
  }
}


/*
 * permute group index: i_perm
 */
void resample3(int *i_perm, int *i_perm_tmp, int n, int *n1)
{
	int i, j, temp, tmp0, tmp1;
  
	// permute step:
	for (i = n - 1; i > 0; --i) {
		// j = rand() % (i + 1);
		j = random_index2(i);
		temp = i_perm[j];
		i_perm[j] = i_perm[i];
		i_perm[i] = temp;
	}
  
	tmp0 = 0;
	tmp1 = 0;
	for (i = 0; i < n; i++) {
		if (i_perm[i] == 1) {
			i_perm_tmp[tmp0++] = i;
		}
		else {
			i_perm_tmp[*n1 + tmp1] = i;
			tmp1++;
		}
	}
}


/* Arrange the N elements of ARRAY in random order.
 Only effective if N is much smaller than RAND_MAX;
 if this may not be the case, use a better random
 number generator. */
void shuffle(int *array, int *N)
{
  // Rprintf("%d", RAND_MAX);  RAND_MAX = 32767;
  int n = *N;
  if (n > 1) 
  {
    int i, j, t;
    for (i = 0; i < n - 1; i++) 
    {
      j = random_index(n, i);
      t = array[j];
      array[j] = array[i];
      array[i] = t;
    }
  }
}


void shuffle_value(double *array, int *N)
{
  // Rprintf("%d", RAND_MAX);  RAND_MAX = 32767;
  int n = *N;
  if (n > 1) 
  {
    int i, j;
    double tmp;
    for (i = 0; i < n - 1; i++) 
    {
      j = random_index(n, i);
      tmp = array[j];
      array[j] = array[i];
      array[i] = tmp;
    }
  }
}


int pending_interrupt() {
  int interrupt_status = 0;
  interrupt_status = pending_interrupt_status();
  return interrupt_status;
}


void print_stop_message()
{
  print_stop_message_internal();
  return;
}