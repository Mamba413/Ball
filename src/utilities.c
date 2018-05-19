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
#include "utilize_cross.h"


void swap(double *x, double *y) {
	double t = *x;
	*x = *y;
	*y = t;
}

void quick_sort_recursive(double *arr, int start, int end) {
	if (start >= end)
		return;
	double mid = arr[end];
	int left = start, right = end - 1;
	while (left < right) {
		while (arr[left] < mid && left < right)
			left++;
		while (arr[right] >= mid && left < right)
			right--;
		swap(&arr[left], &arr[right]);
	}
	if (arr[left] >= arr[end])
		swap(&arr[left], &arr[end]);
	else
		left++;
	if (left)
		quick_sort_recursive(arr, start, left - 1);
	quick_sort_recursive(arr, left + 1, end);
}

/*
quick sort function for finding Max K-1 two sample ball divergence value
*/
void quick_sort(double *arr, int len) {
	quick_sort_recursive(arr, 0, len - 1);
}


double compute_pvalue(double ball_stat_value, double *permuted_stat, int R)
{
	double larger_num = 0.0;
	for (int i = 0; i < R; i++)
	{
		//printf("lastest permute value: %f\n", permuted_stat[i]);
		if (permuted_stat[i] > ball_stat_value)
		{
			larger_num += 1.0;
		}
	}
	double R_double = R;
	double p_value = (1.0 + larger_num) / (1.0 + R_double);
	return(p_value);
}

void Merge(int *permutation, int *source, int *inversion_count, int dim, int n)
{
	int *left = (int *)malloc(n * sizeof(int));
	int *right = (int *)malloc(n * sizeof(int));
	int *left_source = (int *)malloc(n * sizeof(int));
	int *right_source = (int *)malloc(n * sizeof(int));
	int left_index = 0, right_index = 0;
	int i, half_dim = dim / 2, nleft = half_dim, nright = dim - half_dim;
	for (i = 0; i < half_dim; i++) {
		left[i] = permutation[i];
		left_source[i] = source[i];
		right[i] = permutation[i + half_dim];
		right_source[i] = source[i + half_dim];
	}

	if (nleft < nright) {
		right[i] = permutation[i + half_dim];
		right_source[i] = source[i + half_dim];
	}

	for (i = 0; i < dim; i++) {
		if ((left_index < half_dim) && (right_index < dim - half_dim)) {
			if (left[left_index] <= right[right_index]) { // I added "=" in order to support ties
				permutation[i] = left[left_index];
				source[i] = left_source[left_index];
				left_index++;
			}
			else {
				permutation[i] = right[right_index];
				source[i] = right_source[right_index];
				inversion_count[source[i]] += (half_dim - left_index);
				right_index++;
			}
		}
		else {
			if (left_index < half_dim) {
				permutation[i] = left[left_index];
				source[i] = left_source[left_index];
				left_index++;
			}

			if (right_index < dim - half_dim) {
				permutation[i] = right[right_index];
				source[i] = right_source[right_index];
				right_index++;
			}
		}
	}
	free(left);
	free(right);
	free(left_source);
	free(right_source);
}


int Inversions(int *permutation, int *source, int *inversion_count, int dim, int n)
{
	if (dim == 1)
		return 0;
	else {
		Inversions(permutation, source, inversion_count, dim / 2, n);
		Inversions(&permutation[dim / 2], &source[dim / 2], inversion_count, dim - dim / 2, n);
		Merge(permutation, source, inversion_count, dim, n);
	}
	return 0;
}


 /*
 The algorithm like the quick-sort. It finds the upper and lower in turn.
 For each i, the computation complexity is O(n);
 Thus, this function take O(n^2) times.

 Algorithm detail:
 For i=1, ..., n, carry out:
 init the low and upper index: l=1, u=N
 while (l <= u):
 if (z[u] - z[i]) >= (z[i] - z[l]) ==> for pair (i, u), upper index: u, lower index: l ==> update upper index: u = u - 1
 if (z[u] - z[i]) < (z[i] - z[l]) ==> for pair (i, l), upper index: u, lower index: l ==> update lower index: l = l - 1

 Input: n=6, z = [1, 2, 3, 4, 5, 5], zidx = [3, 1, 5, 2, 6, 4], lowzidx = [], higzidx = [];
 Output: lowzidx, higzidx;
 */
void createidx(int *n, int *zidx, double *z, int **lowzidx, int **higzidx)
{
	int i, zi, ileft, iright, jleft, jright, lowpos, higpos;
	double lastval, tmp1, tmp2, tmp;
	for (i = 0; i < *n; i++) {
		zi = zidx[i];
		lowpos = 1;
		higpos = *n;
		ileft = 0;
		iright = *n - 1;
		jleft = 0;
		jright = 0;

		tmp1 = z[iright] - z[i];
		tmp2 = z[i] - z[ileft];
		if (tmp1 > tmp2) {
			lastval = tmp1;
			lowzidx[zi][zidx[iright]] = lowpos;
			higzidx[zi][zidx[iright]] = higpos;
			iright--;
			jright++;
		}
		else {
			lastval = tmp2;
			if (ileft == i) {
				lowzidx[zi][zidx[iright]] = lowpos;
				higzidx[zi][zidx[iright]] = higpos;
				iright--;
				jright++;
			}
			else {
				lowzidx[zi][zidx[ileft]] = lowpos;
				higzidx[zi][zidx[ileft]] = higpos;
				ileft++;
				jleft++;
			}
		}

		while (ileft <= iright) {
			tmp1 = z[iright] - z[i];
			tmp2 = z[i] - z[ileft];
			tmp = MAX(tmp1, tmp2);
			while (lastval == tmp) {
				if (tmp1 > tmp2) {
					lowzidx[zi][zidx[iright]] = lowpos;
					higzidx[zi][zidx[iright]] = higpos;
					iright--;
					jright++;
				}
				else {
					if (ileft == i) {
						lowzidx[zi][zidx[iright]] = lowpos;
						higzidx[zi][zidx[iright]] = higpos;
						iright--;
						jright++;
					}
					else {
						lowzidx[zi][zidx[ileft]] = lowpos;
						higzidx[zi][zidx[ileft]] = higpos;
						ileft++;
						jleft++;
					}
				}
				if (iright < ileft)
					break;
				tmp1 = z[iright] - z[i];
				tmp2 = z[i] - z[ileft];
				tmp = MAX(tmp1, tmp2);
			}
			if (iright < ileft)
				break;
			lowpos += jleft;
			higpos -= jright;
			jleft = 0;
			jright = 0;
			if (tmp1 > tmp2) {
				lastval = tmp1;
				lowzidx[zi][zidx[iright]] = lowpos;
				higzidx[zi][zidx[iright]] = higpos;
				iright--;
				jright++;
			}
			else {
				lastval = tmp2;
				if (ileft == i) {
					lowzidx[zi][zidx[iright]] = lowpos;
					higzidx[zi][zidx[iright]] = higpos;
					iright--;
					jright++;
				}
				else {
					lowzidx[zi][zidx[ileft]] = lowpos;
					higzidx[zi][zidx[ileft]] = higpos;
					ileft++;
					jleft++;
				}
			}
		}
	}
}


void sort(int *n, int *zidx, double *z, int **dzidx)
{
	// the z[i] is the center
	int i, j, zi, ileft, iright, lastpos;
	double lastval, tmp;
	for (i = 0; i < *n; i++) {
		zi = zidx[i];
		j = *n - 1;
		lastpos = *n - 1;
		ileft = 0;
		iright = *n - 1;
		lastval = -1.0;

		while ((ileft != i) || (iright != i)) {
			if (i == ileft) {
				tmp = z[iright] - z[i];
				if (lastval != tmp)
					lastpos = j;
				dzidx[zi][zidx[iright]] = lastpos;
				lastval = tmp;
				iright--;
			}
			else if (i == iright) {
				tmp = z[i] - z[ileft];
				if (lastval != tmp)
					lastpos = j;
				dzidx[zi][zidx[ileft]] = lastpos;
				lastval = tmp;
				ileft++;
			}
			else {
				if (z[i] - z[ileft]>z[iright] - z[i]) {
					tmp = z[i] - z[ileft];
					if (lastval != tmp)
						lastpos = j;
					dzidx[zi][zidx[ileft]] = lastpos;
					lastval = tmp;
					ileft++;
				}
				else {
					tmp = z[iright] - z[i];
					if (lastval != tmp)
						lastpos = j;
					dzidx[zi][zidx[iright]] = lastpos;
					lastval = tmp;
					iright--;
				}
			}
			j--;
		}

		if (lastval == 0)
			dzidx[zi][zi] = lastpos;
		else
			dzidx[zi][zi] = 0;
	}
}


void ranksort2(int n, int **Rxy, double **Dxy, int **Ixy)
{
	int i, j, lastpos = n - 1;
	double lastval;
	for (i = 0; i < n; i++) {
		lastval = -1;
		for (j = n - 1; j >= 0; j--) {
			if (lastval != Dxy[i][j])
				lastpos = j;
			lastval = Dxy[i][j];
			Rxy[i][Ixy[i][j]] = lastpos;
		}
	}
}


void Findx2(int *Rxy, int *Ixy, int *i_perm, int *n1, int *n2, int *Rx)
{
	int j, lastpos, lastval, n, tmp;
	n = *n1 + *n2;

	lastpos = *n1 - 1;
	if (i_perm[Ixy[n - 1]] == 1) {
		tmp = 1;
		lastval = Rxy[Ixy[n - 1]];
	}
	else {
		tmp = 0;
		lastval = -1;
	}
	Rx[Ixy[n - 1]] = lastpos;

	for (j = n - 2; j >= 0; j--) {
		if (i_perm[Ixy[j]] == 1) {
			if (lastval != Rxy[Ixy[j]]) {
				lastpos -= tmp;
				tmp = 0;
			}
			tmp++;
			lastval = Rxy[Ixy[j]];
			Rx[Ixy[j]] = lastpos;
		}
		else {
			if (Rxy[Ixy[j]] == Rxy[Ixy[j + 1]])
				Rx[Ixy[j]] = Rx[Ixy[j + 1]];
			else
				Rx[Ixy[j]] = lastpos - tmp;
		}
	}
}


void Findx(int **Rxy, int **Ixy, int *i_perm, int *n1, int *n2, int **Rx)
{
	int i, n;
	n = *n1 + *n2;
	for (i = 0; i < n; i++)
		Findx2(Rxy[i], Ixy[i], i_perm, n1, n2, Rx[i]);
}


void ranksort3(int n, int *xyidx, double *xy, int **Rxy, int **Ixy)
{
	int i, j, ileft, iright, lastpos;
	double lastval;
	for (i = 0; i < n; i++) {
		lastval = -1;
		ileft = 0;
		iright = n - 1;
		j = n - 1;
		lastpos = n - 1;
		while (ileft<iright) {
			if ((lastval != xy[i] - xy[ileft]) && (lastval != xy[iright] - xy[i]))
				lastpos = j;
			if (ileft == i) {
				lastval = xy[iright] - xy[i];
				Ixy[xyidx[i]][j] = xyidx[iright];
				Rxy[xyidx[i]][xyidx[iright]] = lastpos;
				iright--;
			}
			else if (iright == i) {
				lastval = xy[i] - xy[ileft];
				Ixy[xyidx[i]][j] = xyidx[ileft];
				Rxy[xyidx[i]][xyidx[ileft]] = lastpos;
				ileft++;
			}
			else {
				if (xy[i] - xy[ileft] > xy[iright] - xy[i]) {
					lastval = xy[i] - xy[ileft];
					Ixy[xyidx[i]][j] = xyidx[ileft];
					Rxy[xyidx[i]][xyidx[ileft]] = lastpos;
					ileft++;
				}
				else {
					lastval = xy[iright] - xy[i];
					Ixy[xyidx[i]][j] = xyidx[iright];
					Rxy[xyidx[i]][xyidx[iright]] = lastpos;
					iright--;
				}
			}
			j--;
		}
		Ixy[xyidx[i]][0] = xyidx[i];
		if (lastval == 0)
			Rxy[xyidx[i]][xyidx[i]] = lastpos;
		else
			Rxy[xyidx[i]][xyidx[i]] = 0;
	}
}


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

/*
Apply quicksort algorithm to a, and the index exchange result is recorded in idx.
*/
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

/*
Apply quicksort algorithm to a and b
After sorted, a is increasing while b is decreasing, and idx recode the index exchange result
Example:
Input:
a = [2, 2, 1, 3]
b = [1, 2, 3, 4]
idex = [1, 2, 3, 4]
Output:
a = [1, 2, 2, 3]
b = [3, 1, 2, 4]
c = [3, 1, 2, 4]
*/
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

/*
i_perm takes value 1 to n, i.e., index of sample
This function permute i_perm array to achieve permutation
*/
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