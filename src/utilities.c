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
#include "string.h"
#include "stdio.h"
#include "stdlib.h"
#include "utilize_cross.h"

#ifdef R_BUILD
#include "R.h"
#include "Rinternals.h"
#else

#include "time.h"

#endif

void swap(double *x, double *y) {
    double t = *x;
    *x = *y;
    *y = t;
}

/**
 *
 * @param group_relative_location : [0, 2, 4, 1, 3, 5]
 * @param group :  [0, 1, 2, 0, 1, 2]
 * @param K : 3
 */
void find_group_relative_location(int *group_relative_location, int *group, int *cumsum_size, int num, int K) {
    int *init_group_relative_location = (int *) malloc(K * sizeof(int));
    for (int k = 0; k < K; ++k) {
        init_group_relative_location[k] = 0;
    }
    for (int i = 0; i < num; ++i) {
        for (int j = 0; j < K; ++j) {
            if (group[i] == j) {
                group_relative_location[i] = cumsum_size[j] + init_group_relative_location[j];
                init_group_relative_location[j] += 1;
                break;
            }
        }
    }
}

// input:
// size: an array contain sample size in each group
// cumulate_size: the cumulative sums of size vector
// k: group number
void compute_cumsum_size(int *cumulate_size, int *size, int *k) {
    int i;
    for (i = 0; i < (*k); i++) {
        if (i == 0) {
            cumulate_size[i] = 0;
        } else {
            cumulate_size[i] = cumulate_size[i - 1] + size[i - 1];
        }
    }
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


double compute_pvalue(double ball_stat_value, double *permuted_stat, int R) {
    double larger_num = 0.0;
    for (int i = 0; i < R; i++) {
        if (permuted_stat[i] > ball_stat_value) {
            larger_num += 1.0;
        }
    }
    double p_value = (1.0 + larger_num) / (1.0 + R);
    return p_value;
}

void Merge(int *permutation, int *source, int *inversion_count, int dim, int n) {
    int *left = (int *) malloc(n * sizeof(int));
    int *right = (int *) malloc(n * sizeof(int));
    int *left_source = (int *) malloc(n * sizeof(int));
    int *right_source = (int *) malloc(n * sizeof(int));
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
            } else {
                permutation[i] = right[right_index];
                source[i] = right_source[right_index];
                inversion_count[source[i]] += (half_dim - left_index);
                right_index++;
            }
        } else {
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


int Inversions(int *permutation, int *source, int *inversion_count, int dim, int n) {
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
void createidx(int *n, int *zidx, double *z, int **lowzidx, int **higzidx) {
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
        } else {
            lastval = tmp2;
            if (ileft == i) {
                lowzidx[zi][zidx[iright]] = lowpos;
                higzidx[zi][zidx[iright]] = higpos;
                iright--;
                jright++;
            } else {
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
                } else {
                    if (ileft == i) {
                        lowzidx[zi][zidx[iright]] = lowpos;
                        higzidx[zi][zidx[iright]] = higpos;
                        iright--;
                        jright++;
                    } else {
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
            } else {
                lastval = tmp2;
                if (ileft == i) {
                    lowzidx[zi][zidx[iright]] = lowpos;
                    higzidx[zi][zidx[iright]] = higpos;
                    iright--;
                    jright++;
                } else {
                    lowzidx[zi][zidx[ileft]] = lowpos;
                    higzidx[zi][zidx[ileft]] = higpos;
                    ileft++;
                    jleft++;
                }
            }
        }
    }
}


void sort(int *n, int *zidx, double *z, int **dzidx) {
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
            } else if (i == iright) {
                tmp = z[i] - z[ileft];
                if (lastval != tmp)
                    lastpos = j;
                dzidx[zi][zidx[ileft]] = lastpos;
                lastval = tmp;
                ileft++;
            } else {
                if (z[i] - z[ileft] > z[iright] - z[i]) {
                    tmp = z[i] - z[ileft];
                    if (lastval != tmp)
                        lastpos = j;
                    dzidx[zi][zidx[ileft]] = lastpos;
                    lastval = tmp;
                    ileft++;
                } else {
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

void ranksort2(int n, int **Rxy, double **Dxy, int **Ixy) {
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


/*
This function try to derive the rank of Rx belong to group 1.
tmp, lastpos, lastval are introduced to resolve the situtation when ties appear
*/
void Findx2(int *Rxy, int *Ixy, int *i_perm, int *n1, int *n2, int *Rx) {
    int j, lastpos, lastval, n, tmp;
    n = *n1 + *n2;

    //printf("---------------\n");
    //for (int i = 0; i < n; i++)
    //{
    //	printf("%d ", Rxy[i]);
    //}
    //printf("\n");

    //for (int i = 0; i < n; i++)
    //{
    //	printf("%d ", Ixy[i]);
    //}
    //printf("\n");

    //for (int i = 0; i < n; i++)
    //{
    //	printf("%d ", i_perm[i]);
    //}
    //printf("\n");

    //for (int i = 0; i < n; i++)
    //{
    //	printf("%d ", Rx[i]);
    //}
    //printf("\n");

    lastpos = *n1 - 1;
    Rx[Ixy[n - 1]] = lastpos;
    if (i_perm[Ixy[n - 1]] == 1) {
        tmp = 1;
        lastval = Rxy[Ixy[n - 1]];
    } else {
        tmp = 0;
        lastval = -1;
    }

    for (j = n - 2; j >= 0; j--) {
        //for (int i = 0; i < n; i++)
        //{
        //	printf("%d ", Rx[i]);
        //}
        //printf("\n");
        if (i_perm[Ixy[j]] == 1) {
            if (lastval != Rxy[Ixy[j]]) {
                lastpos -= tmp;
                tmp = 0;
            }
            tmp++;
            lastval = Rxy[Ixy[j]];
            Rx[Ixy[j]] = lastpos;
        } else {
            if (Rxy[Ixy[j]] == Rxy[Ixy[j + 1]])
                Rx[Ixy[j]] = Rx[Ixy[j + 1]];
            else
                Rx[Ixy[j]] = lastpos - tmp;
        }
    }
}


void Findx(int **Rxy, int **Ixy, int *i_perm, int *n1, int *n2, int **Rx) {
    int i, n;
    n = *n1 + *n2;
    for (i = 0; i < n; i++)
        Findx2(Rxy[i], Ixy[i], i_perm, n1, n2, Rx[i]);
}


void ranksort3(int n, int *xyidx, double *xy, int **Rxy, int **Ixy) {
    int i, j, ileft, iright, lastpos;
    double lastval;
    for (i = 0; i < n; i++) {
        lastval = -1;
        ileft = 0;
        iright = n - 1;
        j = n - 1;
        lastpos = n - 1;
        while (ileft < iright) {
            if ((lastval != xy[i] - xy[ileft]) && (lastval != xy[iright] - xy[i]))
                lastpos = j;
            if (ileft == i) {
                lastval = xy[iright] - xy[i];
                Ixy[xyidx[i]][j] = xyidx[iright];
                Rxy[xyidx[i]][xyidx[iright]] = lastpos;
                iright--;
            } else if (iright == i) {
                lastval = xy[i] - xy[ileft];
                Ixy[xyidx[i]][j] = xyidx[ileft];
                Rxy[xyidx[i]][xyidx[ileft]] = lastpos;
                ileft++;
            } else {
                if (xy[i] - xy[ileft] > xy[iright] - xy[i]) {
                    lastval = xy[i] - xy[ileft];
                    Ixy[xyidx[i]][j] = xyidx[ileft];
                    Rxy[xyidx[i]][xyidx[ileft]] = lastpos;
                    ileft++;
                } else {
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
void computeRank(int n, int **Rank) {
    int i, j;
    for (i = 1; i < n; i++)
        for (j = 1; j < n; j++)
            Rank[i][j] += (Rank[i][j - 1] + Rank[i - 1][j] - Rank[i - 1][j - 1]);
}

void initRank(int n, int **Rank, int *xrank, int *yrank, int *i_perm) {
    int i, j;
    for (i = 0; i < n + 1; i++)
        for (j = 0; j < n + 1; j++)
            Rank[i][j] = 0;
    for (i = 0; i < n; i++)
        Rank[xrank[i] + 1][yrank[i_perm[i]] + 1] += 1;
    computeRank(n + 1, Rank);
}

void initRank_bcor(int n, int **Rank, int *xrank, int *yrank) {
    int i, j;
    for (i = 0; i < n + 1; i++)
        for (j = 0; j < n + 1; j++)
            Rank[i][j] = 0;
    for (i = 0; i < n; i++)
        Rank[xrank[i] + 1][yrank[i] + 1] += 1;
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
void ranksort(int *n, int *zrank, double *z, int *zidx) {
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
void quicksort(double *a, int *idx, int l, int u) {
    int i, m, idx_temp;
    double a_temp;
    if (l >= u)
        return;
    m = l;
    for (i = l + 1; i <= u; i++) {
        if (a[i] < a[l]) {
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

    quicksort(a, idx, l, m - 1);
    quicksort(a, idx, m + 1, u);
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
void quicksort2(double *a, double *b, int *idx, int l, int u) {
    int i, m, idx_temp;
    double a_temp;
    if (l >= u)
        return;
    m = l;
    for (i = l + 1; i <= u; i++) {
        if (a[i] < a[l] || (a[i] == a[l] && b[i] > b[l])) {
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

    quicksort2(a, b, idx, l, m - 1);
    quicksort2(a, b, idx, m + 1, u);
}

double **alloc_matrix(int r, int c) {
    /* allocate a matrix with r rows and c columns */
    int i;
    double **matrix;
    matrix = (double **) calloc(r, sizeof(double *));
    for (i = 0; i < r; i++)
        matrix[i] = (double *) calloc(c, sizeof(double));
    return matrix;
}

double ***alloc_3d_matrix(int r, int c, int h) {
    /* allocate a 3D matrix with r rows, c columns, h levels */
    double ***arr3D;
    int i, j;

    arr3D = (double ***) malloc(r * sizeof(double **));

    for (i = 0; i < r; i++) {
        arr3D[i] = (double **) malloc(c * sizeof(double *));
        for (j = 0; j < c; j++) {
            arr3D[i][j] = (double *) malloc(h * sizeof(double));
        }
    }
    return arr3D;
}

int **alloc_int_matrix(int r, int c) {
    /* allocate a matrix with r rows and c columns */
    int i;
    int **matrix;
    matrix = (int **) calloc(r, sizeof(int *));
    for (i = 0; i < r; i++)
        matrix[i] = (int *) calloc(c, sizeof(int));
    return matrix;
}

int ***alloc_3d_int_matrix(int r, int c, int h) {
    /* allocate a 3D matrix with r rows, c columns, h levels */
    int ***arr3D;
    int i, j;

    arr3D = (int ***) malloc(r * sizeof(int **));

    for (i = 0; i < r; i++) {
        arr3D[i] = (int **) malloc(c * sizeof(int *));
        for (j = 0; j < c; j++) {
            arr3D[i][j] = (int *) malloc(h * sizeof(int));
        }
    }
    return arr3D;
}

int ***alloc_int_square_matrix_list(int* size, int number) {
    int ***arr3D = (int ***) malloc(number * sizeof(int **));
    for (int i = 0; i < number; ++i) {
        arr3D[i] = alloc_int_matrix(size[i], size[i]);
    }
    return arr3D;
}

void free_matrix(double **matrix, int r, int c) {
    /* free a matrix with r rows and c columns */
    int i;
    for (i = 0; i < r; i++) {
        free(matrix[i]);
    }
    free(matrix);
}

void free_3d_matrix(double ***arr3D, int r, int c) {
    int i, j;

    for (i = 0; i < r; i++) {
        for (j = 0; j < c; j++) {
            free(arr3D[i][j]);
        }
        free(arr3D[i]);
    }
    free(arr3D);
}

void free_int_matrix(int **matrix, int r, int c) {
    /* free a matrix with r rows and c columns */
    int i;
    for (i = 0; i < r; i++) {
        free(matrix[i]);
    }
    free(matrix);
}

void free_3d_int_matrix(int ***arr3D, int r, int c) {
    int i, j;

    for (i = 0; i < r; i++) {
        for (j = 0; j < c; j++) {
            free(arr3D[i][j]);
        }
        free(arr3D[i]);
    }
    free(arr3D);
}

void free_int_square_matrix_list(int ***arr3d, int* size, int num) {
    for (int i = 0; i < num; ++i) {
        free_int_matrix(arr3d[i], size[i], size[i]);
    }
    free(arr3d);
}

void vector2matrix(double *x, double **y, int N, int d, int isroworder) {
    /* copy a d-variate sample into a matrix, N samples in rows */
    int i, k;
    if (isroworder == 1) {
        for (k = 0; k < d; k++) {
            for (i = 0; i < N; i++) {
                y[i][k] = (*(x + i * d + k));
            }
        }
    } else {
        for (k = 0; k < N; k++) {
            for (i = 0; i < d; i++) {
                y[i][k] = (*(x + k * N + i));
            }
        }
    }
}

void distance2matrix(double *distance, double **distance_matrix, int n) {
    int k = 0;
    for (int i = 0; i < n; ++i) {
        distance_matrix[i][i] = 0.0;
        for (int j = (i + 1); j < n; ++j) {
            distance_matrix[j][i] = distance_matrix[i][j] = distance[k];
            k++;
        }
    }
}


void vector2matrix3d(double *x, double ***y, int r, int c, int h, int isroworder) {
    /* copy a d-variate sample into a matrix, N samples in rows */
    int i, j, k;
    int index = 0;
    if (isroworder == 1) {
        for (k = 0; k < h; k++) {
            for (j = 0; j < c; j++) {
                for (i = 0; i < r; i++) {
                    y[i][j][k] = x[index];
                    index += 1;
                }
            }
        }
    }
}

void distance2matrix3d(double *distance, double ***distance_matrix3d, int n, int v) {
    int s = 0;
    for (int k = 0; k < v; ++k) {
        for (int i = 0; i < n; ++i) {
            distance_matrix3d[i][i][k] = 0.0;
            for (int j = (i + 1); j < n; ++j) {
                distance_matrix3d[j][i][k] = distance_matrix3d[i][j][k] = distance[s];
                s++;
            }
        }
    }
}

void rank_matrix_3d(double ***Dx, int n, int k, int ***Rx) {
    int i, j, h, *x_part_rank;
    double *x_part;
    x_part = (double *) malloc(n * sizeof(double));
    x_part_rank = (int *) malloc(n * sizeof(int));
    for (h = 0; h < k; h++) {
        for (i = 0; i < n; i++) {
            for (j = 0; j < n; j++) {
                x_part[j] = Dx[i][j][h];
            }
            quick_rank(x_part, x_part_rank, n);
            for (j = 0; j < n; j++) {
                Rx[i][j][h] = x_part_rank[j];
            }
        }
    }
    free(x_part);
    free(x_part_rank);
}

/*
 * double aaa[6] = {1.2, 1.3, 1.3, 1.3, 1.3, 0.9};
 * int aaa_r[6] = {0, 0, 0, 0, 0, 0};
 * quick_rank(aaa, aaa_r, 6);
 */
void quick_rank(double *x, int *r, int n) {
    int *x_index, rank_value = n, i_loc, tmp = 1;
    double *x_cpy;
    x_index = (int *) malloc(n * sizeof(int));
    x_cpy = (double *) malloc(n * sizeof(double));
    for (int j = 0; j < n; j++) { x_index[j] = j; }
    for (int j = 0; j < n; j++) { x_cpy[j] = x[j]; }
    quicksort(x_cpy, x_index, 0, n - 1);
    r[x_index[n - 1]] = rank_value;
    for (int i = n - 2; 0 <= i; i--) {
        i_loc = x_index[i];
        if (x[i_loc] == x[x_index[i + 1]]) {
            r[i_loc] = rank_value;
            tmp++;
        } else {
            rank_value = rank_value - tmp;
            r[i_loc] = rank_value;
            tmp = 1;
        }
    }
    free(x_index);
    free(x_cpy);
}

void Euclidean_distance(double *x, double **Dx, int n, int d) {
    /*
     interpret x as an n by d matrix, in row order (n vectors in R^d)
     compute the Euclidean distance matrix Dx
     */
    int i, j, k, p, q;
    double dsum, dif;
    for (i = 1; i < n; i++) {
        Dx[i][i] = 0.0;
        p = i * d;
        for (j = 0; j < i; j++) {
            dsum = 0.0;
            q = j * d;
            for (k = 0; k < d; k++) {
                dif = *(x + p + k) - *(x + q + k);
                dsum += dif * dif;
            }
            Dx[i][j] = Dx[j][i] = sqrt(dsum);
        }
    }
}

void distance(double *x, double *Dx, int *n, int *d) {
    /*
     interpret x as an n by d matrix, in row order (n vectors in R^d)
     compute the Euclidean distance matrix Dx
     */
    int i, j, k, p, q;
    double dsum, dif;
    for (i = 1; i < n[0]; i++) {
        p = i * d[0];
        for (j = 0; j < i; j++) {
            dsum = 0.0;
            q = j * d[0];
            for (k = 0; k < d[0]; k++) {
                dif = *(x + p + k) - *(x + q + k);
                dsum += dif * dif;
            }
            Dx[i * n[0] + j] = Dx[j * n[0] + i] = sqrt(dsum);
        }
    }
}

/*
i_perm takes value 1 to n, i.e., index of sample
This function permute i_perm array to achieve permutation
*/
void resample(int *i_perm, int *i_perm_inv, int *n) {
    int i, j, temp;
#ifdef R_BUILD
    GetRNGstate();
    for (i = *n - 1; i > 0; --i) {
        j = ((int) round(RAND_MAX * unif_rand())) % (i + 1);
        temp = i_perm[j];
        i_perm[j] = i_perm[i];
        i_perm[i] = temp;
    }
    PutRNGstate();
#else
    srand((unsigned) time(NULL));
    for (i = *n - 1; i > 0; --i) {
        j = rand() % (i + 1);
        temp = i_perm[j];
        i_perm[j] = i_perm[i];
        i_perm[i] = temp;
    }
#endif
    for (i = 0; i < *n; ++i) {
        i_perm_inv[i_perm[i]] = i;
    }
}

void shuffle_indicator_inv_matrix(int **i_perm_matrix, int **i_perm_matrix_inv, int *init_perm, int *init_perm_inv,
                                  int num_permutation, int num) {
    int k, temp;
#ifdef R_BUILD
    GetRNGstate();
    for (int i = 0; i < num_permutation; i++) {
        for (int j = num - 1; j > 0; --j) {
            k = ((int) round(RAND_MAX * unif_rand())) % (j + 1);
            temp = init_perm[k];
            init_perm[k] = init_perm[j];
            init_perm[j] = temp;
        }
        for (int j = 0; j < num; ++j) {
            init_perm_inv[init_perm[j]] = j;
        }
        memcpy(i_perm_matrix[i], init_perm, num * sizeof(int));
        memcpy(i_perm_matrix_inv[i], init_perm_inv, num * sizeof(int));
    }
    PutRNGstate();
#else
    srand((unsigned) time(NULL));
    for (int i = 0; i < num_permutation; i++) {
        for (int j = num - 1; j > 0; --j) {
            k = rand() % (j + 1);
            temp = init_perm[k];
            init_perm[k] = init_perm[j];
            init_perm[j] = temp;
        }
        for (int j = 0; j < num; ++j) {
            init_perm_inv[init_perm[j]] = j;
        }
        memcpy(i_perm_matrix[i], init_perm, num * sizeof(int));
        memcpy(i_perm_matrix_inv[i], init_perm_inv, num * sizeof(int));
    }
#endif
}

void resample_matrix(int **i_perm, int *r, int *c) {
    int i, j, k, temp;
    for (k = 0; k < *r; k++) {
        for (i = *c - 1; i > 0; --i) {
            j = random_index2(i);
            temp = i_perm[k][j];
            i_perm[k][j] = i_perm[k][i];
            i_perm[k][i] = temp;
        }
    }
}

void resample_matrix_3d(int ***i_perm, int **init_perm, int *h, int *r, int *c) {
    int j, temp;
#ifdef R_BUILD
    GetRNGstate();
    for (int l = 0; l < *h; ++l) {
        for (int k = 0; k < *r; k++) {
            for (int i = *c - 1; i > 0; --i) {
                j = ((int) round(RAND_MAX * unif_rand())) % (i + 1);
                temp = init_perm[k][j];
                init_perm[k][j] = init_perm[k][i];
                init_perm[k][i] = temp;
            }
        }
        for (int k = 0; k < *r; ++k) {
            memcpy(i_perm[l][k], init_perm[k], (*c * sizeof(int)));
        }
    }
    PutRNGstate();
#else
    srand((unsigned) time(NULL));
    for (int l = 0; l < *h; ++l) {
        for (int k = 0; k < *r; k++) {
            for (int i = *c - 1; i > 0; --i) {
                j = rand() % (i + 1);
                temp = init_perm[k][j];
                init_perm[k][j] = init_perm[k][i];
                init_perm[k][i] = temp;
            }
        }
        for (int k = 0; k < *r; ++k) {
            memcpy(i_perm[l][k], init_perm[k], (*c * sizeof(int)));
        }
    }
#endif
}

void resample2(int *i_perm, int *n) {
    int i, j, temp;
#ifdef R_BUILD
    GetRNGstate();
    for (i = *n - 1; i > 0; --i) {
        j = ((int) round(RAND_MAX * unif_rand())) % (i + 1);
        temp = i_perm[j];
        i_perm[j] = i_perm[i];
        i_perm[i] = temp;
    }
    PutRNGstate();
#else
    srand((unsigned) time(NULL));
    for (i = *n - 1; i > 0; --i) {
        j = rand() % (i + 1);
        temp = i_perm[j];
        i_perm[j] = i_perm[i];
        i_perm[i] = temp;
    }
#endif
}

void resample2_matrix(int **i_perm, int *init_perm, int num_permutation, int n) {
    int i, j, temp;
#ifdef R_BUILD
    GetRNGstate();
    for (int r = 0; r < num_permutation; r++) {
        for (i = n - 1; i > 0; --i) {
            j = ((int) round(RAND_MAX * unif_rand())) % (i + 1);
            temp = init_perm[j];
            init_perm[j] = init_perm[i];
            init_perm[i] = temp;
        }
        memcpy(i_perm[r], init_perm, n * sizeof(int));
    }
    PutRNGstate();
#else
    srand((unsigned) time(NULL));
    for (int r = 0; r < num_permutation; r++) {
        for (i = n - 1; i > 0; --i) {
            j = rand() % (i + 1);
            temp = init_perm[j];
            init_perm[j] = init_perm[i];
            init_perm[i] = temp;
        }
        memcpy(i_perm[r], init_perm, n * sizeof(int));
    }
#endif
}

void shuffle_indicator_matrix(int **i_perm_matrix, int *init_perm, int num_permutation, int num) {
    int k, temp;
#ifdef R_BUILD
    GetRNGstate();
    for (int i = 0; i < num_permutation; i++) {
        for (int j = num - 1; j > 0; --j) {
            k = ((int) round(RAND_MAX * unif_rand())) % (j + 1);
            temp = init_perm[k];
            init_perm[k] = init_perm[j];
            init_perm[j] = temp;
        }
        memcpy(i_perm_matrix[i], init_perm, num * sizeof(int));
    }
    PutRNGstate();
#else
    srand((unsigned) time(NULL));
    for (int i = 0; i < num_permutation; i++) {
        for (int j = num - 1; j > 0; --j) {
            k = rand() % (j + 1);
            temp = init_perm[k];
            init_perm[k] = init_perm[j];
            init_perm[j] = temp;
        }
        memcpy(i_perm_matrix[i], init_perm, num * sizeof(int));
    }
#endif
}

/*
 * permute group index: i_perm
 */
void resample_indicator_label(int *i_perm, int *i_perm_tmp, int n, int *n1) {
    int i, j, temp, tmp0 = 0, tmp1 = 0;
#ifdef R_BUILD
    for (i = n - 1; i > 0; --i) {
        j = ((int) round(RAND_MAX * unif_rand())) % (i + 1);
        temp = i_perm[j];
        i_perm[j] = i_perm[i];
        i_perm[i] = temp;
    }
#else
    srand((unsigned) time(NULL));
    for (i = n - 1; i > 0; --i) {
        j = rand() % (i + 1);
        temp = i_perm[j];
        i_perm[j] = i_perm[i];
        i_perm[i] = temp;
    }
#endif
    for (i = 0; i < n; i++) {
        if (i_perm[i] == 1) {
            i_perm_tmp[tmp0++] = i;
        } else {
            i_perm_tmp[*n1 + tmp1] = i;
            tmp1++;
        }
    }
}

void resample_indicator_label_matrix(int **i_perm_matrix, int **i_perm_tmp_matrix,
                                     int *init_perm, int *init_perm_tmp, int num_permutation, int n, int *n1) {
    int i, j, temp, tmp0, tmp1;
#ifdef R_BUILD
    GetRNGstate();
    for (int k = 0; k < num_permutation; ++k) {
        for (i = n - 1; i > 0; --i) {
            j = ((int) round(RAND_MAX * unif_rand())) % (i + 1);
            temp = init_perm[j];
            init_perm[j] = init_perm[i];
            init_perm[i] = temp;
        }
        memcpy(i_perm_matrix[k], init_perm, n * sizeof(int));
        tmp0 = 0;
        tmp1 = 0;
        for (i = 0; i < n; i++) {
            if (init_perm[i] == 1) {
                init_perm_tmp[tmp0++] = i;
            } else {
                init_perm_tmp[*n1 + tmp1] = i;
                tmp1++;
            }
        }
        memcpy(i_perm_tmp_matrix[k], init_perm_tmp, n * sizeof(int));
    }
    PutRNGstate();
#else
    srand((unsigned) time(NULL));
    for (int k = 0; k < num_permutation; ++k) {
        for (i = n - 1; i > 0; --i) {
            j = rand() % (i + 1);
            temp = init_perm[j];
            init_perm[j] = init_perm[i];
            init_perm[i] = temp;
        }
        memcpy(i_perm_matrix[k], init_perm, n * sizeof(int));
        tmp0 = 0;
        tmp1 = 0;
        for (i = 0; i < n; i++) {
            if (init_perm[i] == 1) {
                init_perm_tmp[tmp0++] = i;
            } else {
                init_perm_tmp[*n1 + tmp1] = i;
                tmp1++;
            }
        }
        memcpy(i_perm_tmp_matrix[k], init_perm_tmp, n * sizeof(int));
    }
#endif
}

int random_index_thread_wrap(int i) {
    return random_index_thread(i);
}

void resample3_thread(int *permuted_arr, int *i_perm, int *i_perm_tmp, int n, int *n1) {
    int i, j, temp, tmp0, tmp1;
    for (i = n - 1; i > 0; --i) {
        j = permuted_arr[i];
        temp = i_perm[j];
        i_perm[j] = i_perm[i];
        i_perm[i] = temp;
    }

    tmp0 = 0;
    tmp1 = 0;
    for (i = 0; i < n; i++) {
        if (i_perm[i] == 1) {
            i_perm_tmp[tmp0++] = i;
        } else {
            i_perm_tmp[*n1 + tmp1] = i;
            tmp1++;
        }
    }
}

/* Arrange the N elements of ARRAY in random order.
 Only effective if N is much smaller than RAND_MAX;
 if this may not be the case, use a better random
 number generator. */
void shuffle(int *array, int *N) {
    int j, tmp;
#ifdef R_BUILD
    GetRNGstate();
    for (int i = *N - 1; i > 0; --i) {
        j = ((int) round(RAND_MAX * unif_rand())) % (i + 1);
        tmp = array[j];
        array[j] = array[i];
        array[i] = tmp;
    }
    PutRNGstate();
#else
    srand((unsigned) time(NULL));
    for (int i = 0; i < *N - 1; i++) {
        j = random_index(*N, i);
        tmp = array[j];
        array[j] = array[i];
        array[i] = tmp;
    }
#endif
}

void shuffle_value(double *array, int *N) {
    int j;
    double tmp;
#ifdef R_BUILD
    GetRNGstate();
    for (int i = *N - 1; i > 0; --i) {
        j = ((int) round(RAND_MAX * unif_rand())) % (i + 1);
        tmp = array[j];
        array[j] = array[i];
        array[i] = tmp;
    }
    PutRNGstate();
#else
    srand((unsigned) time(NULL));
    for (int i = 0; i < *N - 1; i++) {
        j = random_index(*N, i);
        tmp = array[j];
        array[j] = array[i];
        array[i] = tmp;
    }
#endif
}

void shuffle_value_matrix(double **value_matrix, double *init_value, int num_permutation, int num) {
    int k;
    double tmp;
#ifdef R_BUILD
    GetRNGstate();
    for (int i = 0; i < num_permutation; i++) {
        for (int j = num - 1; j > 0; --j) {
            k = ((int) round(RAND_MAX * unif_rand())) % (j + 1);
            tmp = init_value[k];
            init_value[k] = init_value[j];
            init_value[j] = tmp;
        }
        memcpy(value_matrix[i], init_value, num * sizeof(double));
    }
    PutRNGstate();
#else
    int i, j;
    srand((unsigned) time(NULL));
    for (k = 0; k < num_permutation; ++k) {
        for (i = 0; i < num - 1; i++) {
            j = random_index(num, i);
            tmp = init_value[j];
            init_value[j] = init_value[i];
            init_value[i] = tmp;
        }
        memcpy(value_matrix[i], init_value, num * sizeof(double));
    }
#endif
}

int pending_interrupt() {
    int interrupt_status = 0;
    interrupt_status = pending_interrupt_status();
    return interrupt_status;
}

void print_stop_message() {
    print_stop_message_internal();
}