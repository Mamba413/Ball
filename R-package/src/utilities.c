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

#include "utilities.h"
#include "string.h"
#include "stdio.h"
#include "stdlib.h"

#ifdef R_BUILD

#include "Rmath.h"
#include "R.h"
#include "Rinternals.h"

#define RAND_MAX_CONSTANT 2147483647

#else

#include "math.h"
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

/**
 * Coupute the cumulative sample size vector
 * @param size: an array contain sample size in each group
 * @param cumulate_size: the cumulative sums to be computed
 * @param k: group number
 * @example : Input : size = [15, 25, 10], k = 3
 * @example : Outpur : cumsum_size = [0, 15, 40]
 */
void compute_cumsum_size(int *cumsum_size, int *size, int *k) {
    int i;
    for (i = 0; i < (*k); i++) {
        if (i == 0) {
            cumsum_size[i] = 0;
        } else {
            cumsum_size[i] = cumsum_size[i - 1] + size[i - 1];
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

/**
 * quick sort array with size n
*/
void quick_sort(double *array, int n) {
    quick_sort_recursive(array, 0, n - 1);
}

double compute_pvalue(double ball_stat_value, double *permuted_stat, int R) {
    double larger_num = 0.0;
    for (int i = 0; i < R; i++) {
        if (permuted_stat[i] >= ball_stat_value) {
            larger_num += 1.0;
        }
    }
    double p_value = (1.0 + larger_num) / (1.0 + R);
    return p_value;
}

void compute_batch_pvalue(double *ball_stat, const double *permuted_stat, double *p_value,
                          int batch_num, int R) {
    int all_stat_num = R + batch_num;
    int *rank_all = (int *) calloc((size_t) all_stat_num, sizeof(int));
    int *rank_batch = (int *) calloc((size_t) batch_num, sizeof(int));
    memset(rank_all, all_stat_num, sizeof(int));
    memset(rank_batch, batch_num, sizeof(int));
    double *all_stat = (double *) calloc((size_t) all_stat_num, sizeof(double));
    for (int i = 0; i < batch_num; ++i) {
        all_stat[i] = ball_stat[i];
    }
    for (int i = 0; i < R; ++i) {
        all_stat[batch_num + i] = permuted_stat[i];
    }
    quick_rank_min(all_stat, rank_all, all_stat_num);
    quick_rank_min(ball_stat, rank_batch, batch_num);
    double all_prop = 1.0 / (1.0 + R);
    for (int i = 0; i < batch_num; i++) {
        p_value[i] = (1.0 + (R + batch_num - rank_all[i]) - (batch_num - rank_batch[i])) * all_prop;
    }
}

void merge(double *vector, int *index, int *number, int start, int mid, int end) {
    const int left_size = mid - start + 1, right_size = end - mid;
    double left[left_size], right[right_size];
    int left_index[left_size], right_index[right_size];
    int left_merged = 0, right_merged = 0, total_merged = 0;
    for (int i = start; i <= mid; ++i) {
        left[i - start] = vector[i];
        left_index[i - start] = index[i];
    }
    for (int i = mid + 1; i <= end; ++i) {
        right[i - mid - 1] = vector[i];
        right_index[i - mid - 1] = index[i];
    }
    while (left_merged < left_size && right_merged < right_size) {
        if (left[left_merged] < right[right_merged]) {
            number[left_index[left_merged]] += right_merged;
            vector[start + total_merged] = left[left_merged];
            index[start + total_merged] = left_index[left_merged];
            ++left_merged;
            ++total_merged;
        } else {
            vector[start + total_merged] = right[right_merged];
            index[start + total_merged] = right_index[right_merged];
            ++right_merged;
            ++total_merged;
        }
    }
    while (left_merged < left_size) {
        number[left_index[left_merged]] += right_merged;
        vector[start + total_merged] = left[left_merged];
        index[start + total_merged] = left_index[left_merged];
        ++left_merged;
        ++total_merged;
    }
    while (right_merged < right_size) {
        vector[start + total_merged] = right[right_merged];
        index[start + total_merged] = right_index[right_merged];
        ++right_merged;
        ++total_merged;
    }
}

void merge_sort(double *vector, int *index, int *number, int start, int end) {
    if (end - start < 1) return;
    int mid = (start + end) >> 1;
    merge_sort(vector, index, number, start, mid);
    merge_sort(vector, index, number, mid + 1, end);
    merge(vector, index, number, start, mid, end);
}

void count_smaller_number_after_self_solution(double *vector, int *number, const int num) {
    int index[num];
    for (int i = 0; i < num; ++i) {
        index[i] = i;
    }
    merge_sort(vector, index, number, 0, num - 1);
}

void count_smaller_number_after_self_solution2(double *vector, int *index, int *number, const int num) {
    merge_sort(vector, index, number, 0, num - 1);
}

void Merge(int *permutation, int *source, int *inversion_count, int dim, int n) {
    int *left = (int *) malloc(n * sizeof(int));
    int *right = (int *) malloc(n * sizeof(int));
    int *left_source = (int *) malloc(n * sizeof(int));
    int *right_source = (int *) malloc(n * sizeof(int));
    int left_index = 0, right_index = 0;
    int i, half_dim = dim >> 1, nleft = half_dim, nright = dim - half_dim;
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
            if (left[left_index] <= right[right_index]) {
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

void Inversions(int *permutation, int *source, int *inversion_count, int dim, int n) {
    if (dim == 1) {
        return;
    } else {
        Inversions(permutation, source, inversion_count, dim >> 1, n);
        Inversions(&permutation[dim >> 1], &source[dim >> 1], inversion_count, dim - (dim >> 1), n);
        Merge(permutation, source, inversion_count, dim, n);
    }
}

/**
 * Give an univariate array, compute the lower and upper boundary of each ball.
 * @param n : size of array
 * @param zidx : zidx[i] = j means the i-smallest element of z is z[j]
 * @param z : an value array
 * @param lowzidx : the lower boundary to be computed
 * @param higzidx : the upper boundary to be computed
 *
 * @details:
 * The algorithm like the quick-sort. It finds the upper and lower in turn.
 * For each i, the computation complexity is O(n); Thus, this function take O(n^2) times.
 * Algorithm detail:
 * For i=1, ..., n, carry out:
 * init the low and upper index: l=1, u=N
 * while (l <= u):
 * if (z[u] - z[i]) >= (z[i] - z[l]) ==> for pair (i, u), upper index: u, lower index: l ==> update upper index: u = u - 1
 * if (z[u] - z[i]) < (z[i] - z[l]) ==> for pair (i, l), upper index: u, lower index: l ==> update lower index: l = l - 1
 *
 * @refitem https://arxiv.org/pdf/1811.03750.pdf, Section 3.3, Algorithm 2
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

/*
This function try to derive the rank of Rx belong to group 1.
tmp, lastpos, lastval are introduced to resolve the situtation when ties appear
*/
void Findx2(int *Rxy, int *Ixy, int *i_perm, int *n1, int *n2, int *Rx) {
    int j, lastpos, lastval, n, tmp;
    n = *n1 + *n2;

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

/**
 * The rank computation for univariate Ball Covariance
 * @refitem: http://www.jmlr.org/papers/volume17/14-441/14-441.pdf [section 3.1.2]
*/
void computeRank(int n, int **Rank) {
    int i, j;
    for (i = 1; i < n; i++) {
        for (j = 1; j < n; j++) {
            Rank[i][j] += (Rank[i][j - 1] + Rank[i - 1][j] - Rank[i - 1][j - 1]);
        }
    }
}

/**
 * The rank computation for univariate Ball Covariance
 * @refitem: http://www.jmlr.org/papers/volume17/14-441/14-441.pdf [section 3.1.2]
*/
void initRank(int n, int **Rank, int *xrank, int *yrank, int *i_perm) {
    int i, j;
    for (i = 0; i < n + 1; i++) {
        for (j = 0; j < n + 1; j++) {
            Rank[i][j] = 0;
        }
    }
    for (i = 0; i < n; i++) {
        Rank[xrank[i] + 1][yrank[i_perm[i]] + 1] += 1;
    }
    computeRank(n + 1, Rank);
}

void initRank_bcor(int n, int **Rank, int *xrank, int *yrank) {
    int i, j;
    for (i = 0; i < n + 1; i++) {
        for (j = 0; j < n + 1; j++) {
            Rank[i][j] = 0;
        }
    }
    for (i = 0; i < n; i++) {
        Rank[xrank[i] + 1][yrank[i] + 1] += 1;
    }
    computeRank(n + 1, Rank);
}

void quick_rank_max_with_index(const double *x, const int *x_index, int *r, int n) {
    int i_loc, rank_value = n, tmp = 1;
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
}

/**
 * Rank a value vector in a max manner
 * @param x value vector
 * @param r save the rank result
 * @param n the size of vector
 * @example
 * x = {1.2, 1.3, 1.3, 1.3, 1.3, 0.9};
 * r ={0, 0, 0, 0, 0, 0};
 * quick_rank_max(x, r, 6);
 * r = {2, 6, 6, 6, 6, 1};
 */
void quick_rank_max(const double *x, int *r, int n) {
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

/**
 * Rank a value vector in a max manner
 * @param x value vector
 * @param r save the rank result
 * @param n the size of vector
 * @example
 * x = {1.2, 1.3, 1.3, 1.3, 1.3, 0.9};
 * r ={0, 0, 0, 0, 0, 0};
 * quick_rank_min(x, r, 6);
 * r = {2, 3, 3, 3, 3, 1}
 */
void quick_rank_min(const double *x, int *r, int n) {
    int *x_index, rank_value = 1, i_loc, tmp = 1;
    double *x_cpy;
    x_index = (int *) malloc(n * sizeof(int));
    x_cpy = (double *) malloc(n * sizeof(double));
    for (int j = 0; j < n; j++) {
        x_index[j] = j;
        x_cpy[j] = x[j];
    }
    quicksort(x_cpy, x_index, 0, n - 1);
    r[x_index[0]] = 1;
    for (int i = 1; i < n; i++) {
        i_loc = x_index[i];
        if (x[i_loc] == x[x_index[i - 1]]) {
            r[i_loc] = rank_value;
            tmp++;
        } else {
            rank_value = rank_value + tmp;
            r[i_loc] = rank_value;
            tmp = 1;
        }
    }
    free(x_index);
    free(x_cpy);
}

/**
 * compute the rank
 * @param n : the size of zrank
 * @param zrank : a empty array to be computed
 * @param z :
 * @param zidx : the same as the zidx in createidx function
 * @example : Input: n=6, zrank = []; z = [1, 2, 3, 4, 5, 5], zidx = [3, 1, 5, 2, 6, 4];
 * @example : Output: zrank = [2, 4, 1, 5, 3, 5];
*/
void ranksort(int *n, int *zrank, double *z, int *zidx) {
    int i, last_position = 0;
    double last_value = -1.0;

    for (i = *n - 1; i >= 0; i--) {
        if (last_value != z[i]) {
            last_position = i;
        }
        last_value = z[i];
        zrank[zidx[i]] = last_position;
    }
}

/**
 * A matrix version for ranksort
 * @param n the size of
 * @param Rxy a empty 2D array to be computed
 * @param Dxy a 2D value array
 * @param Ixy : Ixy[i][j] = q means the j-smallest value of the Dxy[i] is Dxy[i][j]
 */
void ranksort2(int n, int **Rxy, double **Dxy, int **Ixy) {
    int i, j, lastpos = n - 1;
    double lastval;
    for (i = 0; i < n; i++) {
        lastval = -1;
        for (j = n - 1; j >= 0; j--) {
            if (lastval != Dxy[i][j]) {
                lastpos = j;
            }
            lastval = Dxy[i][j];
            Rxy[i][Ixy[i][j]] = lastpos;
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
            quick_rank_max(x_part, x_part_rank, n);
            for (j = 0; j < n; j++) {
                Rx[i][j][h] = x_part_rank[j];
            }
        }
    }
    free(x_part);
    free(x_part_rank);
}

/**
 * Apply quicksort algorithm to a, and the index exchange result is recorded in idx.
 * @param array an array to be sort in ascending order
 * @param idx an order index array. idx[i] = j means the i smallest value of array is a[j].
*/
void quicksort(double *array, int *idx, int l, int u) {
    int i, m, idx_temp;
    double a_temp;
    if (l >= u)
        return;
    m = l;
    for (i = l + 1; i <= u; i++) {
        if (array[i] < array[l]) {
            ++m;
            idx_temp = idx[m];
            idx[m] = idx[i];
            idx[i] = idx_temp;

            a_temp = array[m];
            array[m] = array[i];
            array[i] = a_temp;
        }
    }
    idx_temp = idx[l];
    idx[l] = idx[m];
    idx[m] = idx_temp;

    a_temp = array[l];
    array[l] = array[m];
    array[m] = a_temp;

    quicksort(array, idx, l, m - 1);
    quicksort(array, idx, m + 1, u);
}

/**
 * @inherit quicksort
*/
void quicksort_int(int *array, int *idx, int l, int u) {
    int i, m, idx_temp;
    int a_temp;
    if (l >= u)
        return;
    m = l;
    for (i = l + 1; i <= u; i++) {
        if (array[i] < array[l]) {
            ++m;
            idx_temp = idx[m];
            idx[m] = idx[i];
            idx[i] = idx_temp;

            a_temp = array[m];
            array[m] = array[i];
            array[i] = a_temp;
        }
    }
    idx_temp = idx[l];
    idx[l] = idx[m];
    idx[m] = idx_temp;

    a_temp = array[l];
    array[l] = array[m];
    array[m] = a_temp;

    quicksort_int(array, idx, l, m - 1);
    quicksort_int(array, idx, m + 1, u);
}

/**
 * Apply quicksort algorithm to a and b. After sorted, a is increasing while b is decreasing.
 * idx recode the index exchange result
 * @example : Input a = [2, 2, 1, 3], b = [1, 2, 3, 4], idx = [1, 2, 3, 4];
 * @example : Output a = [1, 2, 2, 3], b = [3, 1, 2, 4], idx = [3, 1, 2, 4]
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

void quicksort3(double *a, double *b, int *idx, int l, int u) {
    int i, m, idx_temp;
    double a_temp;
    if (l >= u)
        return;
    m = l;
    for (i = l + 1; i <= u; i++) {
        if (a[i] < a[l] || (a[i] == a[l] && b[i] < b[l])) {
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

    quicksort3(a, b, idx, l, m - 1);
    quicksort3(a, b, idx, m + 1, u);
}

double **alloc_matrix(int r, int c) {
    /* allocate a matrix with r rows and c columns */
    int i;
    double **matrix;
    matrix = (double **) calloc((size_t) r, sizeof(double *));
    for (i = 0; i < r; i++)
        matrix[i] = (double *) calloc((size_t) c, sizeof(double));
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
    matrix = (int **) calloc((size_t) r, sizeof(int *));
    for (i = 0; i < r; i++)
        matrix[i] = (int *) calloc((size_t) c, sizeof(int));
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

int ***alloc_int_square_matrix_list(int *size, int number) {
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

void free_int_square_matrix_list(int ***arr3d, int *size, int num) {
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

void Category_distance(const double *x, double **Dx, int n) {
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            Dx[i][j] = x[i] == x[j] ? 0 : 1;
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

/**
 * shuffle i_perm and preserve the inversion in i_perm
 * @param i_perm : array to be shuffled
 * @param i_perm : array convert i_perm to {0, 1, ..., n - 1}
 * @param n : size of i_perm
 *
*/
void resample(int *i_perm, int *i_perm_inv, int *n) {
    int i, j, temp;
#ifdef R_BUILD
    GetRNGstate();
    for (i = *n - 1; i > 0; --i) {
        j = ((int) fround(RAND_MAX_CONSTANT * unif_rand(), 0)) % (i + 1);
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
            k = ((int) fround(RAND_MAX_CONSTANT * unif_rand(), 0)) % (j + 1);
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
#ifdef R_BUILD
    // GetRNGstate();
    int i, j, k, temp;
    // Rprintf("RAND_MAX_CONSTANT: %d", RAND_MAX_CONSTANT);
    for (k = 0; k < *r; k++) {
        for (i = *c - 1; i > 0; --i) {
            j = ((int) fround(RAND_MAX_CONSTANT * unif_rand(), 0)) % (i + 1);
            temp = i_perm[k][j];
            i_perm[k][j] = i_perm[k][i];
            i_perm[k][i] = temp;
        }
    }
    // PutRNGstate();
#else
    srand((unsigned) time(NULL));
    int i, j, k, temp;
    for (k = 0; k < *r; k++) {
        for (i = *c - 1; i > 0; --i) {
            j = (rand()) % (i + 1);
            temp = i_perm[k][j];
            i_perm[k][j] = i_perm[k][i];
            i_perm[k][i] = temp;
        }
    }
#endif
}

void resample_matrix_3d(int ***i_perm, int **init_perm, int *h, int *r, int *c) {
    int j, temp;
#ifdef R_BUILD
    // GetRNGstate();
    for (int l = 0; l < *h; ++l) {
        for (int k = 0; k < *r; k++) {
            for (int i = *c - 1; i > 0; --i) {
                j = ((int) fround(RAND_MAX_CONSTANT * unif_rand(), 0)) % (i + 1);
                temp = init_perm[k][j];
                init_perm[k][j] = init_perm[k][i];
                init_perm[k][i] = temp;
            }
        }
        for (int k = 0; k < *r; ++k) {
            memcpy(i_perm[l][k], init_perm[k], (*c * sizeof(int)));
        }
    }
    // PutRNGstate();
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
        j = ((int) fround(RAND_MAX_CONSTANT * unif_rand(), 0)) % (i + 1);
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
            j = ((int) fround(RAND_MAX_CONSTANT * unif_rand(), 0)) % (i + 1);
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

/*
 * permute group index: i_perm
 */
void resample_indicator_label(int *i_perm, int *i_perm_tmp, int n, int *n1) {
    int i, j, temp, tmp0 = 0, tmp1 = 0;
#ifdef R_BUILD
    GetRNGstate();
    for (i = n - 1; i > 0; --i) {
        j = ((int) fround(RAND_MAX_CONSTANT * unif_rand(), 0)) % (i + 1);
        temp = i_perm[j];
        i_perm[j] = i_perm[i];
        i_perm[i] = temp;
    }
    PutRNGstate();
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
            j = ((int) fround(RAND_MAX_CONSTANT * unif_rand(), 0)) % (i + 1);
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

/* Arrange the N elements of ARRAY in random order.
 Only effective if N is much smaller than RAND_MAX_CONSTANT;
 if this may not be the case, use a better random
 number generator. */
void shuffle(int *array, int *N) {
    int j, tmp;
#ifdef R_BUILD
    GetRNGstate();
    for (int i = *N - 1; i > 0; --i) {
        j = ((int) fround(RAND_MAX_CONSTANT * unif_rand(), 0)) % (i + 1);
        tmp = array[j];
        array[j] = array[i];
        array[i] = tmp;
    }
    PutRNGstate();
#else
    srand((unsigned) time(NULL));
    for (int i = 0; i < *N - 1; i++) {
        j = i + rand() / (RAND_MAX_CONSTANT / (*N - i) + 1);
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
        j = ((int) fround(RAND_MAX_CONSTANT * unif_rand(), 0)) % (i + 1);
        tmp = array[j];
        array[j] = array[i];
        array[i] = tmp;
    }
    PutRNGstate();
#else
    srand((unsigned) time(NULL));
    for (int i = 0; i < *N - 1; i++) {
        j = i + rand() / (RAND_MAX_CONSTANT / (*N - i) + 1);
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
            k = ((int) fround(RAND_MAX_CONSTANT * unif_rand(), 0)) % (j + 1);
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
            j = i + rand() / (RAND_MAX_CONSTANT / (num - i) + 1);
            tmp = init_value[j];
            init_value[j] = init_value[i];
            init_value[i] = tmp;
        }
        memcpy(value_matrix[i], init_value, num * sizeof(double));
    }
#endif
}

#ifdef R_BUILD
void R_CheckUserInterrupt_fn(void *dummy) {
  R_CheckUserInterrupt();
}
#endif

int pending_interrupt(void) {
    int interrupt_status = 0;
#ifdef R_BUILD
    interrupt_status = !(R_ToplevelExec(R_CheckUserInterrupt_fn, NULL));
#endif
    return interrupt_status;
}

void print_stop_message(void) {
#ifdef R_BUILD
    Rprintf("Process stop due to user interruption! \n");
#else
    printf("Process stop due to user interruption! \n");
#endif

}

void declare_gwas_screening(void) {
#ifdef R_BUILD
    Rprintf("=========== Pre-screening SNPs ===========\n");
#else
    printf("=========== Pre-screening SNPs ===========\n");
#endif
}

void declare_gwas_refining(int i, int refine_num) {
#ifdef R_BUILD
    if (i == -1) {
        Rprintf("=========== Refining the p-value of %d SNP ===========\n", refine_num);
    } else {
        Rprintf("Refining SNP... Progress: %d/%d. \n", i, refine_num);
    }
#else
    if (i == -1) {
        printf("=========== Refining the p-value of %d SNP ===========\n", refine_num);
    } else {
        printf("Refining SNP... Progress: %d/%d. \n", i, refine_num);
    }
#endif
}

void print_pvalue(double pvalue) {
#ifdef R_BUILD
    Rprintf("Refined p-value: %11.10f, ", pvalue);
#else
    printf("Refined p-value: %11.10f, ", pvalue);
#endif
}

void print_cost_time(int second) {
    char result[200] = "";
    beautify_time(result, second);
#ifdef R_BUILD
    Rprintf("cost time: %s.\n", result);
#else
    printf("cost time: %s.\n", result);
#endif
}

void int_to_string(char str[], int number) {
    // sprintf(str, "%d", number);
    snprintf(str, 3, "%d", number%100u);
}

void beautify_time(char result[], int seconds) {
    // Add seconds, minutes, hours, days if larger than zero
    char second_str[] = " seconds";
    char out_second_str[50] = "";
    int out_seconds = seconds % 60;
    int_to_string(out_second_str, out_seconds);
    strcat(out_second_str, second_str);
    // minute
    int out_minutes = (seconds / 60) % 60;
    char out_minutes_str[100] = "";
    if (seconds / 60 == 0) {
        strcpy(result, out_second_str);
        return;
    } else if (out_minutes == 1) {
        strcpy(out_minutes_str, "1 minute, ");
        strcat(out_minutes_str, out_second_str);
    } else {
        char minutes[] = " minutes, ";
        int_to_string(out_minutes_str, out_minutes);
        strcat(out_minutes_str, minutes);
        strcat(out_minutes_str, out_second_str);
    }
    // hour
    int out_hours = (seconds / 3600) % 24;
    char out_hours_str[150] = "";
    if (seconds / 3600 == 0) {
        strcpy(result, out_minutes_str);
        return;
    } else if (out_hours == 1) {
        strcpy(out_hours_str, "1 hour, ");
        strcat(out_hours_str, out_minutes_str);
    } else {
        char hours[] = " hours, ";
        int_to_string(out_hours_str, out_hours);
        strcat(out_hours_str, hours);
        strcat(out_hours_str, out_minutes_str);
    }
    // day
    int out_days = (seconds / 86400);
    char out_days_str[200] = "";
    if (out_days == 0) {
        strcpy(result, out_hours_str);
        return;
    } else if (out_days == 1) {
        strcpy(out_days_str, "1 day, ");
        strcat(out_days_str, out_hours_str);
    } else {
        char days[] = " days, ";
        int_to_string(out_days_str, out_days);
        strcat(out_days_str, days);
        strcat(out_days_str, out_hours_str);
    }
    strcpy(result, out_days_str);
}
