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
#ifndef UTILITIES_H_
#define UTILITIES_H_

#define MAX(a,b) (((a)>(b))?(a):(b));
void swap(double *x, double *y);
void quick_sort_recursive(double *arr, int start, int end);
void quick_sort(double *arr, int len);
double compute_pvalue(double ball_stat_value, double *permuted_stat, int R);
void Merge(int *permutation, int *source, int *inversion_count, int dim, int n);
int Inversions(int *permutation, int *source, int *inversion_count, int dim, int n);
void sort(int *n, int *zidx, double *z, int **dzidx);
void createidx(int *n, int *zidx, double *z, int **lowzidx, int **higzidx);
void ranksort2(int n, int **Rxy, double **Dxy, int **Ixy);
void Findx2(int *Rxy, int *Ixy, int *i_perm, int *n1, int *n2, int *Rx);
void Findx(int **Rxy, int **Ixy, int *i_perm, int *n1, int *n2, int **Rx);
void ranksort3(int n, int *xyidx, double *xy, int **Rxy, int **Ixy);
void computeRank(int n, int **Rank);
void initRank(int n, int **Rank, int *xrank, int *yrank, int *i_perm);
void initRank_bcor(int n, int **Rank, int *xrank, int *yrank);
void quicksort(double *a, int *idx, int l, int u);
void quicksort2(double *a, double *b, int *idx, int l, int u);
double **alloc_matrix(int r, int c);
double ***alloc_3d_matrix(int r, int c, int h);
int **alloc_int_matrix(int r, int c);
void free_matrix(double **matrix, int r, int c);
void free_3d_matrix(double ***arr3D, int r, int c);
void free_int_matrix(int **matrix, int r, int c);
void vector2matrix(double *x, double **y, int N, int d, int isroworder);
void vector_2_3dmatrix(double *x, double ***y, int r, int c, int h, int isroworder);
void Euclidean_distance(double *x, double **Dx, int n, int d);
void distance(double *x, double *Dx, int *n, int *d);
void shuffle(int *array, int *N);
void shuffle_value(double *array, int *N);
int pending_interrupt();
void print_stop_message();
void resample(int *i_perm, int *i_perm_inv, int *n);
void resample_matrix(int **i_perm, int *r, int *c);
void resample2(int *i_perm, int *n);
void resample3(int *i_perm, int *i_perm_tmp, int n, int *n1);
int random_index_thread_wrap(int i);
void resample3_thread(int *permuted_arr, int *i_perm, int *i_perm_tmp, int n, int *n1);
void ranksort(int *n, int *zrank, double *z, int *zidx);


#endif /* UTILITIES_H_ */