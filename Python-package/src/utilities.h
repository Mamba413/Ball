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
void count_smaller_number_after_self_solution(double *vector, int *number, int num);
void count_smaller_number_after_self_solution2(double *vector, int *index, int *number, int num);
void find_group_relative_location(int *group_relative_location, int *group, int *cumsum_size, int num, int K);
void compute_cumsum_size(int *cumsum_size, int *size, int *k);
void quick_sort_recursive(double *arr, int start, int end);
void quick_sort(double *array, int n);
double compute_pvalue(double ball_stat_value, double *permuted_stat, int R);
void compute_batch_pvalue(double *ball_stat, const double *permuted_stat, double *p_value, int batch_num, int R);
void Merge(int *permutation, int *source, int *inversion_count, int dim, int n);
void Inversions(int *permutation, int *source, int *inversion_count, int dim, int n);
void sort(int *n, int *zidx, double *z, int **dzidx);
void createidx(int *n, int *zidx, double *z, int **lowzidx, int **higzidx);
void ranksort2(int n, int **Rxy, double **Dxy, int **Ixy);
void Findx2(int *Rxy, int *Ixy, int *i_perm, int *n1, int *n2, int *Rx);
void Findx(int **Rxy, int **Ixy, int *i_perm, int *n1, int *n2, int **Rx);
void ranksort3(int n, int *xyidx, double *xy, int **Rxy, int **Ixy);
void computeRank(int n, int **Rank);
void initRank(int n, int **Rank, int *xrank, int *yrank, int *i_perm);
void initRank_bcor(int n, int **Rank, int *xrank, int *yrank);
void quicksort(double *array, int *idx, int l, int u);
void quicksort_int(int *array, int *idx, int l, int u);
void quicksort2(double *a, double *b, int *idx, int l, int u);
void quicksort3(double *a, double *b, int *idx, int l, int u);
void rank_matrix_3d(double ***Dx, int n, int k, int ***Rx);
void quick_rank_max_with_index(const double *x, const int *x_index, int *r, int n);
void quick_rank_max(const double *x, int *r, int n);
void quick_rank_min(const double *x, int *r, int n);
void quick_rank_max_index(const double *x, int *x_index, int *r, int n);
double **alloc_matrix(int r, int c);
double ***alloc_3d_matrix(int r, int c, int h);
int **alloc_int_matrix(int r, int c);
int ***alloc_3d_int_matrix(int r, int c, int h);
int ***alloc_int_square_matrix_list(int* size, int number);
void free_matrix(double **matrix, int r, int c);
void free_3d_matrix(double ***arr3D, int r, int c);
void free_int_matrix(int **matrix, int r, int c);
void free_3d_int_matrix(int ***arr3D, int r, int c);
void free_int_square_matrix_list(int ***arr3d, int* size, int num);
void vector2matrix(double *x, double **y, int N, int d, int isroworder);
void distance2matrix(double *distance, double **distance_matrix, int n);
void vector2matrix3d(double *x, double ***y, int r, int c, int h, int isroworder);
void Euclidean_distance(double *x, double **Dx, int n, int d);
void Category_distance(const double *x, double **Dx, int n);
void distance(double *x, double *Dx, int *n, int *d);
void shuffle(int *array, int *N);
void shuffle_value(double *array, int *N);
int pending_interrupt();
void print_stop_message();
void resample(int *i_perm, int *i_perm_inv, int *n);
void shuffle_indicator_inv_matrix(int **i_perm_matrix, int **i_perm_matrix_inv, int *init_perm, int *init_perm_inv,
                                  int num_permutation, int num);
void resample_matrix(int **i_perm, int *r, int *c);
void resample2(int *i_perm, int *n);
void resample2_matrix(int **i_perm, int *init_perm, int num_permutation, int n);
void resample_matrix_3d(int ***i_perm, int **init_perm, int *h, int *r, int *c);
void resample_indicator_label(int *i_perm, int *i_perm_tmp, int n, int *n1);
void resample_indicator_label_matrix(int **i_perm_matrix, int **i_perm_tmp_matrix,
                                     int *init_perm, int *init_perm_tmp, int num_permutation, int n, int *n1);
void shuffle_value_matrix(double **value_matrix, double *init_value, int num_permutation, int num);
void resample3(int *i_perm, int *i_perm_tmp, int n, int *n1);
void ranksort(int *n, int *zrank, double *z, int *zidx);
void distance2matrix3d(double *distance, double ***distance_matrix3d, int n, int v);
void declare_gwas_screening();
void declare_gwas_refining(int i, int refine_num);
void print_pvalue(double pvalue);
void print_cost_time(int second);
void beautify_time(char result[], int seconds);

#endif /* UTILITIES_H_ */