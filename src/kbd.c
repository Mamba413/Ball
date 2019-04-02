//
// Created by JinZhu on 2019/3/18.
//

#include "string.h"
#include "stdlib.h"
#include "stdio.h"
#include "utilities.h"

void ball_divergence2(double *bd_stat, int **full_rank, int **sub_rank1, int **sub_rank2, int n1, int n2) {
    double pxx, pxy, pyx, pyy, diff;
    double A_nm = 0.0, C_nm = 0.0;
    double n1_prop = 1.0 / n1, n2_prop = 1.0 / n2;
    for (int i = 0; i < n1; ++i) {
        for (int j = 0; j < n1; ++j) {
            pxx = sub_rank1[i][j] * n1_prop;
            pxy = (full_rank[i][j] - sub_rank1[i][j]) * n2_prop;
            diff = pxx - pxy;
            A_nm += diff * diff;
            // TODO: weighted Ball Divergence
        }
    }
    A_nm *= (n1_prop * n1_prop);
    for (int i = 0; i < n2; ++i) {
        for (int j = 0; j < n2; ++j) {
            pyy = sub_rank2[i][j] * n2_prop;
            pyx = (full_rank[i + n1][j + n1] - sub_rank2[i][j]) * n1_prop;
            diff = pyy - pyx;
            C_nm += diff * diff;
            // TODO: weighted Ball Divergence
        }
    }
    C_nm *= (n2_prop * n2_prop);
    bd_stat[0] = A_nm + C_nm;
    bd_stat[1] = A_nm + C_nm;
}

/**
 *
 * @param sub_rank
 * @param index_matrix
 * @param cumsum_size : [0, 4, 9]
 * @param label : [0, 0, 0, 0, 1, 1, 1, 1, 1, 2, 2, 2, 2]
 * @param num : 13
 * @param max_k : 2
 */
void sub_rank_finder(int ***sub_rank, double **distance_matrix, int **index_matrix, const int *label,
                     const int *group_relative_location,const int *cumsum_size, int num, int max_k) {
    /*
     * g_index: indicator for group
     * s_index: indicator for order
     */
    int s_index, g_index;
    int *init_sub_rank;
    init_sub_rank = (int *) malloc((max_k + 1) * sizeof(int));
    for (int i = 0; i < num; ++i) {
        for (int k = 0; k < (max_k + 1); ++k) {
            init_sub_rank[k] = 1;
        }
        for (int j = 0; j < num; ++j) {
            s_index = index_matrix[i][j];
            g_index = label[s_index];
            if (label[i] == g_index) {
                sub_rank[g_index][group_relative_location[i] - cumsum_size[g_index]][group_relative_location[s_index] -
                                                                                     cumsum_size[g_index]] = init_sub_rank[g_index];
                init_sub_rank[g_index] += 1;
            }
        }
    }
}

/**
 * @inherit sub_rank_finder
 * @param size : [4, 5, 4]
 */
void full_rank_finder(int ***full_rank, double **distance_matrix, int **index_matrix, int *label,
                      int *group_relative_location, int *cumsum_size, int *size, int num, int k_max) {
    /*
     * s_index and g_index have the same meanings in full_rank_finder function
     * row_index: indicator for row
     * col_index: indicator for column
     * i_g_index: indicator for group for the i-th row
     */
    int s_index, g_index, row_index, col_index, i_g_index, upper_index, lower_index;
    int full_rank_matrix_num = ((k_max + 1) * (k_max)) >> 1;
    int *init_full_rank;
    init_full_rank = (int *) malloc(full_rank_matrix_num * sizeof(int));
    for (int i = 0; i < num; ++i) {
        i_g_index = label[i];
        for (int k = 0; k < full_rank_matrix_num; ++k) {
            init_full_rank[k] = 1;
        }
        for (int j = 0; j < num; ++j) {
            s_index = index_matrix[i][j];
            g_index = label[s_index];
            // update pairwise upper rank matrix
            if ((g_index + 1) <= k_max) {
                for (int l_index = (g_index + 1); l_index <= k_max; ++l_index) {
                    if (i_g_index == g_index) {
                        row_index = group_relative_location[i] - cumsum_size[g_index];
                    } else if (i_g_index == l_index) {
                        row_index = group_relative_location[i] - cumsum_size[l_index] + size[g_index];
                    } else {
                        row_index = -1;
                    }
                    if (row_index != -1) {
                        upper_index = (l_index - g_index) + (((((k_max + 1) << 1) - 1 - g_index) * (g_index)) >> 1) - 1;
                        col_index = group_relative_location[s_index] - cumsum_size[g_index];
                        full_rank[upper_index][row_index][col_index] = init_full_rank[upper_index];
                        init_full_rank[upper_index] += 1;
                    }
                }
            }
            // update pairwise lower rank matrix
            if ((g_index - 1) >= 0) {
                for (int l_index = 0; l_index <= (g_index - 1); ++l_index) {
                    if (i_g_index == l_index) {
                        row_index = group_relative_location[i] - cumsum_size[l_index];
                    } else if (i_g_index == g_index) {
                        row_index = group_relative_location[i] - cumsum_size[g_index] + size[l_index];
                    } else {
                        row_index = -1;
                    }
                    if (row_index != -1) {
                        lower_index = (g_index - l_index) + (((((k_max + 1) << 1) - 1 - l_index) * (l_index)) >> 1) - 1;
                        col_index = group_relative_location[s_index] - cumsum_size[g_index] + size[l_index];
                        full_rank[lower_index][row_index][col_index] = init_full_rank[lower_index];
                        init_full_rank[lower_index] += 1;
                    }
                }
            }
        }
    }
}

void ball_divergence_array(double **bd_stat_array, int ***full_rank, int ***sub_rank, int *size, int K) {
    int s = 0;
    for (int i = 0; i < (K - 1); ++i) {
        for (int j = (i + 1); j < K; ++j) {
            ball_divergence2(bd_stat_array[s], full_rank[s], sub_rank[i], sub_rank[j], size[i], size[j]);
            s++;
        }
    }
}

void k_ball_divergence_from_by_sample_ball_divergence(double *kbd_stat, double **bd_stat_array,
                                                      int bd_stat_number, int k) {
    double kbd_maxsum_w0, kbd_maxsum_w1, kbd_max_w0, kbd_max_w1, kbd_sum_w0, kbd_sum_w1;
    double *bd_stat_w0_sum_array = (double *) malloc(k * sizeof(double));
    double *bd_stat_w1_sum_array = (double *) malloc(k * sizeof(double));
    double *bd_stat_w0_array = (double *) malloc(bd_stat_number * sizeof(double));
    double *bd_stat_w1_array = (double *) malloc(bd_stat_number * sizeof(double));
    // compute maxsum version statistic:
    for (int i = 0; i < k; ++i) {
        bd_stat_w0_sum_array[i] = 0;
        bd_stat_w1_sum_array[i] = 0;
    }
    int s = 0;
    for (int l = 0; l < (k - 1); ++l) {
        for (int m = (l + 1); m < k; ++m) {
            bd_stat_w0_sum_array[l] += bd_stat_array[s][0];
            bd_stat_w0_sum_array[m] += bd_stat_array[s][0];
            bd_stat_w1_sum_array[l] += bd_stat_array[s][1];
            bd_stat_w1_sum_array[m] += bd_stat_array[s][1];
            s++;
        }
    }
    quick_sort(bd_stat_w0_sum_array, k);
    quick_sort(bd_stat_w1_sum_array, k);
    kbd_maxsum_w0 = bd_stat_w0_sum_array[k - 1];
    kbd_maxsum_w1 = bd_stat_w1_sum_array[k - 1];
    free(bd_stat_w0_sum_array);
    free(bd_stat_w1_sum_array);

    // compute maximum K-1 version statistic:
    for (int i = 0; i < bd_stat_number; ++i) {
        bd_stat_w0_array[i] = bd_stat_array[i][0];
        bd_stat_w1_array[i] = bd_stat_array[i][1];
    }
    quick_sort(bd_stat_w0_array, bd_stat_number);
    quick_sort(bd_stat_w1_array, bd_stat_number);
    kbd_max_w0 = kbd_max_w1 = 0;
    for (int i = bd_stat_number - 1; i > (bd_stat_number - k); i--) {
        kbd_max_w0 += bd_stat_w0_array[i];
        kbd_max_w1 += bd_stat_w1_array[i];
    }
    free(bd_stat_w0_array);
    free(bd_stat_w1_array);

    // compute summation version:
    kbd_sum_w0 = kbd_sum_w1 = 0.0;
    for (int i = 0; i < bd_stat_number; ++i) {
        kbd_sum_w0 += bd_stat_array[i][0];
        kbd_sum_w1 += bd_stat_array[i][1];
    }

    kbd_stat[0] = kbd_sum_w0;
    kbd_stat[1] = kbd_sum_w1;
    kbd_stat[2] = kbd_max_w0;
    kbd_stat[3] = kbd_max_w1;
    kbd_stat[4] = kbd_maxsum_w0;
    kbd_stat[5] = kbd_maxsum_w1;
}

void KBD3(double *kbd_stat, double *pvalue, double *xy, int *size, int *n, int *k, int *R, int *thread) {
    int s, bd_stat_number = (((*k - 1) * (*k)) >> 1);
    int *cumsum_size = (int *) malloc(*k * sizeof(int));
    compute_cumsum_size(cumsum_size, size, k);
    int *label = (int *) malloc(*n * sizeof(int));
    s = 0;
    for (int i = 0; i < *k; ++i) {
        for (int j = 0; j < size[i]; ++j) {
            label[s++] = i;
        }
    }
    int *pairwise_size = (int *) malloc(bd_stat_number * sizeof(int));
    s = 0;
    for (int i = 0; i < (*k - 1); ++i) {
        for (int j = (i + 1); j < (*k); ++j) {
            pairwise_size[s] = size[i] + size[j];
            s++;
        }
    }
    int *group_relative_location = (int *) malloc(*n * sizeof(int));
    find_group_relative_location(group_relative_location, label, cumsum_size, *n, *k);

    double **distance_matrix = alloc_matrix(*n, *n);
    int **index_matrix = alloc_int_matrix(*n, *n);
    int ***sub_rank = alloc_int_square_matrix_list(size, *k);
    int ***full_rank = alloc_int_square_matrix_list(pairwise_size, bd_stat_number);

    distance2matrix(xy, distance_matrix, *n);
    for (int i = 0; i < *n; i++) {
        for (int j = 0; j < *n; j++) {
            index_matrix[i][j] = j;
        }
    }
    double *distance_matrix_copy = (double *) malloc(*n * sizeof(double));
    for (int i = 0; i < *n; i++) {
        memcpy(distance_matrix_copy, distance_matrix[i], *n * sizeof(double));
        quicksort(distance_matrix_copy, index_matrix[i], 0, *n - 1);
    }
    free(distance_matrix_copy);

    double **bd_stat_array = alloc_matrix(bd_stat_number, 2);
    sub_rank_finder(sub_rank, distance_matrix, index_matrix, label, group_relative_location, cumsum_size, *n, *k - 1);
    full_rank_finder(full_rank, distance_matrix, index_matrix, label, group_relative_location, cumsum_size, size, *n,
                     *k - 1);
    ball_divergence_array(bd_stat_array, full_rank, sub_rank, size, *k);
    k_ball_divergence_from_by_sample_ball_divergence(kbd_stat, bd_stat_array, bd_stat_number, *k);

    if (*R > 0) {
        double *permuted_kbd_sum_w0, *permuted_kbd_sum_w1, *permuted_kbd_max_w0, *permuted_kbd_max_w1, *permuted_kbd_maxsum_w0, *permuted_kbd_maxsum_w1;
        permuted_kbd_sum_w0 = (double *) malloc(*R * sizeof(double));
        permuted_kbd_sum_w1 = (double *) malloc(*R * sizeof(double));
        permuted_kbd_max_w0 = (double *) malloc(*R * sizeof(double));
        permuted_kbd_max_w1 = (double *) malloc(*R * sizeof(double));
        permuted_kbd_maxsum_w0 = (double *) malloc(*R * sizeof(double));
        permuted_kbd_maxsum_w1 = (double *) malloc(*R * sizeof(double));

        int not_parallel = *thread == 1 ? 1 : 0;
        int r;
        if (not_parallel) {
            for (r = 0; r < *R; ++r) {
                double kbd_stat_tmp[6];
                resample2(label, n);
                find_group_relative_location(group_relative_location, label, cumsum_size, *n, *k);
                sub_rank_finder(sub_rank, distance_matrix, index_matrix, label, group_relative_location,
                                cumsum_size, *n, *k - 1);
                full_rank_finder(full_rank, distance_matrix, index_matrix, label, group_relative_location,
                                 cumsum_size, size, *n, *k - 1);
                ball_divergence_array(bd_stat_array, full_rank, sub_rank, size, *k);
                k_ball_divergence_from_by_sample_ball_divergence(kbd_stat_tmp, bd_stat_array, bd_stat_number, *k);
                permuted_kbd_sum_w0[r] = kbd_stat_tmp[0];
                permuted_kbd_sum_w1[r] = kbd_stat_tmp[1];
                permuted_kbd_max_w0[r] = kbd_stat_tmp[2];
                permuted_kbd_max_w1[r] = kbd_stat_tmp[3];
                permuted_kbd_maxsum_w0[r] = kbd_stat_tmp[4];
                permuted_kbd_maxsum_w1[r] = kbd_stat_tmp[5];
            }
        } else {
            int **label_matrix = alloc_int_matrix(*R, *n);
            int **group_relative_location_matrix = alloc_int_matrix(*R, *n);
            resample2_matrix(label_matrix, label, *R, *n);
            for (r = 0; r < *R; ++r) {
                find_group_relative_location(group_relative_location_matrix[r], label_matrix[r], cumsum_size, *n, *k);
            }
#pragma omp parallel
            {
                int j_thread;
                double kbd_stat_thread[6];
                int ***sub_rank_thread = alloc_int_square_matrix_list(size, *k);
                int ***full_rank_thread = alloc_int_square_matrix_list(pairwise_size, bd_stat_number);
                double **bd_stat_array_thread = alloc_matrix(bd_stat_number, 2);
#pragma omp for
                for (j_thread = 0; j_thread < (*R); j_thread++) {
                    sub_rank_finder(sub_rank_thread, distance_matrix, index_matrix,
                                    label_matrix[j_thread], group_relative_location_matrix[j_thread],
                                    cumsum_size, *n, *k - 1);
                    full_rank_finder(full_rank_thread, distance_matrix, index_matrix,
                                     label_matrix[j_thread], group_relative_location_matrix[j_thread],
                                     cumsum_size, size, *n, *k - 1);
                    ball_divergence_array(bd_stat_array_thread, full_rank_thread, sub_rank_thread, size, *k);
                    k_ball_divergence_from_by_sample_ball_divergence(kbd_stat_thread, bd_stat_array_thread,
                                                                     bd_stat_number, *k);
                    permuted_kbd_sum_w0[j_thread] = kbd_stat_thread[0];
                    permuted_kbd_sum_w1[j_thread] = kbd_stat_thread[1];
                    permuted_kbd_max_w0[j_thread] = kbd_stat_thread[2];
                    permuted_kbd_max_w1[j_thread] = kbd_stat_thread[3];
                    permuted_kbd_maxsum_w0[j_thread] = kbd_stat_thread[4];
                    permuted_kbd_maxsum_w1[j_thread] = kbd_stat_thread[5];
                }
                free_int_square_matrix_list(sub_rank_thread, size, *k);
                free_int_square_matrix_list(full_rank_thread, pairwise_size, bd_stat_number);
                free_matrix(bd_stat_array_thread, bd_stat_number, 2);
            };
            r = *R;
            free_int_matrix(label_matrix, *R, *n);
            free_int_matrix(group_relative_location_matrix, *R, *n);
        }
        pvalue[0] = compute_pvalue(kbd_stat[0], permuted_kbd_sum_w0, r);
        pvalue[1] = compute_pvalue(kbd_stat[1], permuted_kbd_sum_w1, r);
        pvalue[2] = compute_pvalue(kbd_stat[2], permuted_kbd_max_w0, r);
        pvalue[3] = compute_pvalue(kbd_stat[3], permuted_kbd_max_w1, r);
        pvalue[4] = compute_pvalue(kbd_stat[4], permuted_kbd_maxsum_w0, r);
        pvalue[5] = compute_pvalue(kbd_stat[5], permuted_kbd_maxsum_w1, r);
        free(permuted_kbd_sum_w0);
        free(permuted_kbd_sum_w1);
        free(permuted_kbd_max_w0);
        free(permuted_kbd_max_w1);
        free(permuted_kbd_maxsum_w0);
        free(permuted_kbd_maxsum_w1);
    }

    free_int_square_matrix_list(full_rank, pairwise_size, bd_stat_number);
    free_int_square_matrix_list(sub_rank, size, *k);
    free_matrix(distance_matrix, *n, *n);
    free_int_matrix(index_matrix, *n, *n);
    free(label);
    free(group_relative_location);
    free(pairwise_size);
    free(cumsum_size);
    free_matrix(bd_stat_array, bd_stat_number, 2);
}