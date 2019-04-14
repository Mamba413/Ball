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
                     const int *group_relative_location, const int *cumsum_size, const int *size,
                     int num, int max_k) {
    /*
     * g_index: indicator for group
     * s_index: indicator for order
     */
    int s_index, g_index;
    for (int i = 0; i < num; ++i) {
        int i_index = label[i], i_cumsum_size = cumsum_size[i_index];
        int i_group_relative_location = group_relative_location[i];
        int rank_value = 1;
        for (int j = 0; j < num; ++j) {
            s_index = index_matrix[i][j];
            g_index = label[s_index];
            if (i_index == g_index) {
                sub_rank[i_index][i_group_relative_location - i_cumsum_size][group_relative_location[s_index] - i_cumsum_size] = rank_value;
                rank_value++;
            }
        }
    }
}

/**
 * @inherit sub_rank_finder
 * @param size : [4, 5, 4]
 */
void full_rank_finder(int ***full_rank, double **distance_matrix, int **index_matrix, const int *label,
                      const int *group_relative_location, const int *cumsum_size, const int *size, int num, int k_max) {
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

void asymptotic_ball_divergence(double *asymptotic_bd_stat, int ***full_rank, int ***sub_rank, int *size,
                                int K, int bd_stat_number) {
    for (int i = 0; i < 2; ++i) {
        asymptotic_bd_stat[i] = 0.0;
    }
    double **bd_stat_array = alloc_matrix(bd_stat_number, 2);
    int s = 0;
    for (int i = 0; i < (K - 1); ++i) {
        for (int j = (i + 1); j < K; ++j) {
            ball_divergence2(bd_stat_array[s], full_rank[s], sub_rank[i], sub_rank[j], size[i], size[j]);
            asymptotic_bd_stat[0] += bd_stat_array[s][0] * size[i] * size[j] / (size[i] + size[j]);
            asymptotic_bd_stat[1] += bd_stat_array[s][1] * size[i] * size[j] / (size[i] + size[j]);
            s++;
        }
    }
}

void compute_pairwise_size(int *pairwise_size, const int *size, const int *k) {
    int s = 0;
    for (int i = 0; i < (*k - 1); ++i) {
        for (int j = (i + 1); j < (*k); ++j) {
            pairwise_size[s++] = size[i] + size[j];
        }
    }
}

/* Comparison function. Receives two generic (void) pointers to the items under comparison. */
int compare_ints(const void *p, const void *q) {
    int x = *(const int *) p;
    int y = *(const int *) q;

    /* Avoid return x - y, which can cause undefined behaviour
       because of signed integer overflow. */
    if (x < y)
        return -1;  // Return -1 if you want ascending, 1 if you want descending order.
    else if (x > y)
        return 1;   // Return 1 if you want ascending, -1 if you want descending order.
    return 0;
}

/* Sort an array of n integers, pointed to by a. */
void sort_ints(int *a, size_t n) {
    qsort(a, n, sizeof *a, &compare_ints);
}

void KBD3(double *kbd_stat, double *pvalue, double *xy, int *size, int *n, int *k, const int *R, const int *thread) {
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
    compute_pairwise_size(pairwise_size, size, k);

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
    sub_rank_finder(sub_rank, distance_matrix, index_matrix, label, group_relative_location, cumsum_size, size, *n,
                    *k - 1);
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
                                cumsum_size, size, *n, *k - 1);
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
                                    cumsum_size, size, *n, *k - 1);
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

void bd_gwas_test(double *bd_stat, double *permuted_bd_stat, double *pvalue, double *xy, const int *snp,
                  const int *n, const int *p, const int *R, const int *nthread) {
#ifdef Ball_OMP_H_
    omp_set_dynamic(0);
    if (*nthread <= 0) {
        omp_set_num_threads(omp_get_num_procs());
    } else {
        omp_set_num_threads(*nthread);
    }
#endif

    int *k_vector = (int *) malloc(*p * sizeof(int));
    int *snp_vector = (int *) malloc(*n * sizeof(int)), **snp_matrix = alloc_int_matrix(*p, *n);
    int **size_list = (int **) malloc(*p * sizeof(int *));
    int **cumsum_size_list = (int **) malloc(*p * sizeof(int *));
    int s = 0;
    for (int i = 0; i < *p; ++i) {
        for (int j = 0; j < *n; ++j) {
            snp_matrix[i][j] = snp[s++];
        }
        memcpy(snp_vector, snp_matrix[i], *n * sizeof(int));
        sort_ints(snp_vector, (size_t) *n);
        k_vector[i] = 1;
        for (int j = 1; j < *n; ++j) {
            if (snp_vector[j] != snp_vector[j - 1]) {
                k_vector[i] += 1;
            }
        }
        size_list[i] = (int *) malloc(k_vector[i] * sizeof(int));
        cumsum_size_list[i] = (int *) malloc(k_vector[i] * sizeof(int));
        cumsum_size_list[i][0] = 0;
        int t = 1;
        for (int j = 1; j < *n; ++j) {
            if (snp_vector[j] != snp_vector[j - 1]) {
                cumsum_size_list[i][t] = j;
                size_list[i][t - 1] = j - cumsum_size_list[i][t - 1];
                t++;
            }
        }
        size_list[i][k_vector[i] - 1] = *n - cumsum_size_list[i][k_vector[i] - 1];
    }
    free(snp_vector);

    int **index_matrix = alloc_int_matrix(*n, *n);
    double **distance_matrix = alloc_matrix(*n, *n);
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

    int *bd_stat_number_vector = (int *) malloc(*p * sizeof(int));
    int **pairwise_size = (int **) malloc(*p * sizeof(int *));
    for (int i = 0; i < *p; ++i) {
        bd_stat_number_vector[i] = (((k_vector[i] - 1) * k_vector[i]) >> 1);
        pairwise_size[i] = (int *) malloc(bd_stat_number_vector[i] * sizeof(int));
        compute_pairwise_size(pairwise_size[i], size_list[i], &k_vector[i]);
    }

    // compute statistic:
    int **snp_group_relative_location = alloc_int_matrix(*p, *n);
    for (int i = 0; i < *p; ++i) {
        find_group_relative_location(snp_group_relative_location[i], snp_matrix[i],
                                     cumsum_size_list[i], *n, k_vector[i]);
    }
    double **asymptotic_bd_stat_array = alloc_matrix(*p, 2);
#pragma omp parallel
    {
        int i_thread;
        int ***sub_rank, ***full_rank;
#pragma omp for
        for (i_thread = 0; i_thread < *p; ++i_thread) {
            sub_rank = alloc_int_square_matrix_list(size_list[i_thread], k_vector[i_thread]);
            full_rank = alloc_int_square_matrix_list(pairwise_size[i_thread], bd_stat_number_vector[i_thread]);
            sub_rank_finder(sub_rank, distance_matrix, index_matrix, snp_matrix[i_thread],
                            snp_group_relative_location[i_thread], cumsum_size_list[i_thread], size_list[i_thread],
                            *n, k_vector[i_thread] - 1);
            full_rank_finder(full_rank, distance_matrix, index_matrix, snp_matrix[i_thread],
                             snp_group_relative_location[i_thread], cumsum_size_list[i_thread], size_list[i_thread],
                             *n, k_vector[i_thread] - 1);
            asymptotic_ball_divergence(asymptotic_bd_stat_array[i_thread], full_rank, sub_rank, size_list[i_thread],
                                       k_vector[i_thread], bd_stat_number_vector[i_thread]);
            free_int_square_matrix_list(sub_rank, size_list[i_thread], k_vector[i_thread]);
            free_int_square_matrix_list(full_rank, pairwise_size[i_thread], bd_stat_number_vector[i_thread]);
        }
    };
    free_int_matrix(snp_matrix, *p, *n);
    free_int_matrix(snp_group_relative_location, *p, *n);
    for (int i = 0; i < *p; ++i) {
        free(size_list[i]);
        free(cumsum_size_list[i]);
        free(pairwise_size[i]);
    }
    free(size_list);
    free(cumsum_size_list);
    free(pairwise_size);
    free(bd_stat_number_vector);
    for (int i = 0; i < *p; ++i) {
        bd_stat[i] = asymptotic_bd_stat_array[i][0];
        bd_stat[i + *p] = asymptotic_bd_stat_array[i][1];
    }

    if (*R > 0) {
        int *k_vector_tmp = (int *) malloc(*p * sizeof(int));
        memcpy(k_vector_tmp, k_vector, *p * sizeof(int));
        sort_ints(k_vector_tmp, (size_t) *p);
        int min_k = k_vector_tmp[0];
        int unique_k_num = 1;
        for (int i = 1; i < *p; ++i) {
            if (k_vector_tmp[i] != k_vector_tmp[i - 1]) {
                unique_k_num++;
            }
        }
        int *unique_k_vector = (int *) malloc(unique_k_num * sizeof(int));
        unique_k_vector[0] = k_vector_tmp[0];
        s = 1;
        for (int i = 1; i < *p; ++i) {
            if (k_vector_tmp[i] != k_vector_tmp[i - 1]) {
                unique_k_vector[s++] = k_vector_tmp[i];
            }
        }
        free(k_vector_tmp);

        // compute permuted statistic when each sample size is balanced
        int batch_size;
        if (*p > 20000) {
            batch_size = 20000 * (*p / 20000);
        } else {
            batch_size = 20000;
        }
        if (*R < batch_size) {
            batch_size = *R;
        }
        int fix_batch_size = batch_size;
        int add_round = (*R % batch_size) == 0 ? 0 : 1;
        int batch_round = (*R / batch_size) + add_round;
        int largeR = *R - (batch_round - add_round) * (batch_size);
        double **permuted_asymptotic_bd_stat_batch = alloc_matrix(batch_size, 2);
        double **permuted_asymptotic_bd_stat_matrix = alloc_matrix((unique_k_num << 1), *R);
        int *label = (int *) malloc(*n * sizeof(int));
        int **label_matrix = alloc_int_matrix(batch_size, *n);
        int **group_relative_location_matrix = alloc_int_matrix(batch_size, *n);
        for (int k = 0; k < unique_k_num; ++k) {
            int row_index = 2 * k;
            int permuted_k = unique_k_vector[k], permuted_bd_stat_number = ((permuted_k - 1) * permuted_k) >> 1;
            int *permuted_size = (int *) malloc(permuted_k * sizeof(int));
            int tmp_cumsum_size = 0;
            for (int j = 0; j < (permuted_k - 1); ++j) {
                permuted_size[j] = *n / permuted_k;
                tmp_cumsum_size += permuted_size[j];
            }
            permuted_size[permuted_k - 1] = *n - tmp_cumsum_size;
            int *permuted_cumsum_size = (int *) malloc(permuted_k * sizeof(int));
            compute_cumsum_size(permuted_cumsum_size, permuted_size, &permuted_k);
            int *permuted_pairwise_size = (int *) malloc(permuted_bd_stat_number * sizeof(int));
            compute_pairwise_size(permuted_pairwise_size, permuted_size, &permuted_k);
            s = 0;
            for (int i = 0; i < permuted_k; ++i) {
                for (int j = 0; j < permuted_size[i]; ++j) {
                    label[s++] = i;
                }
            }
            // use several batch to conduct permutation to prevent memory insufficient
            for (int round = 0; round < batch_round; ++round) {
                if ((round == (batch_round - 1)) && (largeR > 0)) {
                    batch_size = largeR;
                }
                resample2_matrix(label_matrix, label, batch_size, *n);
                for (int r = 0; r < batch_size; ++r) {
                    find_group_relative_location(group_relative_location_matrix[r], label_matrix[r],
                                                 permuted_cumsum_size, *n, permuted_k);
                }
#pragma omp parallel
                {
                    int r_thread;
                    int ***sub_rank_thread = alloc_int_square_matrix_list(permuted_size, permuted_k);
                    int ***full_rank_thread = alloc_int_square_matrix_list(permuted_pairwise_size,
                                                                           permuted_bd_stat_number);
#pragma omp for
                    for (r_thread = 0; r_thread < batch_size; ++r_thread) {
                        find_group_relative_location(group_relative_location_matrix[r_thread], label_matrix[r_thread],
                                                     permuted_cumsum_size, *n, permuted_k);
                        sub_rank_finder(sub_rank_thread, distance_matrix, index_matrix, label_matrix[r_thread],
                                        group_relative_location_matrix[r_thread],
                                        permuted_cumsum_size, permuted_size, *n, permuted_k - 1);
                        full_rank_finder(full_rank_thread, distance_matrix, index_matrix, label_matrix[r_thread],
                                         group_relative_location_matrix[r_thread],
                                         permuted_cumsum_size, permuted_size, *n, permuted_k - 1);
                        asymptotic_ball_divergence(permuted_asymptotic_bd_stat_batch[r_thread], full_rank_thread,
                                                   sub_rank_thread, permuted_size, permuted_k, permuted_bd_stat_number);
                    }
                    free_int_square_matrix_list(sub_rank_thread, permuted_size, permuted_k);
                    free_int_square_matrix_list(full_rank_thread, permuted_pairwise_size, permuted_bd_stat_number);
                };
                int start = round * fix_batch_size;
                for (int r = 0; r < batch_size; ++r) {
                    permuted_asymptotic_bd_stat_matrix[row_index][r + start] = permuted_asymptotic_bd_stat_batch[r][0];
                    permuted_asymptotic_bd_stat_matrix[row_index + 1][r + start] = permuted_asymptotic_bd_stat_batch[r][1];
                    permuted_bd_stat[(row_index * *R) + start + r] = permuted_asymptotic_bd_stat_batch[r][0];
                    permuted_bd_stat[((row_index + 1) * *R) + start + r] = permuted_asymptotic_bd_stat_batch[r][1];
                }
            }
            batch_size = fix_batch_size;
            free(permuted_size);
            free(permuted_cumsum_size);
            free(permuted_pairwise_size);
        }
        free_matrix(permuted_asymptotic_bd_stat_batch, batch_size, 2);
        free_int_matrix(label_matrix, batch_size, *n);
        free_int_matrix(group_relative_location_matrix, batch_size, *n);
        free(unique_k_vector);
        free(label);

        // compute p-value
        for (int i = 0; i < *p; ++i) {
            int row_index = 2 * (k_vector[i] - min_k);
            pvalue[i] = compute_pvalue(asymptotic_bd_stat_array[i][0], permuted_asymptotic_bd_stat_matrix[row_index],
                                       *R);
            pvalue[*p + i] = compute_pvalue(asymptotic_bd_stat_array[i][1],
                                            permuted_asymptotic_bd_stat_matrix[row_index + 1], *R);
        }
        free_matrix(permuted_asymptotic_bd_stat_matrix, (unique_k_num << 1), *R);
    }

    free_matrix(distance_matrix, *n, *n);
    free_int_matrix(index_matrix, *n, *n);
    free_matrix(asymptotic_bd_stat_array, *p, 2);
    free(k_vector);
}