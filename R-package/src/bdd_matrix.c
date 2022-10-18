//
// Created by JinZhu on 2019/10/13.
//

#include "stdio.h"
#include "math.h"
#include "median.h"
#include "utilities.h"
#include "Ball_omp.h"

void bdd_matrix_bias_two_group(double *b_dd, double *x, int *n1_num, int *n2_num, int *nthread) {
    int not_parallel = (*nthread == 1 ? 1 : *nthread);
#ifdef Ball_OMP_H_
    if (not_parallel != 1) {
        omp_set_dynamic(0);
        if (*nthread <= 0) {
            omp_set_num_threads(omp_get_num_procs());
        } else {
            omp_set_num_threads(*nthread);
        }
    }
#endif

    const int n1 = *n1_num;
    const int n2 = *n2_num;
    const int num = n1 + n2;
    double x_t_vec[num];
    double x_n1_vec[n1];
    double x_n2_vec[n2];
    int r_t_vec[num];
    int r_n1_vec[n1];
    int r_n2_vec[n2];
    double **Dxy = alloc_matrix(num, num);
    distance2matrix(x, Dxy, num);

    int **x_rank = alloc_int_matrix(num, num);
    for (int i = 0; i < n1; i++) {
        for (int j = 0; j < n1; j++) {
            x_n1_vec[j] = Dxy[i][j];
        }
        quick_rank_min(x_n1_vec, r_n1_vec, n1);
        for (int j = 0; j < n1; j++) {
            x_rank[i][j] = r_n1_vec[j] - 1;
        }

        for (int j = 0; j < num; j++) {
            x_t_vec[j] = Dxy[i][j];
        }
        quick_rank_min(x_t_vec, r_t_vec, num);
        for (int j = n1; j < num; j++) {
            x_n2_vec[j - n1] = Dxy[i][j];
        }
        quick_rank_min(x_n2_vec, r_n2_vec, n2);
        for (int j = n1; j < num; j++) {
            x_rank[i][j] = r_t_vec[j] - r_n2_vec[j - n1];
        }
    }

    for (int i = n1; i < num; i++) {
        for (int j = n1; j < num; j++) {
            x_n2_vec[j - n1] = Dxy[i][j];
        }
        quick_rank_min(x_n2_vec, r_n2_vec, n2);
        for (int j = n1; j < num; j++) {
            x_rank[i][j] = r_n2_vec[j - n1] - 1;
        }

        for (int j = 0; j < num; j++) {
            x_t_vec[j] = Dxy[i][j];
        }
        quick_rank_min(x_t_vec, r_t_vec, num);
        for (int j = 0; j < n1; j++) {
            x_n1_vec[j] = Dxy[i][j];
        }
        quick_rank_min(x_n1_vec, r_n1_vec, n1);
        for (int j = 0; j < n1; j++) {
            x_rank[i][j] = r_t_vec[j] - r_n1_vec[j];
        }
    }
    free_matrix(Dxy, num, num);

    int c1 = n1 * n1 + n2 * n2;
    double c2 = 1.0 / c1;

    if (not_parallel) {
        int tmp_sum, s = 0;
        for (int i = 0; i < num; i++) {
            for (int j = i; j < num; j++) {
                tmp_sum = 0;
                for (int k = 0; k < num; k++) {
                    tmp_sum += MAX(x_rank[k][i], x_rank[k][j]);
                }
                tmp_sum = c1 - tmp_sum;
                b_dd[s++] = tmp_sum * c2;
            }
        }
    } else {
#pragma omp parallel
        {
            int tmp_sum, s, i, j, k;
#pragma omp for
            for (i = 0; i < num; i++) {
                for (j = i; j < num; j++) {
                    s = i * num - ((i * (i - 1)) >> 1) + j - i;
                    tmp_sum = 0;
                    for (k = 0; k < num; k++) {
                        tmp_sum += MAX(x_rank[k][i], x_rank[k][j]);
                    }
                    tmp_sum = c1 - tmp_sum;
                    b_dd[s] = tmp_sum * c2;
                }
            }
        }
    }
    free_int_matrix(x_rank, num, num);
}

double chi_square_weight(int rank_value, int num) {
    double weight = 0.0;
    double num_double = (double) num;
    double rank_double = (double) rank_value;
    if (rank_value == num) {
        // weight = (num_double * num_double) / (rank_double * (1.0 + num_double - rank_double));
        weight = 0.0;
    } else {
        weight = (num_double * num_double) / (rank_double * (num_double - rank_double));
    }
    return weight;
}

void bdd_matrix_bias(double *b_dd, double *x, int *n, int *nthread, int *weight_type) {
    int not_parallel = (*nthread == 1 ? 1 : *nthread);
#ifdef Ball_OMP_H_
    if (not_parallel != 1) {
        omp_set_dynamic(0);
        if (*nthread <= 0) {
            omp_set_num_threads(omp_get_num_procs());
        } else {
            omp_set_num_threads(*nthread);
        }
    }
#endif

    const int num = *n;
    double x_vec[num];
    int r_vec[num];
    int **x_rank = alloc_int_matrix(num, num);
    double **Dxy = alloc_matrix(num, num);
    distance2matrix(x, Dxy, num);

    for (int i = 0; i < num; i++) {
        for (int j = 0; j < num; j++) {
            x_vec[j] = Dxy[i][j];
        }
        quick_rank_max(x_vec, r_vec, num);
        for (int j = 0; j < num; j++) {
            x_rank[i][j] = r_vec[j];
        }
    }

    double median_value = 0.0;
    if (*weight_type == 4) {
        int s = 0;
        const int Dxy_vec_len = (((*n) * (*n - 1)) >> 1);
        double Dxy_vec[Dxy_vec_len];
        for (int u = 0; u < (*n - 1); u++)
        {
            for (int v = u + 1; v < *n; v++)
            {
                Dxy_vec[s++] = Dxy[u][v];
            }
        }
        median_value = find_median(Dxy_vec, (((*n) * (*n - 1)) >> 1));
        // printf("median distance: %f\n", median_value);
    }
    double num_double = (double) (*n);
    double **weight = alloc_matrix(num, num);
    for (int i = 0; i < num; i++) {
        for (int j = 0; j < num; j++) {
            double rank_double = (double) x_rank[i][j];
            if (*weight_type == 1) {
                weight[i][j] = 1.0;
            } else if (*weight_type == 2) {
                weight[i][j] = num_double / rank_double;
            } else if (*weight_type == 3) {
                weight[i][j] = chi_square_weight(x_rank[i][j], num);
            } else if (*weight_type == 4) {
                weight[i][j] = exp(-0.5 * Dxy[i][j] / median_value);
            }
        }
    }

    double w_vec[num];
    double ws_vec[num];
    double total_weight_sum[num];
    double **weight_sum = alloc_matrix(num, num);
    for (int i = 0; i < num; i++) {
        total_weight_sum[i] = 0.0;
        for (int j = 0; j < num; j++) {
            x_vec[j] = Dxy[i][j];
            w_vec[j] = weight[i][j];
            ws_vec[j] = 0.0;
            total_weight_sum[i] += w_vec[j];
        }

        quick_rank_min_weight_sum(x_vec, w_vec, ws_vec, num);
        for (int j = 0; j < num; j++) {
            weight_sum[i][j] = ws_vec[j];
        }
    }
    free_matrix(Dxy, num, num);

    double c1 = 0.0;
    for (int i = 0; i < num; i++) {
        c1 += total_weight_sum[i];
    }
    double c2 = 1.0 / (num * num);

    if (not_parallel) {
        int s = 0;
        double tmp_sum;
        for (int i = 0; i < num; i++) {
            for (int j = i; j < num; j++) {
                tmp_sum = 0.0;
                for (int k = 0; k < num; k++) {
                    tmp_sum += MAX(weight_sum[k][i], weight_sum[k][j]);
                }
                b_dd[s++] = (c1 - tmp_sum) * c2;
            }
        }
    } else {
#pragma omp parallel
        {
            int s, i, j, k;
            double tmp_sum;
#pragma omp for
            for (i = 0; i < num; i++) {
                for (j = i; j < num; j++) {
                    s = i * num - ((i * (i - 1)) >> 1) + j - i;
                    tmp_sum = 0.0;
                    for (k = 0; k < num; k++) {
                        tmp_sum += MAX(weight_sum[k][i], weight_sum[k][j]);
                    }
                    b_dd[s] = (c1 - tmp_sum) * c2;
                }
            }
        };
    }

    free_int_matrix(x_rank, num, num);
    free_matrix(weight, num, num);
    free_matrix(weight_sum, num, num);
}

void joint_kernel_matrix_bias(double *kernel, double *x, int *k, int *n, int *nthread, int *weight_type) {
    int not_parallel = (*nthread == 1 ? 1 : *nthread);
#ifdef Ball_OMP_H_
    if (not_parallel != 1) {
        omp_set_dynamic(0);
        if (*nthread <= 0) {
            omp_set_num_threads(omp_get_num_procs());
        } else {
            omp_set_num_threads(*nthread);
        }
    }
#endif

    double ***Dx;

    // vector to matrix
    Dx = alloc_3d_matrix(*n, *n, *k);

    // convert x to 3D distance matrix:
    distance2matrix3d(x, Dx, *n, *k);

    const int num = *n;

    double c1 = 1.0;
    double c2 = 1.0 / (num * num);

    if (not_parallel) {
        int u = 0, tmp_sum, inclusion_count;
        for (int i = 0; i < num; i++) {
            for (int j = i; j < num; j++) {
                inclusion_count = 0;
                for (int s = 0; s < num; s++) {
                    for (int t = 0; t < num; t++) {
                        tmp_sum = 1;
                        for (int m = 0; m < *k; m++) {
                            // tmp_sum *= MAX(Dx[s][i][m], Dx[s][j][m]) <= Dx[s][t][m];
                            tmp_sum *= (Dx[s][i][m] <= Dx[s][t][m]) * (Dx[s][j][m] <= Dx[s][t][m]);
                            if (tmp_sum == 0) {
                                continue;
                            }
                        }
                        inclusion_count += tmp_sum;
                    }
                }
                kernel[u++] = c1 * c2 * inclusion_count;
            }
        }
    } else {
#pragma omp parallel
        {
            int i, j, s, t, m;
            int u = 0, tmp_sum, inclusion_count;
#pragma omp for
            for (i = 0; i < num; i++) {
                for (j = i; j < num; j++) {
                    u = i * num - ((i * (i - 1)) >> 1) + j - i;
                    inclusion_count = 0;
                    for (int s = 0; s < num; s++) {
                        for (int t = 0; t < num; t++) {
                            tmp_sum = 1;
                            for (int m = 0; m < *k; m++) {
                                tmp_sum *= (Dx[s][i][m] <= Dx[s][t][m]) * (Dx[s][j][m] <= Dx[s][t][m]);
                                if (tmp_sum == 0) {
                                    continue;
                                }
                            }
                            inclusion_count += tmp_sum;
                        }
                    }
                    kernel[u++] = c1 * c2 * inclusion_count;
                }
            }
        };
    }

    free_3d_matrix(Dx, *n, *n);
}

void cross_kernel_matrix_bias_crude(double *kernel, double *x, int *k, int *n, int *nthread, int *weight_type) {
    int not_parallel = (*nthread == 1 ? 1 : *nthread);
#ifdef Ball_OMP_H_
    if (not_parallel != 1) {
        omp_set_dynamic(0);
        if (*nthread <= 0) {
            omp_set_num_threads(omp_get_num_procs());
        } else {
            omp_set_num_threads(*nthread);
        }
    }
#endif

    // convert x to 3D distance matrix:
    double ***Dx;
    Dx = alloc_3d_matrix(*n, *n, *k);
    distance2matrix3d(x, Dx, *n, *k);

    const int num = *n;
    const int K = *k;
    double c1 = 1.0;
    double c2 = 1.0 / (num * num);
    int seq_pow_num[K]; 
    for (int k_tmp = 0; k_tmp < K; k_tmp++)
    {
        if (k_tmp == 0) {
            seq_pow_num[0] = 1;
        } else {
            seq_pow_num[k_tmp] = seq_pow_num[k_tmp - 1] * num;
        }
    }
    int combination_type_num = num * seq_pow_num[K - 1];  // num^K
    if (not_parallel) {
        int s = 0;
        for (int u = 0; u < num; u++) {
            // find the balls that include u.
            int index1[num * num], index2[num * num];
            int complete_inclusion_count = 0;
            for (int i = 0; i < num; i++) {
                for (int j = 0; j < num; j++) {
                    int complete_inclusion = 1;
                    for (int k_tmp = 0; k_tmp < K; k_tmp++) {
                        if (Dx[i][u][k_tmp] > Dx[i][j][k_tmp]) {
                            complete_inclusion = 0; 
                            continue;
                        }
                    }
                    if (complete_inclusion) {
                        index1[complete_inclusion_count] = i;
                        index2[complete_inclusion_count] = j;
                        complete_inclusion_count++;
                    }
                }
            }

            // compute the kernel values
            int combination_tuple[K];
            for (int i = 0; i < combination_type_num; i++) {
                int cross_inclusion_count = 0;
                for (int k_tmp = 0; k_tmp < K; k_tmp++) {
                    combination_tuple[K - k_tmp - 1] = (i / seq_pow_num[k_tmp]) % num;
                }
                for (int j = 0; j < complete_inclusion_count; j++) {
                    int marginal_inclusion = 1;
                    for (int k_tmp = 0; k_tmp < K; k_tmp++) {
                        if (Dx[index1[j]][combination_tuple[k_tmp]][k_tmp] > Dx[index1[j]][index2[j]][k_tmp]) {
                            marginal_inclusion = 0;
                            continue;
                        }
                    }
                    cross_inclusion_count += marginal_inclusion;
                }
                kernel[s++] = c2 * ((double) cross_inclusion_count); 
            }
        }
    }

    free_3d_matrix(Dx, *n, *n);
}