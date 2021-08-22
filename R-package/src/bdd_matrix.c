//
// Created by JinZhu on 2019/10/13.
//

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

void bdd_matrix_bias(double *b_dd, double *x, int *n, int *nthread) {
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
        quick_rank_min(x_vec, r_vec, num);
        for (int j = 0; j < num; j++) {
            x_rank[i][j] = r_vec[j];
        }
    }
    free_matrix(Dxy, num, num);


    int c1 = num * (num + 1);
    double c2 = 1.0 / (num * num);

    if (not_parallel) {
        int tmp_sum, s = 0;
        for (int i = 0; i < num; i++) {
            for (int j = i; j < num; j++) {
                tmp_sum = 0;
                for (int k = 0; k < num; k++) {
                    tmp_sum += MAX(x_rank[k][i], x_rank[k][j]);
                }
                b_dd[s++] = (c1 - tmp_sum) * c2;
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
                    b_dd[s] = (c1 - tmp_sum) * c2;
                }
            }
        };
    }

    free_int_matrix(x_rank, num, num);
}