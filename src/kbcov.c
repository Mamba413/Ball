#include "Ball_omp.h"

#include "stdio.h"
#include "stdlib.h"
#include "math.h"
#include "utilities.h"


void K_Ball_Covariance(double *kbcov_stat, double ***Dx, int ***Rx, int **i_perm, int *n, int *k) {
    int i, j, t, var_index, compute_joint_index;
    double p_all, p_prod, p_inv_prod;
    double *p_k_array;
    double bcov_weight0 = 0.0, bcov_weight_prob = 0.0, bcov_weight_hhg = 0.0, bcov_fixed_ball = 0.0;
    int hhg_include_ball_index;
    double hhg_ball_num = 0.0, pow_n_k = pow(*n, *k);

    p_k_array = (double *) malloc(*k * sizeof(double));

    for (i = 0; i < (*n); i++) {
        for (j = 0; j < (*n); j++) {
            // init the computation for each Ball
            p_all = *n;
            for (var_index = 0; var_index < (*k); var_index++) { p_k_array[var_index] = 0.0; }
            p_prod = 1.0, p_inv_prod = 1.0;
            hhg_include_ball_index = 1;

            // to compute P_{i, j}^{\mu_{1}}, ..., P_{i, j}^{\mu_{*k}} (marginal):
            compute_joint_index = 1;
            for (var_index = 0; var_index < (*k); var_index++) {
                p_k_array[var_index] = Rx[i_perm[var_index][i]][i_perm[var_index][j]][var_index];
                if (compute_joint_index == 1) {
                    compute_joint_index = p_k_array[var_index] == 0 ? 0 : 1;
                }
            }
            // to compute P_{i, j}^{\mu_{1}, ..., \mu_{*k}} (joint):
            if (compute_joint_index) {
                for (t = 0; t < (*n); t++) {
                    for (var_index = 0; var_index < (*k); var_index++) {
                        if (Dx[i_perm[var_index][i]][i_perm[var_index][j]][var_index] <
                            Dx[i_perm[var_index][i]][i_perm[var_index][t]][var_index]) {
                            p_all -= 1.0;
                            break;
                        }
                    }
                }
            } else {
                p_all = 0.0;
            }

            for (var_index = 0; var_index < (*k); var_index++) {
                p_prod *= p_k_array[var_index];
                if (hhg_include_ball_index == 1) {
                    hhg_include_ball_index = (p_k_array[var_index] > 2 && p_k_array[var_index] < *n) ? 1 : 0;
                }
                p_inv_prod *= (*n - p_k_array[var_index]);
            }
            p_all = p_all / *n;
            p_prod = p_prod / pow_n_k;
            bcov_fixed_ball = (p_all - p_prod) * (p_all - p_prod);

            bcov_weight0 += bcov_fixed_ball;
            bcov_weight_prob += bcov_fixed_ball / p_prod;
            if (hhg_include_ball_index == 1) {
                p_inv_prod = p_inv_prod / pow_n_k;
                bcov_weight_hhg += bcov_fixed_ball / (p_prod * p_inv_prod);
                hhg_ball_num += 1.0;
            }
        }
    }
    kbcov_stat[0] = bcov_weight0 / (1.0 * (*n) * (*n));
    kbcov_stat[1] = bcov_weight_prob / (1.0 * (*n) * (*n));
    kbcov_stat[2] = bcov_weight_hhg / (hhg_ball_num);

    free(p_k_array);
}

void kbcov_test_single(double *kbcov_stat, double *pvalue, double *x, int *k, int *n, int *R, int *dst) {
    int i, j, **i_perm, ***Rx;
    // create a 3D matrix to store distance matrix:
    double ***Dx;

    Dx = alloc_3d_matrix(*n, *n, *k);
    Rx = alloc_3d_int_matrix(*n, *n, *k);
    i_perm = alloc_int_matrix(*k, *n);

    // convert x to 3D distance matrix:
    if (*dst == 1) { vector_2_3dmatrix(x, Dx, *n, *n, *k, 1); }

    // compute rank 3D distance matrix row by row and one by one
    rank_matrix_3d(Dx, *n, *k, Rx);

    // create 2D permute index matrix:
    for (j = 0; j < *k; j++) {
        for (i = 0; i < *n; i++) { i_perm[j][i] = i; }
    }

    // compute kbcov value
    K_Ball_Covariance(kbcov_stat, Dx, Rx, i_perm, n, k);
    //printf("KBcov: %f \n", kbcov_stat[0]);
    //printf("KBcov-Probability: %f \n", kbcov_stat[1]);
    //printf("KBcov-HHG: %f \n", kbcov_stat[2]);

    // Hypothesis test
    if (*R > 0) {
        double bcov_tmp[3], *permuted_bcov_weight0, *permuted_bcov_weight_prob, *permuted_bcov_weight_hhg;
        permuted_bcov_weight0 = (double *) malloc(*R * sizeof(double));
        permuted_bcov_weight_prob = (double *) malloc(*R * sizeof(double));
        permuted_bcov_weight_hhg = (double *) malloc(*R * sizeof(double));

        for (i = 0; i < *R; i++) {
            // stop permutation if user stop it manually:
            if (pending_interrupt()) {
                print_stop_message();
                break;
            }
            resample_matrix(i_perm, k, n);
            K_Ball_Covariance(bcov_tmp, Dx, Rx, i_perm, n, k);
            permuted_bcov_weight0[i] = bcov_tmp[0];
            permuted_bcov_weight_prob[i] = bcov_tmp[1];
            permuted_bcov_weight_hhg[i] = bcov_tmp[2];
        }
        pvalue[0] = compute_pvalue(kbcov_stat[0], permuted_bcov_weight0, i);
        pvalue[1] = compute_pvalue(kbcov_stat[1], permuted_bcov_weight_prob, i);
        pvalue[2] = compute_pvalue(kbcov_stat[2], permuted_bcov_weight_hhg, i);

        //printf("%f, %f, %f\n", pvalue[0], pvalue[1], pvalue[2]);
        free(permuted_bcov_weight0);
        free(permuted_bcov_weight_prob);
        free(permuted_bcov_weight_hhg);
    }

    free_3d_matrix(Dx, *n, *n);
    free_3d_int_matrix(Rx, *n, *n);
    free_int_matrix(i_perm, *k, *n);
}

void kbcov_test_parallel(double *kbcov_stat, double *pvalue, double *x, int *k, int *n, int *R, int *dst) {
    int i, j, **i_perm, ***Rx;
    // create a 3D matrix to store distance matrix:
    double ***Dx;

    Dx = alloc_3d_matrix(*n, *n, *k);
    Rx = alloc_3d_int_matrix(*n, *n, *k);
    i_perm = alloc_int_matrix(*k, *n);

    // convert x to 3D distance matrix:
    if (*dst == 1) { vector_2_3dmatrix(x, Dx, *n, *n, *k, 1); }

    // compute rank 3D distance matrix row by row and one by one
    rank_matrix_3d(Dx, *n, *k, Rx);

    // create 2D permute index matrix:
    for (j = 0; j < *k; j++) {
        for (i = 0; i < *n; i++) { i_perm[j][i] = i; }
    }

    // compute kbcov value
    K_Ball_Covariance(kbcov_stat, Dx, Rx, i_perm, n, k);

    // permutation test
    if (*R > 0) {
        double *permuted_bcov_weight0, *permuted_bcov_weight_prob, *permuted_bcov_weight_hhg;
        permuted_bcov_weight0 = (double *) malloc(*R * sizeof(double));
        permuted_bcov_weight_prob = (double *) malloc(*R * sizeof(double));
        permuted_bcov_weight_hhg = (double *) malloc(*R * sizeof(double));

#pragma omp parallel
        {
            int i_thread;
            double bcov_tmp[3];
#pragma omp for
            for (i_thread = 0; i_thread < *R; i_thread++) {
#pragma omp critical
                {
                    resample_matrix(i_perm, k, n);
                }
                K_Ball_Covariance(bcov_tmp, Dx, Rx, i_perm, n, k);
                permuted_bcov_weight0[i_thread] = bcov_tmp[0];
                permuted_bcov_weight_prob[i_thread] = bcov_tmp[1];
                permuted_bcov_weight_hhg[i_thread] = bcov_tmp[2];
            }
        }

        pvalue[0] = compute_pvalue(kbcov_stat[0], permuted_bcov_weight0, *R);
        pvalue[1] = compute_pvalue(kbcov_stat[1], permuted_bcov_weight_prob, *R);
        pvalue[2] = compute_pvalue(kbcov_stat[2], permuted_bcov_weight_hhg, *R);

        free(permuted_bcov_weight0);
        free(permuted_bcov_weight_prob);
        free(permuted_bcov_weight_hhg);
    }
    free_3d_matrix(Dx, *n, *n);
    free_3d_int_matrix(Rx, *n, *n);
    free_int_matrix(i_perm, *k, *n);
}

void kbcov_test(double *kbcov_stat, double *pvalue, double *x, int *k, int *n, int *R, int *dst, int *thread) {
#ifdef Ball_OMP_H_
    omp_set_num_threads(*thread);
#endif

    if (*thread > 1) {
        kbcov_test_parallel(kbcov_stat, pvalue, x, k, n, R, dst);
    } else {
        kbcov_test_single(kbcov_stat, pvalue, x, k, n, R, dst);
    }
}