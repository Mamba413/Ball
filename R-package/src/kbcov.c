#include "Ball_omp.h"

#include "stdio.h"
#include "stdlib.h"
#include "utilities.h"

#ifdef R_BUILD
#include "R.h"
#include "Rinternals.h"
#else

#include "time.h"

#endif

void K_U_Ball_Covariance(double *kbcov_stat, double **x, int **i_perm, int *n, int *k) {
    // TODO...
}

void K_Ball_Covariance(double *kbcov_stat, double ***Dx, int ***Rx, int **i_perm, const int *n, const int *k) {
    int i, j, t, var_index, compute_joint_index;
    double p_all, p_prod, p_inv_prod, at_least_point_number;
    double *p_k_array;
    double bcov_weight0 = 0.0, bcov_weight_prob = 0.0, bcov_weight_hhg = 0.0, bcov_fixed_ball;
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

            at_least_point_number = i == j ? 1 : 2;
            // to compute P_{i, j}^{\mu_{1}}, ..., P_{i, j}^{\mu_{*k}} (marginal):
            compute_joint_index = 1;
            for (var_index = 0; var_index < (*k); var_index++) {
                p_k_array[var_index] = Rx[i_perm[var_index][i]][i_perm[var_index][j]][var_index];
                if (compute_joint_index == 1) {
                    // Note: for each variable, if there is any only contain the at-least-point-number, we can derive that the joint number directly.
                    compute_joint_index = p_k_array[var_index] == at_least_point_number ? 0 : 1;
                }
            }
            // Compute P_{i, j}^{\mu_{1}, ..., \mu_{*k}} (joint), O(n^3):
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
                p_all = at_least_point_number;
            }

            for (var_index = 0; var_index < (*k); var_index++) {
                p_prod *= p_k_array[var_index];
                if (hhg_include_ball_index == 1) {
                    hhg_include_ball_index = (p_k_array[var_index] < *n) ? 1 : 0;
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
    kbcov_stat[2] = hhg_ball_num > 0 ? (bcov_weight_hhg / hhg_ball_num) : 0.0;

    free(p_k_array);
}

void KBCOV(double *kbcov_stat, double *pvalue, double *x, int *k, int *n, int *R, const int *thread) {
    int **i_perm, ***Rx;
    // create a 3D matrix to store distance matrix:
    double ***Dx;

    Dx = alloc_3d_matrix(*n, *n, *k);
    Rx = alloc_3d_int_matrix(*n, *n, *k);
    i_perm = alloc_int_matrix(*k, *n);

    // convert x to 3D distance matrix:
    distance2matrix3d(x, Dx, *n, *k);

    // compute rank 3D distance matrix row by row and one by one
    rank_matrix_3d(Dx, *n, *k, Rx);

    // create 2D permute index matrix:
    for (int j = 0; j < *k; j++) {
        for (int i = 0; i < *n; i++) { i_perm[j][i] = i; }
    }

    // compute kbcov value
    K_Ball_Covariance(kbcov_stat, Dx, Rx, i_perm, n, k);

    // permutation test
    if (*R > 0) {
        double *permuted_bcov_weight0, *permuted_bcov_weight_prob, *permuted_bcov_weight_hhg;
        permuted_bcov_weight0 = (double *) malloc(*R * sizeof(double));
        permuted_bcov_weight_prob = (double *) malloc(*R * sizeof(double));
        permuted_bcov_weight_hhg = (double *) malloc(*R * sizeof(double));

        int not_parallel = *thread == 1 ? 1 : 0;
        int r;
        if (not_parallel) {
            double bcov_tmp[3];
            for (r = 0; r < *R; r++) {
                // stop permutation if user stop it manually:
                if (pending_interrupt()) {
                    print_stop_message();
                    break;
                }
                resample_matrix(i_perm, k, n);
                // for (int j = 0; j < *k; j++) {
                //     for (int i = 0; i < *n; i++) { Rprintf("%d, ", i_perm[j][i]); }
                //     Rprintf("\n");
                // }
                K_Ball_Covariance(bcov_tmp, Dx, Rx, i_perm, n, k);
                permuted_bcov_weight0[r] = bcov_tmp[0];
                permuted_bcov_weight_prob[r] = bcov_tmp[1];
                permuted_bcov_weight_hhg[r] = bcov_tmp[2];
            }
        } else {
            int ***i_perm_3d_matrix = alloc_3d_int_matrix(*R, *k, *n);
            resample_matrix_3d(i_perm_3d_matrix, i_perm, R, k, n);
#pragma omp parallel
            {
                int i_thread;
                double bcov_tmp[3];
#pragma omp for
                for (i_thread = 0; i_thread < *R; i_thread++) {
                    K_Ball_Covariance(bcov_tmp, Dx, Rx, i_perm_3d_matrix[i_thread], n, k);
                    permuted_bcov_weight0[i_thread] = bcov_tmp[0];
                    permuted_bcov_weight_prob[i_thread] = bcov_tmp[1];
                    permuted_bcov_weight_hhg[i_thread] = bcov_tmp[2];
                }
            }
            free_3d_int_matrix(i_perm_3d_matrix, *R, *k);
            r = *R;
        }

        pvalue[0] = compute_pvalue(kbcov_stat[0], permuted_bcov_weight0, r);
        pvalue[1] = compute_pvalue(kbcov_stat[1], permuted_bcov_weight_prob, r);
        pvalue[2] = compute_pvalue(kbcov_stat[2], permuted_bcov_weight_hhg, r);

        free(permuted_bcov_weight0);
        free(permuted_bcov_weight_prob);
        free(permuted_bcov_weight_hhg);
    }
    free_3d_matrix(Dx, *n, *n);
    free_3d_int_matrix(Rx, *n, *n);
    free_int_matrix(i_perm, *k, *n);
}

/**
 * @todo
 */
void UKBCOV(double *kbcov_stat, double *pvalue, double *x, int *k, int *n, int *R, int *thread) {

}

void kbcov_test(double *kbcov_stat, double *pvalue, double *x, int *k, int *n, int *R, const int *dst,
                int *thread) {
    int not_parallel = *thread == 1 ? 1 : 0;
#ifdef Ball_OMP_H_
    if (not_parallel != 1) {
        omp_set_dynamic(0);
        if (*thread <= 0) {
            omp_set_num_threads(omp_get_num_procs());
        } else {
            omp_set_num_threads(*thread);
        }
    }
#endif

#ifdef R_BUILD
    GetRNGstate();
#else
    srand((unsigned) time(NULL));
#endif

    if (*dst == 1) {
        KBCOV(kbcov_stat, pvalue, x, k, n, R, thread);
    } else {
        UKBCOV(kbcov_stat, pvalue, x, k, n, R, thread);
    }
    
#ifdef R_BUILD
    PutRNGstate();
#endif
}
