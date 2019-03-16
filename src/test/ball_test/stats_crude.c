//
// Created by JinZhu on 2019/3/15.
//


#include "stdio.h"
#include "stdlib.h"
#include "math.h"

/**
 * A crude way to implement Ball Divergence
 * @param bd_stat The computation result for Ball Divergence will be save in here
 * @param Dx Distance matrix for two group
 * @param n1 sample size of first group
 * @param n2 sample size of second group
 */
void Ball_Divergence_Crude(double *bd_stat, double **Dx, int n1, int n2) {
    double A_nm = 0.0, C_nm = 0.0, x_count, y_count;
    for (int i = 0; i < n1; ++i) {
        for (int j = 0; j < n1; ++j) {
            x_count = y_count = 0.0;
            for (int k = 0; k < n1; ++k) {
                if (Dx[i][k] <= Dx[i][j]) {
                    x_count++;
                }
            }
            for (int k = n1; k < (n2 + n1); ++k) {
                if (Dx[i][k] <= Dx[i][j]) {
                    y_count++;
                }
            }
            A_nm += pow((x_count / n1 - y_count / n2), 2);
        }
    }
    A_nm /= pow(n1, 2);

    for (int i = n1; i < (n2 + n1); ++i) {
        for (int j = n1; j < (n2 + n1); ++j) {
            x_count = y_count = 0.0;
            for (int k = 0; k < n1; ++k) {
                if (Dx[i][k] <= Dx[i][j]) {
                    x_count++;
                }
            }
            for (int k = n1; k < (n2 + n1); ++k) {
                if (Dx[i][k] <= Dx[i][j]) {
                    y_count++;
                }
            }
            C_nm += pow(x_count / n1 - y_count / n2, 2);
        }
    }
    C_nm /= pow(n2, 2);

    bd_stat[0] = A_nm + C_nm;
}

/**
 * A crude way to implement Ball Covariance
 * @param bcov_stat The computation result for bcov_stat will be save in here
 * @param Dx Distance matrix for variable x
 * @param Dy Distance matrix for variable y
 * @param n sample size
 */
void Ball_Covariance_Crude(double *bcov_stat, double **Dx, double **Dy, int n) {
    double x_count, y_count, joint_count, bcov_weight0 = 0.0;
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            y_count = x_count = joint_count = 0.0;
            for (int k = 0; k < n; ++k) {
                if (Dx[i][k] <= Dx[i][j] && Dy[i][k] <= Dy[i][j]) {
                    joint_count++;
                }
                if (Dx[i][k] <= Dx[i][j]) {
                    x_count++;
                }
                if (Dy[i][k] <= Dy[i][j]) {
                    y_count++;
                }
            }
            bcov_weight0 += pow((joint_count / (n) - x_count * y_count / pow(n, 2)), 2);
        }
    }
    bcov_weight0 /= pow(n, 2);
    bcov_stat[0] = bcov_weight0;
}

/**
 * A crude way to implement K Ball Covariance
 * @param kbcov_stat The computation result for kbcov_stat will be save in here
 * @param Dx Distance matrix
 * @param n sample size
 * @param k variable number
 */
void K_Ball_Covariance_Crude(double *kbcov_stat, double ***Dx, int n, int k) {
    int i, j, t, var_index;
    double p_all, p_prod, p_inv_prod;
    double *p_k_array;
    double bcov_weight0 = 0.0, bcov_weight_prob = 0.0, bcov_weight_hhg = 0.0, bcov_fixed_ball = 0.0;
    int hhg_include_ball_index;
    double hhg_ball_num = 0.0, pow_n_k = pow(n, k);

    p_k_array = (double *) malloc(k * sizeof(double));

    for (i = 0; i < (n); i++) {
        for (j = 0; j < (n); j++) {
            // init the computation for each Ball
            p_all = n;
            for (var_index = 0; var_index < (k); var_index++) {
                p_k_array[var_index] = 0.0;
            }
            p_prod = 1.0, p_inv_prod = 1.0;
            hhg_include_ball_index = 1;

            // compute the count
            for (t = 0; t < (n); t++) {
                // compute P_{i, j}^{\mu_{1}, ..., \mu_{k}}:
                for (var_index = 0; var_index < (k); var_index++) {
                    if (Dx[i][j][var_index] < Dx[i][t][var_index]) {
                        p_all -= 1.0;
                        break;
                    }
                }
                // compute P_{i, j}^{\mu_{1}}, ..., P_{i, j}^{\mu_{k}}:
                for (var_index = 0; var_index < (k); var_index++) {
                    if (Dx[i][j][var_index] >= Dx[i][t][var_index]) {
                        p_k_array[var_index] += 1.0;
                    }
                }
            }

            for (var_index = 0; var_index < (k); var_index++) {
                p_prod *= p_k_array[var_index];
                if (hhg_include_ball_index == 1) {
                    if (p_k_array[var_index] > 2 && p_k_array[var_index] < n) {
                        hhg_include_ball_index = 1;
                    } else {
                        hhg_include_ball_index = 0;
                    }
                }
                p_inv_prod *= (n - p_k_array[var_index]);
            }

            p_all = p_all / n;
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
    kbcov_stat[0] = bcov_weight0 / (1.0 * (n) * (n));
    kbcov_stat[1] = bcov_weight_prob / (1.0 * (n) * (n));
    kbcov_stat[2] = bcov_weight_hhg / (hhg_ball_num);

    free(p_k_array);
}