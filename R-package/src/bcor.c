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

#include "math.h"
#include "string.h"
#include "stdlib.h"
#include "stdio.h"
#include "utilities.h"
#include "bcor.h"
#include "Ball_omp.h"

#ifdef R_BUILD
#include "R.h"
#endif

/**
 *
 * @param sub_rank
 * @param index_matrix
 * @param cumsum_size : [0, 4, 9]
 * @param label : [0, 0, 0, 0, 1, 1, 1, 1, 1, 2, 2, 2, 2]
 * @param num : 13
 * @param max_k : 2
 */
void Ball_Correlation_KSample_NoTies(double *bcor_stat, double **margin_prop, double **distance_matrix,
                                     int **index_matrix, const int *label, const int *size, int num) {
    /*
     * g_index: indicator for group
     * s_index: indicator for order
     */
    int hhg_ball_num = 0;
    int s_index, g_index;
    double inv_num = 1.0 / num;
    double pxy, px, py, bcov_fixed_ball;
    double bcov_weight0 = 0.0, bcov_weight_prob = 0.0, bcov_weight_hhg = 0.0;
    double bcov_weight0_x = 0.0, bcov_weight_prob_x = 0.0, bcov_weight0_y = 0.0, bcov_weight_prob_y = 0.0;
    for (int i = 0; i < num; ++i) {
        int i_index = label[i];
        int rank_value = 1;
        px = size[i_index] * inv_num;
        for (int j = 0; j < num; ++j) {
            s_index = index_matrix[i][j];
            g_index = label[s_index];
            py = margin_prop[i][s_index];
            if (i_index == g_index) {
                pxy = rank_value * inv_num;

                // compute BCov(X, Y)
                bcov_fixed_ball = pxy - px * py;
                bcov_fixed_ball *= bcov_fixed_ball;
                bcov_weight0 += bcov_fixed_ball;
                bcov_weight_prob += bcov_fixed_ball / (px * py);
                if (px != 1 && py != 1) {
                    bcov_weight_hhg += bcov_fixed_ball / ((px) * (1.0 - px) * (py) * (1.0 - py));
                    hhg_ball_num += 1;
                }
                // compute BCov(X, X)
                bcov_weight0_x += px * px * (1.0 - px) * (1.0 - px);
                bcov_weight_prob_x += (1.0 - px) * (1.0 - px);
                // compute BCov(Y, Y)
                bcov_weight0_y += py * py * (1.0 - py) * (1.0 - py);
                bcov_weight_prob_y += (1.0 - py) * (1.0 - py);

                rank_value++;
            } else {
                // compute BCov(Y, Y)
                bcov_weight0_y += py * py * (1.0 - py) * (1.0 - py);
                bcov_weight_prob_y += (1.0 - py) * (1.0 - py);
            }
        }
    }

    if ((bcov_weight_prob_x * bcov_weight_prob_y) > 0) {
        bcor_stat[0] = bcov_weight0 / (sqrt(bcov_weight0_x * bcov_weight0_y));
        bcor_stat[1] = bcov_weight_prob / (sqrt(bcov_weight_prob_x * bcov_weight_prob_y));
        bcor_stat[2] = bcov_weight_hhg / (hhg_ball_num);
    } else {
        bcor_stat[0] = bcor_stat[1] = bcor_stat[2] = 0;
    }
}

void Ball_Correlation_KSample(double *bcor_stat, double **margin_prop, double **distance_matrix,
                              int **index_matrix, const int *label, const int *size, int num) {
    int hhg_ball_num = 0;
    int s_index, g_index;
    double inv_num = 1.0 / num;
    double pxy, px, py, bcov_fixed_ball;
    double bcov_weight0 = 0.0, bcov_weight_prob = 0.0, bcov_weight_hhg = 0.0;
    double bcov_weight0_x = 0.0, bcov_weight_prob_x = 0.0, bcov_weight0_y = 0.0, bcov_weight_prob_y = 0.0;
    for (int i = 0; i < num; ++i) {
        int i_index = label[i];
        int rank_value = size[i_index];
        double tmp = -1;
        int tmp_rank = 0;
        for (int j = 0; j < num; ++j) {
            s_index = index_matrix[i][num - j - 1];
            g_index = label[s_index];
            py = margin_prop[i][s_index];
            if (i_index == g_index) {
                if (distance_matrix[i][s_index] != tmp) {
                    pxy = rank_value * inv_num;
                    tmp_rank = rank_value;
                } else {
                    pxy = tmp_rank * inv_num;
                }
                tmp = distance_matrix[i][s_index];
                rank_value--;

                // Compute BCov(X, Y)
                px = size[i_index] * inv_num;
                bcov_fixed_ball = pxy - px * py;
                bcov_fixed_ball *= bcov_fixed_ball;
                bcov_weight0 += bcov_fixed_ball;
                bcov_weight_prob += bcov_fixed_ball / (px * py);
                if (px != 1 && py != 1) {
                    bcov_weight_hhg += bcov_fixed_ball / ((px) * (1.0 - px) * (py) * (1.0 - py));
                    hhg_ball_num += 1;
                }
                // compute BCov(X, X)
                bcov_weight_prob_x += (1.0 - px) * (1.0 - px);
                bcov_weight0_x += px * px * (1.0 - px) * (1.0 - px);
                // compute BCov(Y, Y)
                bcov_weight_prob_y += (1.0 - py) * (1.0 - py);
                bcov_weight0_y += py * py * (1.0 - py) * (1.0 - py);
            } else {
                // compute BCov(Y, Y)
                bcov_weight0_y += py * py * (1.0 - py) * (1.0 - py);
                bcov_weight_prob_y += (1.0 - py) * (1.0 - py);
            }
        }
    }

    if ((bcov_weight_prob_x * bcov_weight_prob_y) > 0) {
        bcor_stat[0] = bcov_weight0 / (sqrt(bcov_weight0_x * bcov_weight0_y));
        bcor_stat[1] = bcov_weight_prob / (sqrt(bcov_weight_prob_x * bcov_weight_prob_y));
        bcor_stat[2] = bcov_weight_hhg / (hhg_ball_num);
    } else {
        bcor_stat[0] = bcor_stat[1] = bcor_stat[2] = 0;
    }
}

void Ball_Correlation(double *bcor_stat, const int *n, double **Dx, double **Dy, int **xidx, int **yidx) {
    int i, j, k, pi, src, lastpos, *yrank, *isource, *icount, *xy_index, *xy_temp, **xyidx;
    double pxy, px, py, lastval, *xx_cpy, *yy_cpy;
    double bcov_fixed_ball, bcov_weight0 = 0.0, bcov_weight_prob = 0.0, bcov_weight_hhg = 0.0;
    double bcov_weight0_x = 0.0, bcov_weight_prob_x = 0.0;
    double bcov_weight0_y = 0.0, bcov_weight_prob_y = 0.0;
    double hhg_ball_num = 0.0;

    yrank = (int *) malloc(*n * sizeof(int));
    isource = (int *) malloc(*n * sizeof(int));
    icount = (int *) malloc(*n * sizeof(int));
    xy_index = (int *) malloc(*n * sizeof(int));
    xy_temp = (int *) malloc(*n * sizeof(int));
    xyidx = alloc_int_matrix(*n, *n);
    xx_cpy = (double *) malloc(*n * sizeof(double));
    yy_cpy = (double *) malloc(*n * sizeof(double));

    for (i = 0; i < *n; i++) {
        for (j = 0; j < *n; j++) {
            xyidx[i][j] = j;
        }
    }

    for (i = 0; i < (*n); i++) {
        memcpy(xx_cpy, Dx[i], *n * sizeof(double));
        for (j = 0; j < *n; j++) {
            yy_cpy[j] = Dy[i][j];
        }
        quicksort2(xx_cpy, yy_cpy, xyidx[i], 0, *n - 1);
    }
    free(xx_cpy);
    free(yy_cpy);

    for (i = 0; i < (*n); i++) {
        pi = i;
        lastval = 0;
        lastpos = -1;
        for (j = *n - 1, k = *n - 1; j >= 1; --j, --k) {
            k -= (yidx[pi][k] == pi);
            if (lastpos == -1 || Dy[pi][yidx[pi][k]] != lastval) {
                lastval = Dy[pi][yidx[pi][k]];
                lastpos = j;
            }
            src = yidx[pi][k];
            src -= (src > i);
            yrank[src] = lastpos;
        }

        for (j = 0, k = 0; j < *n - 1; ++j, ++k) {
            k += (xyidx[i][k] == i);
            src = xyidx[i][k];
            src -= (src > i);
            xy_index[j] = yrank[src];
            isource[j] = j;
            icount[j] = 0;
            xy_temp[j] = xy_index[j];
        }
        Inversions(xy_temp, isource, icount, *n - 1, *n);
        lastval = 0;
        lastpos = -1;

        for (j = *n - 2, k = *n - 1; j >= 0; --j, --k) {
            k -= (xidx[i][k] == i);
            if (lastpos == -1 || Dx[i][xidx[i][k]] != lastval) {
                lastval = Dx[i][xidx[i][k]];
                lastpos = j;
            }

            pxy = lastpos - icount[j] + 2;
            px = lastpos + 2;
            py = xy_index[j] + 1;
            px /= (*n);
            py /= (*n);
            pxy /= (*n);
            // compute BCov(X, Y)
            bcov_fixed_ball = pxy - px * py;
            bcov_fixed_ball *= bcov_fixed_ball;
            bcov_weight0 += bcov_fixed_ball;
            bcov_weight_prob += bcov_fixed_ball / (px * py);
            if (px != 1 && py != 1) {
                bcov_weight_hhg += bcov_fixed_ball / ((px) * (1.0 - px) * (py) * (1.0 - py));
                hhg_ball_num += 1;
            }
            // compute BCov(X, X)
            bcov_weight_prob_x += (1.0 - px) * (1.0 - px);
            bcov_weight0_x += px * px * (1.0 - px) * (1.0 - px);
            // compute BCov(Y, Y)
            bcov_weight_prob_y += (1.0 - py) * (1.0 - py);
            bcov_weight0_y += py * py * (1.0 - py) * (1.0 - py);
        }
        pxy = 0;
        px = 0;
        py = 0;
        for (j = 0; j < *n; j++) {
            if (Dx[i][xidx[i][j]] == 0) {
                px += 1;
                if (Dy[pi][xidx[i][j]] == 0) {
                    pxy += 1;
                    py += 1;
                }
            } else if (Dy[pi][xidx[i][j]] == 0) {
                py += 1;
            }
        }
        px /= (*n);
        py /= (*n);
        pxy /= (*n);
        bcov_fixed_ball = pow(pxy - px * py, 2);
        bcov_weight0 += bcov_fixed_ball;
        bcov_weight_prob += bcov_fixed_ball / (px * py);
        if (px != 1 && py != 1) {
            bcov_weight_hhg += bcov_fixed_ball / ((px) * (1.0 - px) * (py) * (1.0 - py));
            hhg_ball_num += 1;
        }
        // compute BCov(X, X)
        bcov_weight_prob_x += (1.0 - px) * (1.0 - px);
        bcov_weight0_x += px * px * (1.0 - px) * (1.0 - px);
        // compute BCov(Y, Y)
        bcov_weight_prob_y += (1.0 - py) * (1.0 - py);
        bcov_weight0_y += py * py * (1.0 - py) * (1.0 - py);
    }
    if ((bcov_weight_prob_x * bcov_weight_prob_y) > 0) {
        bcor_stat[0] = bcov_weight0 / (sqrt(bcov_weight0_x * bcov_weight0_y));
        bcor_stat[1] = bcov_weight_prob / (sqrt(bcov_weight_prob_x * bcov_weight_prob_y));
        bcor_stat[2] = bcov_weight_hhg / (hhg_ball_num);
    } else {
        bcor_stat[0] = bcor_stat[1] = bcor_stat[2] = 0;
    }

    free_int_matrix(xyidx, *n, *n);
    free(isource);
    free(icount);
    free(xy_index);
    free(yrank);
    free(xy_temp);
}

void Ball_Correlation_NoTies(double *bcor_stat, const int *n, int **y_within_ball, int **xidx, double **Dy) {
    double px, py, pxy, n_prop = 1.0 / *n;
    double bcov_weight0_x = 0.0, bcov_weight_prob_x = 0.0;
    double bcov_weight0_y = 0.0, bcov_weight_prob_y = 0.0;
    int sorted_j, *inv_count = malloc(*n * sizeof(int));
    int *y_count_vec = malloc(*n * sizeof(int));
    double *dy_vec = malloc(*n * sizeof(double)), *dy_vec_tmp = malloc(*n * sizeof(double));
    double bcov_fixed_ball, bcov_weight0 = 0.0, bcov_weight_prob = 0.0, bcov_weight_hhg = 0.0, hhg_ball_num = 0.0;

    for (int i = 0; i < *n; ++i) {
        for (int j = 0; j < *n; ++j) {
            y_count_vec[j] = y_within_ball[i][j];
            dy_vec_tmp[j] = Dy[i][j];
            inv_count[j] = 0;
        }
        for (int j = 0; j < *n; ++j) {
            dy_vec[j] = dy_vec_tmp[xidx[i][j]];
        }
        // Count of the number after self:
        count_smaller_number_after_self_solution(dy_vec, inv_count, *n);
        // Compute Ball Covariance:
        for (int j = 0; j < *n; ++j) {
            sorted_j = xidx[i][j];
            px = (j + 1) * n_prop;
            py = y_count_vec[sorted_j];
            pxy = (py - inv_count[j]) * n_prop;
            py *= n_prop;
            bcov_fixed_ball = pxy - px * py;
            bcov_fixed_ball *= bcov_fixed_ball;
            bcov_weight0 += bcov_fixed_ball;
            bcov_weight_prob += bcov_fixed_ball / (px * py);
            if (px != 1 && py != 1) {
                bcov_weight_hhg += bcov_fixed_ball / ((px) * (1.0 - px) * (py) * (1.0 - py));
                hhg_ball_num += 1;
            }
            // compute BCov(X, X)
            bcov_weight_prob_x += (1.0 - px) * (1.0 - px);
            bcov_weight0_x += px * px * (1.0 - px) * (1.0 - px);
            // compute BCov(Y, Y)
            bcov_weight_prob_y += (1.0 - py) * (1.0 - py);
            bcov_weight0_y += py * py * (1.0 - py) * (1.0 - py);
        }
    }
    if ((bcov_weight_prob_x * bcov_weight_prob_y) > 0) {
        bcor_stat[0] = bcov_weight0 / (sqrt(bcov_weight0_x * bcov_weight0_y));
        bcor_stat[1] = bcov_weight_prob / (sqrt(bcov_weight_prob_x * bcov_weight_prob_y));
        bcor_stat[2] = bcov_weight_hhg / (hhg_ball_num);
    } else {
        bcor_stat[0] = bcor_stat[1] = bcor_stat[2] = 0;
    }

    free(inv_count);
    free(dy_vec);
    free(dy_vec_tmp);
    free(y_count_vec);
}

void U_Ball_Correlation(double *bcor_stat, int *n, double *x, int *yrank, int **lowyidx, int **higyidx) {
    int **Rank, **lowxidx, **higxidx, *xidx, *xrank;
    xidx = (int *) malloc(*n * sizeof(int));
    xrank = (int *) malloc(*n * sizeof(int));
    Rank = alloc_int_matrix((*n) + 1, (*n) + 1);
    lowxidx = alloc_int_matrix(*n, *n);
    higxidx = alloc_int_matrix(*n, *n);

    for (int k = 0; k < *n; k++) {
        xidx[k] = k;
    }
    quicksort(x, xidx, 0, *n - 1);
    ranksort(n, xrank, x, xidx);
    createidx(n, xidx, x, lowxidx, higxidx);
    initRank_bcor(*n, Rank, xrank, yrank);
    free(xrank);
    free(xidx);

    int i, j, pi, pj;
    double px, py, pxy;
    double bcov_fixed_ball, bcov_weight0 = 0.0, bcov_weight_prob = 0.0, bcov_weight_hhg = 0.0;
    double bcov_weight0_x = 0.0, bcov_weight_prob_x = 0.0;
    double bcov_weight0_y = 0.0, bcov_weight_prob_y = 0.0;
    double hhg_ball_num = 0.0;
    for (i = 0; i < *n; i++) {
        for (j = 0; j < *n; j++) {
            pi = i;
            pj = j;
            px = higxidx[i][j] - lowxidx[i][j] + 1;
            py = higyidx[pi][pj] - lowyidx[pi][pj] + 1;
            pxy = Rank[higxidx[i][j]][higyidx[pi][pj]] + Rank[lowxidx[i][j] - 1][lowyidx[pi][pj] - 1];
            pxy -= (Rank[higxidx[i][j]][lowyidx[pi][pj] - 1] + Rank[lowxidx[i][j] - 1][higyidx[pi][pj]]);

            pxy /= (*n);
            px /= (*n);
            py /= (*n);
            bcov_fixed_ball = pow(pxy - px * py, 2);
            bcov_weight0 += bcov_fixed_ball;
            bcov_weight_prob += bcov_fixed_ball / (px * py);
            if (px != 1 && py != 1) {
                bcov_weight_hhg += bcov_fixed_ball / ((px) * (1.0 - px) * (py) * (1.0 - py));
                hhg_ball_num += 1;
            }
            bcov_weight_prob_x += (1.0 - px) * (1.0 - px);
            bcov_weight_prob_y += (1.0 - py) * (1.0 - py);
            bcov_weight0_x += px * px * (1.0 - px) * (1.0 - px);
            bcov_weight0_y += py * py * (1.0 - py) * (1.0 - py);
        }
    }
    if ((bcov_weight_prob_x * bcov_weight_prob_y) > 0) {
        bcor_stat[0] = bcov_weight0 / (sqrt(bcov_weight0_x * bcov_weight0_y));
        bcor_stat[1] = bcov_weight_prob / (sqrt(bcov_weight_prob_x * bcov_weight_prob_y));
        bcor_stat[2] = bcov_weight_hhg / (hhg_ball_num);
    } else {
        bcor_stat[0] = bcor_stat[1] = bcor_stat[2] = 0;
    }

    free_int_matrix(Rank, (*n) + 1, (*n) + 1);
    free_int_matrix(lowxidx, *n, *n);
    free_int_matrix(higxidx, *n, *n);
}

void _bcor_test(double *bcor_stat, double *y, double *x, int *x_number, int *feature_number,
                int *n, int *p) {
    int i, j, y_ties = 0;
    double **Dy = alloc_matrix(*n, *n);
    if (*p == 0) {
        distance2matrix(y, Dy, *n);
    } else {
        Euclidean_distance(y, Dy, *n, *p);
    }

    int **yidx = alloc_int_matrix(*n, *n);
    for (i = 0; i < *n; i++) {
        for (j = 0; j < *n; j++) {
            yidx[i][j] = j;
        }
    }

    // sort distance matrix y:
    double *y_cpy = (double *) malloc(*n * sizeof(double));
    int **y_within_ball = alloc_int_matrix(*n, *n);
    for (i = 0; i < (*n); i++) {
        memcpy(y_cpy, Dy[i], *n * sizeof(double));
        quicksort(y_cpy, yidx[i], 0, *n - 1);
        if (!y_ties) {
            for (j = 1; j < *n; ++j) {
                if (y_cpy[j] == y_cpy[j - 1]) {
                    y_ties = 1;
                }
            }
        }
    }
    free(y_cpy);
    if (y_ties) {
        for (i = 0; i < *n; ++i) {
            quick_rank_max_with_index(Dy[i], yidx[i], y_within_ball[i], *n);
        }
    }

#pragma omp parallel
    {
        int i_feature_thread, i_thread, j_thread, k_thread, stop_index_thread, x_size_thread, x_ties = 0;
        double bcorsis_stat_tmp[3];
        int **y_within_ball_thread, **xidx;
        double **Dx, *x_cpy, *x_thread;
        // main loop for calculate bcor statistic
#pragma omp for
        for (i_feature_thread = 0; i_feature_thread < *feature_number; i_feature_thread++) {
            // extract value to x_thread
            k_thread = 0;
            i_thread = (i_feature_thread) * (*n);
            x_size_thread = x_number[i_feature_thread] * (*n);
            stop_index_thread = (i_feature_thread) * (*n) + x_size_thread;
            x_thread = (double *) malloc(x_size_thread * sizeof(double));
            while (i_thread < stop_index_thread) {
                x_thread[k_thread] = x[i_thread];
                k_thread++;
                i_thread++;
            }

            Dx = alloc_matrix(*n, *n);
            x_cpy = (double *) malloc(*n * sizeof(double));
            xidx = alloc_int_matrix(*n, *n);
            Euclidean_distance(x_thread, Dx, *n, x_number[i_feature_thread]);
            for (i_thread = 0; i_thread < *n; i_thread++) {
                for (j_thread = 0; j_thread < *n; j_thread++) {
                    xidx[i_thread][j_thread] = j_thread;
                }
            }
            for (i_thread = 0; i_thread < (*n); i_thread++) {
                memcpy(x_cpy, Dx[i_thread], *n * sizeof(double));
                quicksort(x_cpy, xidx[i_thread], 0, *n - 1);
                if (!x_ties) {
                    for (j_thread = 1; j_thread < *n; ++j_thread) {
                        if (x_cpy[j_thread] == x_cpy[j_thread - 1]) {
                            x_ties = 1;
                        }
                    }
                }
            }
            free(x_thread);
            free(x_cpy);

            if (y_ties) {
                if (x_ties) {
                    Ball_Correlation(bcorsis_stat_tmp, n, Dx, Dy, xidx, yidx);
                } else {
                    Ball_Correlation_NoTies(bcorsis_stat_tmp, n, y_within_ball, xidx, Dy);
                }
            } else {
                y_within_ball_thread = alloc_int_matrix(*n, *n);
                for (i_thread = 0; i_thread < *n; ++i_thread) {
                    quick_rank_max_with_index(Dx[i_thread], xidx[i_thread], y_within_ball_thread[i_thread], *n);
                }
                Ball_Correlation_NoTies(bcorsis_stat_tmp, n, y_within_ball_thread, yidx, Dx);
                free_int_matrix(y_within_ball_thread, *n, *n);
            }
            bcor_stat[3 * i_feature_thread + 0] = bcorsis_stat_tmp[0];
            bcor_stat[3 * i_feature_thread + 1] = bcorsis_stat_tmp[1];
            bcor_stat[3 * i_feature_thread + 2] = bcorsis_stat_tmp[2];
            free_matrix(Dx, *n, *n);
            free_int_matrix(xidx, *n, *n);
        }
    }

    free_matrix(Dy, *n, *n);
    free_int_matrix(yidx, *n, *n);
    free_int_matrix(y_within_ball, *n, *n);
}

void _u_bcor_test(double *bcor_stat, double *y, double *x, int *f_number, int *n) {
    int i, *yidx, *yrank, **lowyidx, **higyidx;

    yidx = (int *) malloc(*n * sizeof(int));
    yrank = (int *) malloc(*n * sizeof(int));
    lowyidx = alloc_int_matrix(*n, *n);
    higyidx = alloc_int_matrix(*n, *n);

    for (i = 0; i < *n; i++) {
        yidx[i] = i;
    }
    quicksort(y, yidx, 0, *n - 1);
    ranksort(n, yrank, y, yidx);
    createidx(n, yidx, y, lowyidx, higyidx);
    free(yidx);

#pragma omp parallel
    {
        int f_thread, i_thread, s_thread, start_index;
        double *x_thread;
        double bcorsis_stat_tmp[3];

        x_thread = (double *) malloc((*n) * sizeof(double));

        // main loop for calculate bcor statistic
#pragma omp for
        for (f_thread = 0; f_thread < *f_number; f_thread++) {
            // extract value to x_thread
            start_index = f_thread * (*n);
            for (i_thread = 0; i_thread < (*n); i_thread++) {
                x_thread[i_thread] = x[start_index + i_thread];
            }

            U_Ball_Correlation(bcorsis_stat_tmp, n, x_thread, yrank, lowyidx, higyidx);
            for (s_thread = 0; s_thread < 3; s_thread++) {
                bcor_stat[3 * f_thread + s_thread] = bcorsis_stat_tmp[s_thread];
            }
        }
        free(x_thread);
    }

    free(yrank);
    free_int_matrix(lowyidx, *n, *n);
    free_int_matrix(higyidx, *n, *n);
}

void _bcor_stat(double *bcor_stat, double *y, double *x, const int *n) {
    double **Dx, **Dy, *x_cpy, *y_cpy;
    int **xidx, **yidx;
    Dx = alloc_matrix(*n, *n);
    Dy = alloc_matrix(*n, *n);
    xidx = alloc_int_matrix(*n, *n);
    yidx = alloc_int_matrix(*n, *n);

    x_cpy = (double *) malloc((*n) * sizeof(double));
    y_cpy = (double *) malloc((*n) * sizeof(double));
    distance2matrix(x, Dx, *n);
    distance2matrix(y, Dy, *n);

    int s, t;
    for (s = 0; s < *n; s++) {
        for (t = 0; t < *n; t++) {
            xidx[s][t] = t;
            yidx[s][t] = t;
        }
    }

    for (s = 0; s < (*n); s++) {
        // copy site to x_cpy and y_cpy
        memcpy(x_cpy, Dx[s], *n * sizeof(double));
        memcpy(y_cpy, Dy[s], *n * sizeof(double));
        quicksort(x_cpy, xidx[s], 0, *n - 1);
        quicksort(y_cpy, yidx[s], 0, *n - 1);
    }
    free(x_cpy);
    free(y_cpy);

    int i, j, k, pi, src, lastpos, *yrank, *isource, *icount, *xy_index, *xy_temp, **xyidx;
    double pxy, px, py, lastval, *xx_cpy, *yy_cpy;
    double bcov_fixed_ball, bcov_weight0 = 0.0, bcov_weight_prob = 0.0, bcov_weight_hhg = 0.0;
    double bcov_weight0_x = 0.0, bcov_weight_prob_x = 0.0;
    double bcov_weight0_y = 0.0, bcov_weight_prob_y = 0.0;
    double hhg_ball_num = 0.0;

    yrank = (int *) malloc(*n * sizeof(int));
    isource = (int *) malloc(*n * sizeof(int));
    icount = (int *) malloc(*n * sizeof(int));
    xy_index = (int *) malloc(*n * sizeof(int));
    xy_temp = (int *) malloc(*n * sizeof(int));
    xyidx = alloc_int_matrix(*n, *n);
    xx_cpy = (double *) malloc(*n * sizeof(double));
    yy_cpy = (double *) malloc(*n * sizeof(double));

    for (i = 0; i < *n; i++)
        for (j = 0; j < *n; j++)
            xyidx[i][j] = j;

    for (i = 0; i < (*n); i++) {
        memcpy(xx_cpy, Dx[i], *n * sizeof(double));
        for (j = 0; j < *n; j++)
            yy_cpy[j] = Dy[i][j];
        quicksort2(xx_cpy, yy_cpy, xyidx[i], 0, *n - 1);
    }
    free(xx_cpy);
    free(yy_cpy);

    for (i = 0; i < (*n); i++) {
        pi = i;
        lastval = 0;
        lastpos = -1;
        for (j = *n - 1, k = *n - 1; j >= 1; --j, --k) {
            k -= (yidx[pi][k] == pi);
            if (lastpos == -1 || Dy[pi][yidx[pi][k]] != lastval) {
                lastval = Dy[pi][yidx[pi][k]];
                lastpos = j;
            }
            src = yidx[pi][k];
            src -= (src > i);
            yrank[src] = lastpos;
        }

        for (j = 0, k = 0; j < *n - 1; ++j, ++k) {
            k += (xyidx[i][k] == i);
            src = xyidx[i][k]; // NOTE: k may be different than from the line above
            src -= (src > i);
            xy_index[j] = yrank[src];
            isource[j] = j;
            icount[j] = 0;
            xy_temp[j] = xy_index[j];
        }
        Inversions(xy_temp, isource, icount, *n - 1, *n);
        lastval = 0;
        lastpos = -1;

        for (j = *n - 2, k = *n - 1; j >= 0; --j, --k) {
            k -= (xidx[i][k] == i);
            if (lastpos == -1 || Dx[i][xidx[i][k]] != lastval) {
                lastval = Dx[i][xidx[i][k]];
                lastpos = j;
            }

            pxy = lastpos - icount[j] + 2;
            px = lastpos + 2;
            py = xy_index[j] + 1;
            px /= (*n);
            py /= (*n);
            pxy /= (*n);
            // compute BCov(X, Y)
            bcov_fixed_ball = pow(pxy - px * py, 2);
            bcov_weight0 += bcov_fixed_ball;
            bcov_weight_prob += bcov_fixed_ball / (px * py);
            if (px != 1 && py != 1) {
                bcov_weight_hhg += bcov_fixed_ball / ((px) * (1.0 - px) * (py) * (1.0 - py));
                hhg_ball_num += 1;
            }
            // compute BCov(X, X)
            bcov_weight_prob_x += (1.0 - px) * (1.0 - px);
            bcov_weight0_x += px * px * (1.0 - px) * (1.0 - px);

            // compute BCov(Y, Y)
            bcov_weight_prob_y += (1.0 - py) * (1.0 - py);
            bcov_weight0_y += py * py * (1.0 - py) * (1.0 - py);
        }
        pxy = 0;
        px = 0;
        py = 0;
        for (j = 0; j < *n; j++) {
            if (Dx[i][xidx[i][j]] == 0) {
                px += 1;
                if (Dy[pi][xidx[i][j]] == 0) {
                    pxy += 1;
                    py += 1;
                }
            } else if (Dy[pi][xidx[i][j]] == 0)
                py += 1;
        }
        px /= (*n);
        py /= (*n);
        pxy /= (*n);
        // compute BCov(X, Y)
        bcov_fixed_ball = pow(pxy - px * py, 2);
        bcov_weight0 += bcov_fixed_ball;
        bcov_weight_prob += bcov_fixed_ball / (px * py);
        if (px != 1 && py != 1) {
            bcov_weight_hhg += bcov_fixed_ball / ((px) * (1.0 - px) * (py) * (1.0 - py));
            hhg_ball_num += 1;
        }
        // compute BCov(X, X)
        bcov_weight_prob_x += (1.0 - px) * (1.0 - px);
        bcov_weight0_x += px * px * (1.0 - px) * (1.0 - px);

        // compute BCov(Y, Y)
        bcov_weight_prob_y += (1.0 - py) * (1.0 - py);
        bcov_weight0_y += py * py * (1.0 - py) * (1.0 - py);
    }
    if ((bcov_weight_prob_x * bcov_weight_prob_y) > 0) {
        bcor_stat[0] = bcov_weight0 / (sqrt(bcov_weight0_x * bcov_weight0_y));
        bcor_stat[1] = bcov_weight_prob / (sqrt(bcov_weight_prob_x * bcov_weight_prob_y));
        bcor_stat[2] = bcov_weight_hhg / (hhg_ball_num);
    } else {
        bcor_stat[0] = bcor_stat[1] = bcor_stat[2] = 0;
    }

    // free memory
    free_int_matrix(xidx, *n, *n);
    free_int_matrix(yidx, *n, *n);
    free_matrix(Dx, *n, *n);
    free_matrix(Dy, *n, *n);
    free(isource);
    free(icount);
    free(xy_index);
    free(yrank);
    free(xy_temp);
    free_int_matrix(xyidx, *n, *n);
}

void _bcor_k_sample(double *bcor_stat, double *xy, const double *snp, const int *p,
                    const int *n, const int *y_p) {
    int *snp_vector = (int *) malloc(*n * sizeof(int)), *snp_index = (int *) malloc(*n * sizeof(int));
    int **snp_matrix = alloc_int_matrix(*p, *n);
    int **size_list = (int **) malloc(*p * sizeof(int *));
    int *cumsum_size_list;
    int s = 0;
    for (int i = 0; i < *p; ++i) {
        for (int j = 0; j < *n; ++j) {
            snp_vector[j] = (int) snp[s++];
            snp_index[j] = j;
        }
        quicksort_int(snp_vector, snp_index, 0, *n - 1);
        int u = 0;
        snp_matrix[i][snp_index[0]] = u;
        for (int j = 1; j < *n; ++j) {
            if (snp_vector[j] != snp_vector[j - 1]) {
                u++;
            }
            snp_matrix[i][snp_index[j]] = u;
        }
        u++;
        size_list[i] = (int *) malloc(u * sizeof(int));
        cumsum_size_list = (int *) malloc(u * sizeof(int));
        cumsum_size_list[0] = 0;
        int t = 1;
        for (int j = 1; j < *n; ++j) {
            if (snp_vector[j] != snp_vector[j - 1]) {
                cumsum_size_list[t] = j;
                size_list[i][t - 1] = j - cumsum_size_list[t - 1];
                t++;
            }
        }
        size_list[i][u - 1] = *n - cumsum_size_list[u - 1];
        free(cumsum_size_list);
    }
    free(snp_vector);
    free(snp_index);

    // pre-computing marginal ball
    double **distance_matrix = alloc_matrix(*n, *n);
    if (*y_p == 0) {
        distance2matrix(xy, distance_matrix, *n);
    } else {
        Euclidean_distance(xy, distance_matrix, *n, *y_p);
    }

    int **index_matrix = alloc_int_matrix(*n, *n);
    for (int i = 0; i < *n; i++) {
        for (int j = 0; j < *n; j++) {
            index_matrix[i][j] = j;
        }
    }
    double *distance_matrix_copy = (double *) malloc(*n * sizeof(double));
    int *marginal_rank = (int *) malloc(*n * sizeof(int));
    double **marginal_prop = alloc_matrix(*n, *n);
    int ties = 0;
    double inv_n = 1.0 / *n;
    for (int i = 0; i < *n; i++) {
        memcpy(distance_matrix_copy, distance_matrix[i], *n * sizeof(double));
        quicksort(distance_matrix_copy, index_matrix[i], 0, *n - 1);
        quick_rank_max_with_index(distance_matrix[i], index_matrix[i], marginal_rank, *n);
        for (int j = 0; j < *n; ++j) {
            marginal_prop[i][j] = marginal_rank[j] * inv_n;
        }
        if (!ties) {
            for (int j = 1; j < *n; ++j) {
                if (distance_matrix_copy[j] == distance_matrix_copy[j - 1]) {
                    ties = 1;
                    break;
                }
            }
        }
    }
    free(distance_matrix_copy);
    free(marginal_rank);

    void (*ball_correlation_value)(double *, double **, double **, int **, const int *, const int *, int);
    if (ties) {
        ball_correlation_value = &Ball_Correlation_KSample;
    } else {
        ball_correlation_value = &Ball_Correlation_KSample_NoTies;
    }

    // compute statistic:
    double **bcor_stat_array = alloc_matrix(*p, 3);
#pragma omp parallel
    {
        int i_thread;
#pragma omp for
        for (i_thread = 0; i_thread < *p; ++i_thread) {
            ball_correlation_value(bcor_stat_array[i_thread], marginal_prop, distance_matrix, index_matrix,
                                   snp_matrix[i_thread], size_list[i_thread], *n);
        }
    }

    free_matrix(marginal_prop, *n, *n);
    free_matrix(distance_matrix, *n, *n);
    free_int_matrix(index_matrix, *n, *n);
    free_int_matrix(snp_matrix, *p, *n);
    for (int i = 0; i < *p; ++i) {
        free(size_list[i]);
    }
    free(size_list);

    // return result:
    for (int i = 0; i < *p; ++i) {
        bcor_stat[3 * i + 0] = bcor_stat_array[i][0];
        bcor_stat[3 * i + 1] = bcor_stat_array[i][1];
        bcor_stat[3 * i + 2] = bcor_stat_array[i][2];
    }
    free_matrix(bcor_stat_array, *p, 3);
}

/**
 * R API function
 * @param bcorsis_stat : bcor statistics or p-value for screening
 * @param y : response
 * @param x : covariate
 * @param x_number : if x_number = [1, 2, 2, 1, 1], then x[:, 1], x[:, 2:3], x[:, 4:5], x[:, 6], x[:, 7] will be used to compute ball correlation
 * @param f_number : the number of covariate
 * @param size : sample size of each group
 * @param n : total sample size
 * @param p : dimensionality of response variable
 * @param k : group number
 * @param dst_y : whether y should be recompute as distance
 * @param dst_x : whether x should be recompute as distance
 * @param nthread : control the number threads used to compute statistics
 */
void bcor_test(double *bcorsis_stat, double *y, double *x, int *x_number, int *f_number,
               int *n, int *p, int *k, int *dst_y, int *dst_x, int *nthread) {
#ifdef Ball_OMP_H_
    omp_set_dynamic(0);
    if (*nthread <= 0) {
        omp_set_num_threads(omp_get_num_procs());
    } else {
        omp_set_num_threads(*nthread);
    }
#endif
    if (*k <= 1) {
        if (*dst_y == 1) {
            if (*dst_x == 0) {
                _bcor_test(bcorsis_stat, y, x, x_number, f_number, n, p);
            } else {
                _bcor_stat(bcorsis_stat, y, x, n);
            }
        } else {
            _u_bcor_test(bcorsis_stat, y, x, f_number, n);
        }
    } else {
        _bcor_k_sample(bcorsis_stat, y, x, f_number, n, p);
    }

}
