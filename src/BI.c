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
#include "stdlib.h"
#include "string.h"
#include "stdio.h"
#include "utilize_R.h"
#include "utilities.h"
#include "BI.h"
#include "Ball_omp.h"

/**
 * @param y_within_ball : each row is the rank vector of this row (allow ties in each row)
 * @param xidx : the order index of each row (not allow ties in each row)
 * @param Dy : distance matrix (not allow ties)
 * @param i_perm : a permutation of {0, 1, ..., n-1}
 */
void Ball_Information_NoTies(double *bcov_stat, const int *n, int **y_within_ball,
                             int **xidx, double **Dy, const int *i_perm) {
    double px, py, pxy, n_prop = 1.0 / *n;
    int sorted_j, *inv_count = malloc(*n * sizeof(int));
    int *y_count_vec = malloc(*n * sizeof(int));
    double *dy_vec = malloc(*n * sizeof(double)), *dy_vec_tmp = malloc(*n * sizeof(double));
    double bcov_fixed_ball, bcov_weight0 = 0.0, bcov_weight_prob = 0.0, bcov_weight_hhg = 0.0, hhg_ball_num = 0.0;

    for (int i = 0; i < *n; ++i) {
        int ith_perm = i_perm[i];
        for (int j = 0; j < *n; ++j) {
            y_count_vec[j] = y_within_ball[ith_perm][i_perm[j]];
            dy_vec_tmp[j] = Dy[ith_perm][i_perm[j]];
            inv_count[j] = 0;
        }
        for (int j = 0; j < *n; ++j) {
            dy_vec[j] = dy_vec_tmp[xidx[i][j]];
        }
        // Count of the number after self:
        count_smaller_number_after_self_solution(dy_vec, inv_count, *n);
        // Compute Ball Covariance:
        for (int j = 0; j < *n; ++j) {
            sorted_j = xidx[i][j];   // if xidx have ties, this step is wrong.
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
        }
    }
    bcov_stat[0] = bcov_weight0 / (1.0 * (*n) * (*n));
    bcov_stat[1] = bcov_weight_prob / (1.0 * (*n) * (*n));
    bcov_stat[2] = hhg_ball_num > 0 ? (bcov_weight_hhg / hhg_ball_num) : 0.0;
    free(inv_count);
    free(dy_vec);
    free(dy_vec_tmp);
    free(y_count_vec);
}

void Ball_Information2(double *bcov_stat, const int *n, double **Dx, double **Dy, int **xrank, int **yrank,
                       const int *i_perm, int *i_perm_inv) {
    int i, j, *xyidx, *perm_y_rank;
    double pxy, px, py, *xx_cpy, *yy_cpy, *yy_cpy1;
    double bcov_fixed_ball, bcov_weight0 = 0.0, bcov_weight_prob = 0.0, bcov_weight_hhg = 0.0, hhg_ball_num = 0.0;
    double n_prob = 1.0 / (*n);

    xx_cpy = (double *) malloc(*n * sizeof(double));
    yy_cpy = (double *) malloc(*n * sizeof(double));
    yy_cpy1 = (double *) malloc(*n * sizeof(double));
    perm_y_rank = (int *) malloc(*n * sizeof(int));
    xyidx = (int *) malloc(*n * sizeof(int));

    int *right_smaller = (int *) malloc(*n * sizeof(int));
    for (i = 0; i < (*n); i++) {
        memcpy(xx_cpy, Dx[i], *n * sizeof(double));
        for (j = 0; j < *n; j++) {
            yy_cpy[j] = Dy[i_perm[i]][i_perm[j]];
            perm_y_rank[j] = yrank[i_perm[i]][i_perm[j]];
        }
        for (int k = 0; k < *n; k++) {
            xyidx[k] = k;
            right_smaller[k] = 0;
        }
        quicksort3(xx_cpy, yy_cpy, xyidx, 0, *n - 1);
        memcpy(yy_cpy1, yy_cpy, *n * sizeof(double));
        count_smaller_number_after_self_solution(yy_cpy1, right_smaller, *n);

        double tmp_x = -1, tmp_y = -1;
        int smaller, tmp_smaller = right_smaller[*n - 1];
        for (int k = (*n) - 1; k >= 0; k--) {
            smaller = right_smaller[k];
            if (xx_cpy[k] != tmp_x || (xx_cpy[k] == tmp_x && yy_cpy[k] != tmp_y)) {
                tmp_x = xx_cpy[k];
                tmp_y = yy_cpy[k];
                tmp_smaller = smaller;
            } else {
                right_smaller[k] = tmp_smaller;
            }
        }

        for (int k = 0; k < (*n); k++) {
            px = xrank[i][xyidx[k]] * n_prob;
            py = perm_y_rank[xyidx[k]] * n_prob;
            pxy = (perm_y_rank[xyidx[k]] - right_smaller[k]) * n_prob;
            bcov_fixed_ball = pxy - px * py;
            bcov_fixed_ball *= bcov_fixed_ball;
            bcov_weight0 += bcov_fixed_ball;
            bcov_weight_prob += bcov_fixed_ball / (px * py);
            if (xrank[i][xyidx[k]] != (*n) && perm_y_rank[xyidx[k]] != (*n)) {
                bcov_weight_hhg += bcov_fixed_ball / ((px) * (1.0 - px) * (py) * (1.0 - py));
                hhg_ball_num += 1;
            }
        }
    }

    bcov_stat[0] = bcov_weight0 / (1.0 * (*n) * (*n));
    bcov_stat[1] = bcov_weight_prob / (1.0 * (*n) * (*n));
    bcov_stat[2] = hhg_ball_num > 0 ? (bcov_weight_hhg / hhg_ball_num) : 0.0;

    free(right_smaller);
    free(xx_cpy);
    free(yy_cpy);
    free(yy_cpy1);
    free(perm_y_rank);
    free(xyidx);
}

void Ball_Information(double *bcov_stat, int *n, double **Dx, double **Dy, int **xidx, int **yidx, int *i_perm,
                      int *i_perm_inv) {
    int i, j, k, pi, src, lastpos, *yrank, *isource, *icount, *xy_index, *xy_temp, **xyidx;
    double pxy, px, py, lastval, *xx_cpy, *yy_cpy;
    double bcov_fixed_ball, bcov_weight0 = 0.0, bcov_weight_prob = 0.0, bcov_weight_hhg = 0.0, hhg_ball_num = 0.0;

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
            yy_cpy[j] = Dy[i_perm[i]][i_perm[j]];
        }
        quicksort2(xx_cpy, yy_cpy, xyidx[i], 0, *n - 1);
    }
    free(xx_cpy);
    free(yy_cpy);

    for (i = 0; i < (*n); i++) {
        pi = i_perm[i];
        lastval = 0;
        lastpos = -1;
        for (j = *n - 1, k = *n - 1; j >= 1; --j, --k) {
            k -= (yidx[pi][k] == pi);
            if (lastpos == -1 || Dy[pi][yidx[pi][k]] != lastval) {
                lastval = Dy[pi][yidx[pi][k]];
                lastpos = j;
            }
            src = i_perm_inv[yidx[pi][k]];
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
            bcov_fixed_ball = pxy - px * py;
            bcov_fixed_ball *= bcov_fixed_ball;
            bcov_weight0 += bcov_fixed_ball;
            bcov_weight_prob += bcov_fixed_ball / (px * py);
            if (px != 1 && py != 1) {
                bcov_weight_hhg += bcov_fixed_ball / ((px) * (1.0 - px) * (py) * (1.0 - py));
                hhg_ball_num += 1;
            }
        }
        pxy = 0;
        px = 0;
        py = 0;
        for (j = 0; j < *n; j++) {
            if (Dx[i][xidx[i][j]] == 0) {
                px += 1;
                if (Dy[pi][i_perm[xidx[i][j]]] == 0) {
                    pxy += 1;
                    py += 1;
                }
            } else if (Dy[pi][i_perm[xidx[i][j]]] == 0) {
                py += 1;
            }
        }
        px /= (*n);
        py /= (*n);
        pxy /= (*n);
        bcov_fixed_ball = pxy - px * py;
        bcov_fixed_ball *= bcov_fixed_ball;
        bcov_weight0 += bcov_fixed_ball;
        bcov_weight_prob += bcov_fixed_ball / (px * py);
        if (px != 1 && py != 1) {
            bcov_weight_hhg += bcov_fixed_ball / ((px) * (1.0 - px) * (py) * (1.0 - py));
            hhg_ball_num += 1;
        }
    }
    bcov_stat[0] = bcov_weight0 / (1.0 * (*n) * (*n));
    bcov_stat[1] = bcov_weight_prob / (1.0 * (*n) * (*n));
    bcov_stat[2] = hhg_ball_num > 0 ? (bcov_weight_hhg / hhg_ball_num) : 0.0;

    // free memory
    free(isource);
    free(icount);
    free(xy_index);
    free(yrank);
    free(xy_temp);
    free_int_matrix(xyidx, *n, *n);
}

void BI(double *bcov, double *pvalue, double *x, double *y, int *n, int *R, int *thread) {
    int i, j, **xidx, **yidx, *i_perm, *i_perm_inv, x_ties = 0, y_ties = 0, xy_ties = 0;
    int **x_within_ball = alloc_int_matrix(*n, *n);
    int **y_within_ball = alloc_int_matrix(*n, *n);
    double **Dx, **Dy, *x_cpy, *y_cpy;

    Dx = alloc_matrix(*n, *n);
    Dy = alloc_matrix(*n, *n);
    xidx = alloc_int_matrix(*n, *n);
    yidx = alloc_int_matrix(*n, *n);
    i_perm = (int *) malloc(*n * sizeof(int));
    i_perm_inv = (int *) malloc(*n * sizeof(int));
    x_cpy = (double *) malloc(*n * sizeof(double));
    y_cpy = (double *) malloc(*n * sizeof(double));

    // convert distance vector to distance matrix:
    distance2matrix(x, Dx, *n);
    distance2matrix(y, Dy, *n);

    for (i = 0; i < *n; i++) {
        for (j = 0; j < *n; j++) {
            xidx[i][j] = j;
            yidx[i][j] = j;
        }
        i_perm[i] = i;
        i_perm_inv[i] = i;
    }

    // sort each row of Dx and Dy, time complexity: O(n^2 * log{n})
    // xidx, yidx is index of each row after sorted
    for (i = 0; i < (*n); i++) {
        memcpy(x_cpy, Dx[i], *n * sizeof(double));
        memcpy(y_cpy, Dy[i], *n * sizeof(double));
        quicksort(x_cpy, xidx[i], 0, *n - 1);
        quicksort(y_cpy, yidx[i], 0, *n - 1);
        if (!x_ties) {
            for (j = 1; j < *n; ++j) {
                if (x_cpy[j] == x_cpy[j - 1]) {
                    x_ties = 1;
                }
            }
        }
        if (!y_ties) {
            for (j = 1; j < *n; ++j) {
                if (y_cpy[j] == y_cpy[j - 1]) {
                    y_ties = 1;
                }
            }
        }
    }
    xy_ties = x_ties && y_ties;
    free(x_cpy);
    free(y_cpy);

    if (xy_ties) {
        for (i = 0; i < *n; ++i) {
            quick_rank_max_with_index(Dx[i], xidx[i], x_within_ball[i], *n);
            quick_rank_max_with_index(Dy[i], yidx[i], y_within_ball[i], *n);
        }
        Ball_Information2(bcov, n, Dx, Dy, x_within_ball, y_within_ball, i_perm, i_perm_inv);
//        Ball_Information(bcov, n, Dx, Dy, xidx, yidx, i_perm, i_perm_inv);
    } else if (x_ties) {
        for (i = 0; i < *n; ++i) {
            quick_rank_max_with_index(Dx[i], xidx[i], y_within_ball[i], *n);
        }
        Ball_Information_NoTies(bcov, n, y_within_ball, yidx, Dx, i_perm);
    } else if (y_ties) {
        for (i = 0; i < *n; ++i) {
            quick_rank_max_with_index(Dy[i], yidx[i], y_within_ball[i], *n);
        }
        Ball_Information_NoTies(bcov, n, y_within_ball, xidx, Dy, i_perm);
    } else {
        for (i = 0; i < *n; ++i) {
            for (j = 0; j < *n; ++j) {
                y_within_ball[i][yidx[i][j]] = j + 1;
            }
        }
        Ball_Information_NoTies(bcov, n, y_within_ball, xidx, Dy, i_perm);
    }

    if (*R > 0) {
        double *permuted_bcov_weight0, *permuted_bcov_weight_prob, *permuted_bcov_weight_hhg;
        permuted_bcov_weight0 = (double *) malloc(*R * sizeof(double));
        permuted_bcov_weight_prob = (double *) malloc(*R * sizeof(double));
        permuted_bcov_weight_hhg = (double *) malloc(*R * sizeof(double));
        int not_parallel = *thread == 1 ? 1 : 0;
        if (xy_ties) {
            if (not_parallel) {
                double bcov_tmp[3];
                for (i = 0; i < *R; i++) {
                    // stop permutation if user stop it manually:
                    if (pending_interrupt()) {
                        print_stop_message();
                        break;
                    }
                    resample2(i_perm, n);
                    Ball_Information2(bcov_tmp, n, Dx, Dy, x_within_ball, y_within_ball, i_perm, i_perm_inv);
//                    resample(i_perm, i_perm_inv, n);
//                    Ball_Information(bcov_tmp, n, Dx, Dy, xidx, yidx, i_perm, i_perm_inv);
                    permuted_bcov_weight0[i] = bcov_tmp[0];
                    permuted_bcov_weight_prob[i] = bcov_tmp[1];
                    permuted_bcov_weight_hhg[i] = bcov_tmp[2];
                }
            } else {
                int **i_perm_matrix, **i_perm_inv_matrix;
                i_perm_matrix = alloc_int_matrix(*R, *n);
                i_perm_inv_matrix = alloc_int_matrix(*R, *n);
                shuffle_indicator_inv_matrix(i_perm_matrix, i_perm_inv_matrix, i_perm, i_perm_inv, *R, *n);
#pragma omp parallel
                {
                    int i_thread;
                    double bcov_tmp_thread[3];
#pragma omp for
                    for (i_thread = 0; i_thread < (*R); i_thread++) {
                        Ball_Information2(bcov_tmp_thread, n, Dx, Dy, x_within_ball, y_within_ball,
                                          i_perm_matrix[i_thread], i_perm_inv_matrix[i_thread]);
//                        Ball_Information(bcov_tmp_thread, n, Dx, Dy, xidx, yidx,
//                                         i_perm_matrix[i_thread], i_perm_inv_matrix[i_thread]);
                        permuted_bcov_weight0[i_thread] = bcov_tmp_thread[0];
                        permuted_bcov_weight_prob[i_thread] = bcov_tmp_thread[1];
                        permuted_bcov_weight_hhg[i_thread] = bcov_tmp_thread[2];
                    }
                }
                i = *R;
                free_int_matrix(i_perm_matrix, *R, *n);
                free_int_matrix(i_perm_inv_matrix, *R, *n);
            }
        } else {
            if (not_parallel) {
                double bcov_tmp[3];
                for (i = 0; i < *R; i++) {
                    // stop permutation if user stop it manually:
                    if (pending_interrupt()) {
                        print_stop_message();
                        break;
                    }
                    resample2(i_perm, n);
                    if (x_ties) {
                        Ball_Information_NoTies(bcov_tmp, n, y_within_ball, yidx, Dx, i_perm);
                    } else {
                        Ball_Information_NoTies(bcov_tmp, n, y_within_ball, xidx, Dy, i_perm);
                    }
                    permuted_bcov_weight0[i] = bcov_tmp[0];
                    permuted_bcov_weight_prob[i] = bcov_tmp[1];
                    permuted_bcov_weight_hhg[i] = bcov_tmp[2];
                }
            } else {
                int **i_perm_matrix = alloc_int_matrix(*R, *n);
                resample2_matrix(i_perm_matrix, i_perm, *R, *n);
#pragma omp parallel
                {
                    int i_thread;
                    double bcov_tmp_thread[3];
#pragma omp for
                    for (i_thread = 0; i_thread < (*R); i_thread++) {
                        if (x_ties) {
                            Ball_Information_NoTies(bcov_tmp_thread, n, y_within_ball, yidx, Dx,
                                                    i_perm_matrix[i_thread]);
                        } else {
                            Ball_Information_NoTies(bcov_tmp_thread, n, y_within_ball, xidx, Dy,
                                                    i_perm_matrix[i_thread]);
                        }
                        permuted_bcov_weight0[i_thread] = bcov_tmp_thread[0];
                        permuted_bcov_weight_prob[i_thread] = bcov_tmp_thread[1];
                        permuted_bcov_weight_hhg[i_thread] = bcov_tmp_thread[2];
                    }
                }
                i = *R;
                free_int_matrix(i_perm_matrix, *R, *n);
            }
        }

        pvalue[0] = compute_pvalue(bcov[0], permuted_bcov_weight0, i);
        pvalue[1] = compute_pvalue(bcov[1], permuted_bcov_weight_prob, i);
        pvalue[2] = compute_pvalue(bcov[2], permuted_bcov_weight_hhg, i);

        free(permuted_bcov_weight0);
        free(permuted_bcov_weight_prob);
        free(permuted_bcov_weight_hhg);
    }

    free_matrix(Dx, *n, *n);
    free_matrix(Dy, *n, *n);
    free_int_matrix(xidx, *n, *n);
    free_int_matrix(yidx, *n, *n);
    free(i_perm);
    free(i_perm_inv);
    free_int_matrix(x_within_ball, *n, *n);
    free_int_matrix(y_within_ball, *n, *n);
}

void U_Ball_Information(double *bcov_stat, int *n, int **Rank,
                        int **lowxidx, int **higxidx, int **lowyidx, int **higyidx,
                        int *i_perm) {
    int i, j, pi, pj;
    double px, py, pxy;
    double bcov_weight0 = 0.0, bcov_weight_prob = 0.0, bcov_weight_hhg = 0.0, bcov_fixed_ball = 0.0;
    double hhg_ball_num = 0.0;
    for (i = 0; i < *n; i++) {
        for (j = 0; j < *n; j++) {
            pi = i_perm[i];
            pj = i_perm[j];
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
        }
    }
    bcov_stat[0] = bcov_weight0 / (1.0 * (*n) * (*n));
    bcov_stat[1] = bcov_weight_prob / (1.0 * (*n) * (*n));
    bcov_stat[2] = hhg_ball_num > 0 ? (bcov_weight_hhg / hhg_ball_num) : 0.0;
}

void UBI(double *bcov, double *pvalue, double *x, double *y, int *n, int *R, int *thread) {
    int i, j, *xidx, *yidx, *xrank, *yrank, *i_perm, **Rank, **lowxidx, **higxidx, **lowyidx, **higyidx;

    xidx = (int *) malloc(*n * sizeof(int));
    yidx = (int *) malloc(*n * sizeof(int));
    xrank = (int *) malloc(*n * sizeof(int));
    yrank = (int *) malloc(*n * sizeof(int));
    i_perm = (int *) malloc(*n * sizeof(int));
    Rank = alloc_int_matrix((*n) + 1, (*n) + 1);
    lowxidx = alloc_int_matrix(*n, *n);
    higxidx = alloc_int_matrix(*n, *n);
    lowyidx = alloc_int_matrix(*n, *n);
    higyidx = alloc_int_matrix(*n, *n);

    for (i = 0; i < *n; i++) {
        xidx[i] = i;
        yidx[i] = i;
        i_perm[i] = i;
    }

    // first step: sort the x, y and obtain index of x, y (after the x, y have been sorted)
    quicksort(x, xidx, 0, *n - 1);
    quicksort(y, yidx, 0, *n - 1);
    ranksort(n, xrank, x, xidx);
    ranksort(n, yrank, y, yidx);
    createidx(n, xidx, x, lowxidx, higxidx);
    createidx(n, yidx, y, lowyidx, higyidx);

    initRank(*n, Rank, xrank, yrank, i_perm);
    U_Ball_Information(bcov, n, Rank, lowxidx, higxidx, lowyidx, higyidx, i_perm);
    if (*R > 0) {
        double bcov_tmp[3], *permuted_bcov_weight0, *permuted_bcov_weight_prob, *permuted_bcov_weight_hhg;
        permuted_bcov_weight0 = (double *) malloc(*R * sizeof(double));
        permuted_bcov_weight_prob = (double *) malloc(*R * sizeof(double));
        permuted_bcov_weight_hhg = (double *) malloc(*R * sizeof(double));

        int not_parallel = *thread == 1 ? 1 : 0;
        if (not_parallel) {
            for (j = 0; j < *R; j++) {
                // stop permutation if user stop it manually:
                if (pending_interrupt()) {
                    print_stop_message();
                    break;
                }
                resample2(i_perm, n);
                initRank(*n, Rank, xrank, yrank, i_perm);
                U_Ball_Information(bcov_tmp, n, Rank, lowxidx, higxidx, lowyidx, higyidx, i_perm);
                permuted_bcov_weight0[j] = bcov_tmp[0];
                permuted_bcov_weight_prob[j] = bcov_tmp[1];
                permuted_bcov_weight_hhg[j] = bcov_tmp[2];
            }
        } else {
            int **i_perm_thread = alloc_int_matrix(*R, *n);
            resample2_matrix(i_perm_thread, i_perm, *R, *n);
#pragma omp parallel
            {
                int j_thread, **Rank_thread;
                double bcov_tmp_thread[3];
                Rank_thread = alloc_int_matrix((*n) + 1, (*n) + 1);
#pragma omp for
                for (j_thread = 0; j_thread < (*R); j_thread++) {
                    initRank(*n, Rank_thread, xrank, yrank, i_perm_thread[j_thread]);
                    U_Ball_Information(bcov_tmp_thread, n, Rank_thread, lowxidx, higxidx, lowyidx, higyidx,
                                       i_perm_thread[j_thread]);
                    permuted_bcov_weight0[j_thread] = bcov_tmp_thread[0];
                    permuted_bcov_weight_prob[j_thread] = bcov_tmp_thread[1];
                    permuted_bcov_weight_hhg[j_thread] = bcov_tmp_thread[2];
                }
                free_int_matrix(Rank_thread, (*n) + 1, (*n) + 1);
            }
            free_int_matrix(i_perm_thread, *R, *n);
            j = *R;
        }

        pvalue[0] = compute_pvalue(bcov[0], permuted_bcov_weight0, j);
        pvalue[1] = compute_pvalue(bcov[1], permuted_bcov_weight_prob, j);
        pvalue[2] = compute_pvalue(bcov[2], permuted_bcov_weight_hhg, j);

        free(permuted_bcov_weight0);
        free(permuted_bcov_weight_prob);
        free(permuted_bcov_weight_hhg);
    }

    free_int_matrix(Rank, (*n) + 1, (*n) + 1);
    free_int_matrix(lowxidx, *n, *n);
    free_int_matrix(higxidx, *n, *n);
    free_int_matrix(lowyidx, *n, *n);
    free_int_matrix(higyidx, *n, *n);
    free(xidx);
    free(yidx);
    free(xrank);
    free(yrank);
    free(i_perm);
}

void bcov_test(double *bcov, double *pvalue, double *x, double *y, int *n, int *R, int *dst, int *thread) {
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
    if ((*dst)) {
        BI(bcov, pvalue, x, y, n, R, thread);
    } else {
        UBI(bcov, pvalue, x, y, n, R, thread);
    }
}
