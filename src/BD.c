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
#include "stdio.h"
#include "BD.h"
#include "utilities.h"
#include "Ball_omp.h"

void Ball_Divergence(double *bd_stat, int **Rxy, int **Rx, int *i_perm_tmp, int *n1, int *n2) {
    int i, j, n;
    double TS_weight0 = 0, SS_weight0 = 0, TS_weight1 = 0.0, SS_weight1 = 0.0;
    double p1, p2, p3, ans;
    n = *n1 + *n2;

    // Calculate A_{ij}^{X} and A_{ij}^{Y}:
    for (i = 0; i < *n1; i++) {
        for (j = 0; j < *n1; j++) {
            p1 = Rx[i_perm_tmp[i]][i_perm_tmp[j]] + 1;
            p2 = Rxy[i_perm_tmp[i]][i_perm_tmp[j]] - p1 + 1;
            p3 = (p1 + p2) / n;
            if (p3 * (1 - p3) == 0) { continue; }
            //TS++;
            //TS += pow(fabs(p1/(*n1)-p2/(*n2)),*gamma)/(pow(p3,*alpha)*pow(1-p3,*beta));
            ans = p1 / (*n1) - p2 / (*n2);
            TS_weight0 += (ans * ans);
            // TODO: Weight Ball Divergence:
            TS_weight1 += (ans * ans) * 1.0;
        }
    }
    // Calculate C_{kl}^{X} and C_{kl}^{Y}:
    for (i = *n1; i < n; i++) {
        for (j = *n1; j < n; j++) {
            p1 = Rx[i_perm_tmp[i]][i_perm_tmp[j]] + 1;
            p2 = Rxy[i_perm_tmp[i]][i_perm_tmp[j]] - p1 + 1;
            p3 = (p1 + p2) / n;
            if (p3 * (1 - p3) == 0) { continue; }

            //SS++;
            //SS += pow(fabs(p1/(*n1)-p2/(*n2)),*gamma)/(pow(p3,*alpha)*pow(1-p3,*beta));
            ans = p1 / (*n1) - p2 / (*n2);
            SS_weight0 += (ans * ans);
            // TODO: Weight Ball Divergence:
            SS_weight1 += (ans * ans) * 1.0;
        }
    }
    bd_stat[0] = TS_weight0 / (1.0 * (*n1) * (*n1)) + SS_weight0 / (1.0 * (*n2) * (*n2));
    bd_stat[1] = TS_weight1 / (1.0 * (*n1) * (*n1)) + SS_weight1 / (1.0 * (*n2) * (*n2));
}

/**
 * @deprecated
 */
void Ball_Divergence_parallel(double *bd_stat, int **Rxy, int **Rx, int *i_perm_tmp, int *n1, int *n2, int *nthread) {
    int n;
    n = *n1 + *n2;
    double TS_weight0 = 0, SS_weight0 = 0, TS_weight1 = 0.0, SS_weight1 = 0.0;
#ifdef Ball_OMP_H_
    omp_set_num_threads(*nthread);
#endif
#pragma omp parallel
    {
        int i, j;
        double ans;
        double TS_value_w0 = 0.0, SS_value_w0 = 0.0, TS_value_w1 = 0.0, SS_value_w1 = 0.0;
        double p1, p2, p3;

        // Calculate A_{ij}^{X} and A_{ij}^{Y}
        // TODO: further optimized by : pragma omp for reduction ...
#pragma omp for
        for (i = 0; i < *n1; i++) {
            for (j = 0; j < *n1; j++) {
                p1 = Rx[i_perm_tmp[i]][i_perm_tmp[j]] + 1;
                p2 = Rxy[i_perm_tmp[i]][i_perm_tmp[j]] - p1 + 1;
                p3 = (p1 + p2) / n;
                if (p3 * (1 - p3) == 0) { continue; }
                ans = p1 / (*n1) - p2 / (*n2);
                TS_value_w0 += (ans * ans);
                TS_value_w1 += (ans * ans) * 1.0;
            }
        }
#pragma omp critical
        {
            TS_weight0 += TS_value_w0;
            TS_weight1 += TS_value_w1;
        }

        // Calculate C_{ij}^{X} and C_{ij}^{Y}:
#pragma omp for
        for (i = *n1; i < n; i++) {
            for (j = *n1; j < n; j++) {
                p1 = Rx[i_perm_tmp[i]][i_perm_tmp[j]] + 1;
                p2 = Rxy[i_perm_tmp[i]][i_perm_tmp[j]] - p1 + 1;
                p3 = (p1 + p2) / n;
                if (p3 * (1 - p3) == 0) { continue; }
                ans = p1 / (*n1) - p2 / (*n2);
                SS_value_w0 += (ans * ans);
                SS_value_w1 += (ans * ans) * 1.0;
            }
        }
#pragma omp critical
        {
            SS_weight0 += SS_value_w0;
            SS_weight1 += SS_value_w1;
        }

    }
    bd_stat[0] = TS_weight0 / (1.0 * (*n1) * (*n1)) + SS_weight0 / (1.0 * (*n2) * (*n2));
    bd_stat[1] = TS_weight1 / (1.0 * (*n1) * (*n1)) + SS_weight1 / (1.0 * (*n2) * (*n2));
}

/**
 * @deprecated
 */
void Ball_Divergence_wrapper(double *bd_stat, int **Rxy, int **Rx, int *i_perm_tmp, int *n1, int *n2, int *nthread) {
    if (*nthread > 1) {
        Ball_Divergence_parallel(bd_stat, Rxy, Rx, i_perm_tmp, n1, n2, nthread);
    } else {
        Ball_Divergence(bd_stat, Rxy, Rx, i_perm_tmp, n1, n2);
    }
}


/*
 * Two sample ball divergence test, based on permutation technique (for multivariate data)
 * input:
 * xy: vectorized distance matrix calculate used original data
 */
void BD(double *bd, double *pvalue, double *xy, int *n1, int *n2, int *R, int *thread) {
    //  computes TST(x,y)
    int i, j, n;
    /*
    Dxy: distance matrix
    Ixy: each row corresponding to sample index
    Rxy: each row corresponding the rank in each row (all data)
    Rx: each row corresponding the rank in each row (group 1)
    i_perm: group indicator (value 1 corresponding to group 1, and value 0 corresponding to group 0)
    i_perm_tmp: index of sample, value: 0, 1, 2, ..., n - 1;
    */
    int **Ixy, **Rxy, **Rx, *i_perm, *i_perm_tmp;
    double **Dxy;

    n = *n1 + *n2;
    Dxy = alloc_matrix(n, n);
    Ixy = alloc_int_matrix(n, n);
    Rx = alloc_int_matrix(n, n);
    Rxy = alloc_int_matrix(n, n);
    i_perm = (int *) malloc(n * sizeof(int));
    i_perm_tmp = (int *) malloc(n * sizeof(int));

    // get vectorized distance matrix Dxy:
    distance2matrix(xy, Dxy, n);

    // compute ball divergence:
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            Ixy[i][j] = j;
        }
    }
    for (i = 0; i < n; i++) {
        i_perm[i] = i < (*n1) ? 1 : 0;
        i_perm_tmp[i] = i;
    }
    for (i = 0; i < n; i++) {
        quicksort(Dxy[i], Ixy[i], 0, n - 1);
    }
    ranksort2(n, Rxy, Dxy, Ixy);
    Findx(Rxy, Ixy, i_perm, n1, n2, Rx);
    Ball_Divergence(bd, Rxy, Rx, i_perm_tmp, n1, n2);
    free_matrix(Dxy, n, n);

    if (*R > 0) {
        double *permuted_bd_w0, *permuted_bd_w1;
        permuted_bd_w0 = (double *) malloc(*R * sizeof(double));
        permuted_bd_w1 = (double *) malloc(*R * sizeof(double));
        int not_parallel = *thread == 1 ? 1 : 0;
#ifdef Ball_OMP_H_
        if (not_parallel != 1) {
            omp_set_dynamic(0);
            if (*thread < 0) {
                omp_set_num_threads(omp_get_num_procs());
            } else {
                omp_set_num_threads(*thread);
            }
        }
#endif
        if (not_parallel) {
            double bd_tmp[2];
#ifdef R_BUILD
            GetRNGstate();
#endif
            for (i = 0; i < *R; i++) {
                // stop permutation if user stop it manually:
                if (pending_interrupt()) {
                    print_stop_message();
                    break;
                }
                resample_indicator_label(i_perm, i_perm_tmp, n, n1);
                Findx(Rxy, Ixy, i_perm, n1, n2, Rx);
                Ball_Divergence(bd_tmp, Rxy, Rx, i_perm_tmp, n1, n2);
                permuted_bd_w0[i] = bd_tmp[0];
                permuted_bd_w1[i] = bd_tmp[1];
            }
#ifdef R_BUILD
            PutRNGstate();
#endif
            free(i_perm);
            free(i_perm_tmp);
        } else {
            int **i_perm_matrix, **i_perm_tmp_matrix;
            i_perm_matrix = alloc_int_matrix(*R, n);
            i_perm_tmp_matrix = alloc_int_matrix(*R, n);
            resample_indicator_label_matrix(i_perm_matrix, i_perm_tmp_matrix, i_perm, i_perm_tmp, *R, n, n1);
#pragma omp parallel
            {
                int **Rx_thread = alloc_int_matrix(n, n);
                int i_thread;
                double bd_tmp[2];
#pragma omp for
                for (i_thread = 0; i_thread < (*R); i_thread++) {
                    Findx(Rxy, Ixy, i_perm_matrix[i_thread], n1, n2, Rx_thread);
                    Ball_Divergence(bd_tmp, Rxy, Rx_thread, i_perm_tmp_matrix[i_thread], n1, n2);
                    permuted_bd_w0[i_thread] = bd_tmp[0];
                    permuted_bd_w1[i_thread] = bd_tmp[1];
                }
                free_int_matrix(Rx_thread, n, n);
            }
            free_int_matrix(i_perm_matrix, *R, n);
            free_int_matrix(i_perm_tmp_matrix, *R, n);
            i = *R;
        }
        pvalue[0] = compute_pvalue(bd[0], permuted_bd_w0, i);
        pvalue[1] = compute_pvalue(bd[1], permuted_bd_w1, i);

        free(permuted_bd_w0);
        free(permuted_bd_w1);
        free_int_matrix(Ixy, n, n);
        free_int_matrix(Rxy, n, n);
        free_int_matrix(Rx, n, n);
    }
}

/**
 * Note:
 * I have tried to use following to avoid the lock during the parallelism. However, this method is not efficient.
 * int **random_permute_matrix;
 * random_permute_matrix = alloc_int_matrix(n, n);
 * for (i = 0; i < *R; i++) {
 *      for (j = n - 1; j > 0; j--) {
 *          random_permute_matrix[i][j] = random_index_thread_wrap(j);
 *      }
 *      random_permute_matrix[i][0] = 0;
 * }
 * ....
 * resample3_thread(random_permute_matrix[i_thread], i_perm_thread, i_perm_tmp_thread, n, n1);
 * ....
 * free_int_matrix(random_permute_matrix, n, n);
 * @deprecated
 */
void BD_parallel(double *bd, double *pvalue, double *xy, int *n1, int *n2, int *p, int *R, int *nthread) {
    int i, j, n;
    int **Ixy, **Rxy, **Rx, *i_perm, *i_perm_tmp;
    double **Dxy;

    n = *n1 + *n2;
    Dxy = alloc_matrix(n, n);
    Ixy = alloc_int_matrix(n, n);
    Rx = alloc_int_matrix(n, n);
    Rxy = alloc_int_matrix(n, n);
    i_perm = (int *) malloc(n * sizeof(int));
    i_perm_tmp = (int *) malloc(n * sizeof(int));

    // get vectorized distance matrix Dxy:
    vector2matrix(xy, Dxy, n, n, 1);

    for (i = 0; i < n; i++)
        for (j = 0; j < n; j++)
            Ixy[i][j] = j;

    for (i = 0; i < n; i++) {
        i_perm[i] = i < (*n1) ? 1 : 0;
        i_perm_tmp[i] = i;
    }
    for (i = 0; i < n; i++) {
        quicksort(Dxy[i], Ixy[i], 0, n - 1);
    }
    ranksort2(n, Rxy, Dxy, Ixy);
    free_matrix(Dxy, n, n);
    Findx(Rxy, Ixy, i_perm, n1, n2, Rx);
    Ball_Divergence(bd, Rxy, Rx, i_perm_tmp, n1, n2);

    free_int_matrix(Rx, n, n);
    free(i_perm);
    free(i_perm_tmp);

    if (*R > 0) {
        double *permuted_bd_w0, *permuted_bd_w1;
        permuted_bd_w0 = (double *) malloc(*R * sizeof(double));
        permuted_bd_w1 = (double *) malloc(*R * sizeof(double));

#ifdef Ball_OMP_H_
        omp_set_num_threads(*nthread);
#endif
        // Init parallel
#pragma omp parallel
        {
            int **Rx_thread, *i_perm_thread, *i_perm_tmp_thread;
            int k, i_thread;
            double *bd_tmp;
            Rx_thread = alloc_int_matrix(n, n);
            i_perm_thread = (int *) malloc(n * sizeof(int));
            i_perm_tmp_thread = (int *) malloc(n * sizeof(int));

            for (k = 0; k < n; k++) {
                i_perm_thread[k] = k < (*n1) ? 1 : 0;
                i_perm_tmp_thread[k] = k;
            }

#pragma omp for
            for (i_thread = 0; i_thread < (*R); i_thread++) {
                bd_tmp = (double *) malloc(2 * sizeof(double));
#pragma omp critical
                {
                    resample_indicator_label(i_perm_thread, i_perm_tmp_thread, n, n1);
                }
                Findx(Rxy, Ixy, i_perm_thread, n1, n2, Rx_thread);
                Ball_Divergence(bd_tmp, Rxy, Rx_thread, i_perm_tmp_thread, n1, n2);
                permuted_bd_w0[i_thread] = bd_tmp[0];
                permuted_bd_w1[i_thread] = bd_tmp[1];
                free(bd_tmp);
            }
            free(i_perm_thread);
            free(i_perm_tmp_thread);
            free_int_matrix(Rx_thread, n, n);
        }

        pvalue[0] = compute_pvalue(bd[0], permuted_bd_w0, *R);
        pvalue[1] = compute_pvalue(bd[1], permuted_bd_w1, *R);

        free(permuted_bd_w0);
        free(permuted_bd_w1);
    }
    // free matrix:
    free_int_matrix(Ixy, n, n);
    free_int_matrix(Rxy, n, n);
}


/**
 * Two sample ball divergence test, based on permutation technique (for univariate data)
 * input:
 * xy: original data
 */
void UBD(double *bd, double *pvalue, double *xy, int *n1, int *n2, int *R, int *thread) {
    //  computes TST(x,y)
    int i, j, n;
    int **Ixy, **Rxy, **Rx, *i_perm, *i_perm_tmp, *xyidx;

    n = *n1 + *n2;
    Ixy = alloc_int_matrix(n, n);
    Rx = alloc_int_matrix(n, n);
    Rxy = alloc_int_matrix(n, n);
    i_perm = (int *) malloc(n * sizeof(int));
    i_perm_tmp = (int *) malloc(n * sizeof(int));
    xyidx = (int *) malloc(n * sizeof(int));
    for (i = 0; i < n; i++) {
        xyidx[i] = i;
        for (j = 0; j < n; j++) { Ixy[i][j] = j; }
    }

    for (i = 0; i < n; i++) {
        i_perm[i] = i < *n1 ? 1 : 0;
        i_perm_tmp[i] = i;
    }

    quicksort(xy, xyidx, 0, n - 1);
    ranksort3(n, xyidx, xy, Rxy, Ixy);
    free(xyidx);
    Findx(Rxy, Ixy, i_perm, n1, n2, Rx);
    Ball_Divergence(bd, Rxy, Rx, i_perm_tmp, n1, n2);

    if (*R > 0) {
        double *permuted_bd_w0, *permuted_bd_w1;
        permuted_bd_w0 = (double *) malloc(*R * sizeof(double));
        permuted_bd_w1 = (double *) malloc(*R * sizeof(double));

        int not_parallel = *thread == 1 ? 1 : 0;
#ifdef Ball_OMP_H_
        if (not_parallel != 1) {
            omp_set_dynamic(0);
            if (*thread < 0) {
                omp_set_num_threads(omp_get_num_procs());
            } else {
                omp_set_num_threads(*thread);
            }
        }
#endif

        if (not_parallel) {
            double bd_tmp[2];
            for (i = 0; i < *R; i++) {
                if (pending_interrupt()) {
                    print_stop_message();
                    break;
                }
                resample_indicator_label(i_perm, i_perm_tmp, n, n1);
                Findx(Rxy, Ixy, i_perm, n1, n2, Rx);
                Ball_Divergence(bd_tmp, Rxy, Rx, i_perm_tmp, n1, n2);
                permuted_bd_w0[i] = bd_tmp[0];
                permuted_bd_w1[i] = bd_tmp[1];
            }
        } else {
            int **i_perm_matrix, **i_perm_tmp_matrix;
            i_perm_matrix = alloc_int_matrix(*R, n);
            i_perm_tmp_matrix = alloc_int_matrix(*R, n);
            resample_indicator_label_matrix(i_perm_matrix, i_perm_tmp_matrix, i_perm, i_perm_tmp, *R, n, n1);

#pragma omp parallel
            {
                int i_thread, **Rx_thread;
                double bd_tmp_thread[2];
                Rx_thread = alloc_int_matrix(n, n);
#pragma omp for
                for (i_thread = 0; i_thread < (*R); i_thread++) {
                    Findx(Rxy, Ixy, i_perm_matrix[i_thread], n1, n2, Rx_thread);
                    Ball_Divergence(bd_tmp_thread, Rxy, Rx_thread, i_perm_tmp_matrix[i_thread], n1, n2);
                    permuted_bd_w0[i_thread] = bd_tmp_thread[0];
                    permuted_bd_w1[i_thread] = bd_tmp_thread[1];
                }
                free_int_matrix(Rx_thread, n, n);
            }
            free_int_matrix(i_perm_matrix, *R, n);
            free_int_matrix(i_perm_tmp_matrix, *R, n);
            i = *R;
        }

        pvalue[0] = compute_pvalue(bd[0], permuted_bd_w0, i);
        pvalue[1] = compute_pvalue(bd[1], permuted_bd_w1, i);
        free(permuted_bd_w0);
        free(permuted_bd_w1);
    }

    free_int_matrix(Ixy, n, n);
    free_int_matrix(Rxy, n, n);
    free_int_matrix(Rx, n, n);
    free(i_perm);
    free(i_perm_tmp);
}

/**
 * @deprecated
 */
void UBD_parallel(double *bd, double *pvalue, double *xy, int *n1, int *n2, int *R, int *nthread) {
    //  computes TST(x,y)
    int i, j, n;
    int **Ixy, **Rxy, **Rx, *i_perm, *i_perm_tmp, *xyidx;

    n = *n1 + *n2;
    Ixy = alloc_int_matrix(n, n);
    Rx = alloc_int_matrix(n, n);
    Rxy = alloc_int_matrix(n, n);
    i_perm = (int *) malloc(n * sizeof(int));
    i_perm_tmp = (int *) malloc(n * sizeof(int));
    xyidx = (int *) malloc(n * sizeof(int));

    for (i = 0; i < n; i++) {
        xyidx[i] = i;
        for (j = 0; j < n; j++)
            Ixy[i][j] = j;
    }

    for (i = 0; i < n; i++) {
        i_perm[i] = i < *n1 ? 1 : 0;
        i_perm_tmp[i] = i;
    }

    quicksort(xy, xyidx, 0, n - 1);
    ranksort3(n, xyidx, xy, Rxy, Ixy);
    Findx(Rxy, Ixy, i_perm, n1, n2, Rx);
    Ball_Divergence(bd, Rxy, Rx, i_perm_tmp, n1, n2);

    free_int_matrix(Rx, n, n);
    free(i_perm);
    free(i_perm_tmp);
    free(xyidx);

    if (*R > 0) {
        double *permuted_bd_w0, *permuted_bd_w1;
        permuted_bd_w0 = (double *) malloc(*R * sizeof(double));
        permuted_bd_w1 = (double *) malloc(*R * sizeof(double));

#ifdef Ball_OMP_H_
        omp_set_num_threads(*nthread);
#endif
#pragma omp parallel
        {
            int **Rx_thread, *i_perm_thread, *i_perm_tmp_thread;
            int k, i_thread;
            double *bd_tmp;
            Rx_thread = alloc_int_matrix(n, n);
            i_perm_thread = (int *) malloc(n * sizeof(int));
            i_perm_tmp_thread = (int *) malloc(n * sizeof(int));

#pragma omp critical
            {
                for (k = 0; k < n; k++) {
                    i_perm_thread[k] = k < (*n1) ? 1 : 0;
                    i_perm_tmp_thread[k] = k;
                }
            }
#pragma omp for
            for (i_thread = 0; i_thread < (*R); i_thread++) {
                bd_tmp = (double *) malloc(2 * sizeof(double));
#pragma omp critical
                {
                    resample_indicator_label(i_perm_thread, i_perm_tmp_thread, n, n1);
                }
                Findx(Rxy, Ixy, i_perm_thread, n1, n2, Rx_thread);
                Ball_Divergence(bd_tmp, Rxy, Rx_thread, i_perm_tmp_thread, n1, n2);
                permuted_bd_w0[i_thread] = bd_tmp[0];
                permuted_bd_w1[i_thread] = bd_tmp[1];
                free(bd_tmp);
            }
            free(i_perm_thread);
            free(i_perm_tmp_thread);
            free_int_matrix(Rx_thread, n, n);
        }
        pvalue[0] = compute_pvalue(bd[0], permuted_bd_w0, *R);
        pvalue[1] = compute_pvalue(bd[1], permuted_bd_w1, *R);
        free(permuted_bd_w0);
        free(permuted_bd_w1);
    }
    //
    free_int_matrix(Ixy, n, n);
    free_int_matrix(Rxy, n, n);
}


/*
 * return ball divergence value of two sample with sample size n1, n2;
 * permutation procedure not consider at all
 * input:
 * xy: vectorized distance matrix
 */
void bd_value(double *bd_stat, double *xy, int *n1, int *n2) {
    int i, j, n;
    int **Ixy, **Rxy, **Rx, *i_perm;
    double **Dxy;
    double p1, p2, p3, ans, TS_weight0 = 0, SS_weight0 = 0, TS_weight1 = 0.0, SS_weight1 = 0.0;

    n = *n1 + *n2;
    Dxy = alloc_matrix(n, n);
    Ixy = alloc_int_matrix(n, n);
    Rx = alloc_int_matrix(n, n);
    Rxy = alloc_int_matrix(n, n);
    i_perm = (int *) malloc(n * sizeof(int));
    vector2matrix(xy, Dxy, n, n, 1);
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            Ixy[i][j] = j;
        }
    }

    for (i = 0; i < n; i++) {
        i_perm[i] = i < (*n1) ? 1 : 0;
    }
    for (i = 0; i < n; i++) {
        quicksort(Dxy[i], Ixy[i], 0, n - 1);
    }
    ranksort2(n, Rxy, Dxy, Ixy);
    free_matrix(Dxy, n, n);
    Findx(Rxy, Ixy, i_perm, n1, n2, Rx);
    free_int_matrix(Ixy, n, n);
    free(i_perm);

    // Calculate A_{ij}^{X} and A_{ij}^{Y}:
    for (i = 0; i < *n1; i++) {
        for (j = 0; j < *n1; j++) {
            p1 = Rx[i][j] + 1;
            p2 = Rxy[i][j] - p1 + 1;
            p3 = (p1 + p2) / n;
            if (p3 * (1 - p3) == 0) { continue; }
            ans = p1 / (*n1) - p2 / (*n2);
            TS_weight0 += (ans * ans);
            // TODO: Weight Ball Divergence:
            TS_weight1 += (ans * ans) * 1.0;
        }
    }
    // Calculate C_{kl}^{X} and C_{kl}^{Y}:
    for (i = *n1; i < n; i++) {
        for (j = *n1; j < n; j++) {
            p1 = Rx[i][j] + 1;
            p2 = Rxy[i][j] - p1 + 1;
            p3 = (p1 + p2) / n;
            if (p3 * (1 - p3) == 0) { continue; }
            ans = p1 / (*n1) - p2 / (*n2);
            SS_weight0 += (ans * ans);
            // TODO: Weight Ball Divergence:
            SS_weight1 += (ans * ans) * 1.0;
        }
    }
    bd_stat[0] = TS_weight0 / (1.0 * (*n1) * (*n1)) + SS_weight0 / (1.0 * (*n2) * (*n2));
    bd_stat[1] = TS_weight1 / (1.0 * (*n1) * (*n1)) + SS_weight1 / (1.0 * (*n2) * (*n2));

    // free:
    free_int_matrix(Rxy, n, n);
    free_int_matrix(Rx, n, n);
}

// input:
// size: an array contain sample size in each group
// cumulate_size: the cumulative sums of size vector  
// k: group number
void compute_cumulate_size(int *cumulate_size, int *size, int *k) {
    int i;
    for (i = 0; i < (*k); i++) {
        if (i == 0) {
            cumulate_size[i] = 0;
        } else {
            cumulate_size[i] = cumulate_size[i - 1] + size[i - 1];
        }
    }
}


// need to fetch sub-matrix and vectorized it
// xy: vectorized distance matrix
// ij_dst: vectorized distance matrix of group i and group j
// cumulate_size: the cumulative sums of size vector  
// size: an array contain sample size in each group
// *p < *q is needed and *p == i, *q == j
void get_ij_dst(double *xy, double *ij_dst, int *cumulate_size, int *size, int *n, int *p, int *q) {
    int k = 0, k1 = 0, k2 = 0;
    int i = *p;
    int j = *q;
    int n1 = size[i], n2 = size[j];
    int num = n1 + n2;
    // for group1:
    int index_ii = (*n) * cumulate_size[i] + cumulate_size[i];
    int index_ij = (*n) * cumulate_size[i] + cumulate_size[j];
    for (k1 = 0; k1 < n1; k1++) {
        for (k2 = 0; k2 < num; k2++) {
            if (k2 < n1) {
                ij_dst[k] = xy[index_ii + k2];
            } else {
                ij_dst[k] = xy[index_ij + k2 - n1];
            }
            k = k + 1;
        }
        index_ii = index_ii + (*n);
        index_ij = index_ij + (*n);
    }
    // for group2:
    index_ii = (*n) * cumulate_size[j] + cumulate_size[i];
    index_ij = (*n) * cumulate_size[j] + cumulate_size[j];
    for (k1 = 0; k1 < n2; k1++) {
        for (k2 = 0; k2 < num; k2++) {
            if (k2 < n1) {
                ij_dst[k] = xy[index_ii + k2];
            } else {
                ij_dst[k] = xy[index_ij + k2 - n1];
            }
            k = k + 1;
        }
        index_ii = index_ii + (*n);
        index_ij = index_ij + (*n);
    }
}


// return a distance matrix after permute index is given
// input: 
// xy: vectorized distance matrix
// new_xy: vectorized distance matrix store after sample was permuted
// index: permutation index
// N: sample size
void permute_dst(double *xy, double *new_xy, int *index, int *N) {
    int n = (*N);
    int row_index = 0, col_index = 0, exchange_index1 = 0, exchange_index2 = 0;
    for (int j = 0; j < (*N); j++) {
        for (int k = 0; k < (*N); k++) {
            row_index = index[j];
            col_index = index[k];
            exchange_index2 = row_index * n + col_index;
            new_xy[exchange_index1] = xy[exchange_index2];
            exchange_index1 = exchange_index1 + 1;
        }
    }
}


/**
 * @deprecated
 * Compute K-sample ball divergence
 * @param kbd_stat : K-samples ball divergence statistic
 * @param xy : vectorized distance matrix
 * @param size : an array contain sample size in each group
 * @param n : sample size
 * @param k : group number
 */
void kbd_value(double *kbd_stat, double *xy, int *size, int *n, int *k) {
    int K = *k;
    int i, j, s = 0, t = 0;
    int *cumulate_size;
    double *ij_dst, *bd_stat_w0_array, *bd_stat_w1_array, *bd_stat_w0_part_sum_array, *bd_stat_w1_part_sum_array;
    int two_group_size, tmp_value1 = 0, tmp_value2 = 0;
    int *n1 = &tmp_value1;
    int *n2 = &tmp_value2;
    double bd_stat_value[2];
    double kbd_stat_value_sum_w0 = 0.0, kbd_stat_value_sum_w1 = 0.0;
    double kbd_stat_value_max_w0 = 0.0, kbd_stat_value_max_w1 = 0.0;
    double kbd_stat_value_max1_w0, kbd_stat_value_max1_w1;
    int bd_stat_number = K * (K - 1) / 2;

    bd_stat_w0_array = (double *) malloc(bd_stat_number * sizeof(double));
    bd_stat_w1_array = (double *) malloc(bd_stat_number * sizeof(double));

    bd_stat_w0_part_sum_array = (double *) malloc(K * sizeof(double));
    bd_stat_w1_part_sum_array = (double *) malloc(K * sizeof(double));

    cumulate_size = (int *) malloc(K * sizeof(int));
    compute_cumulate_size(cumulate_size, size, k);

    // compute two type of KBD:
    for (i = 0; i < K; i++) {
        bd_stat_w0_part_sum_array[i] = 0.0;
        bd_stat_w1_part_sum_array[i] = 0.0;
        for (t = 0; t < (i - 1); t++) {
            bd_stat_w0_part_sum_array[i] += bd_stat_w0_array[(t + 1) * i - t * (t + 1) / 2];
            bd_stat_w1_part_sum_array[i] += bd_stat_w1_array[(t + 1) * i - t * (t + 1) / 2];
        }
        for (j = (i + 1); j < K; j++) {
            tmp_value1 = size[i];
            tmp_value2 = size[j];
            two_group_size = tmp_value1 + tmp_value2;
            two_group_size = two_group_size * two_group_size;
            ij_dst = (double *) malloc(two_group_size * sizeof(double));
            get_ij_dst(xy, ij_dst, cumulate_size, size, n, &i, &j);
            bd_value(bd_stat_value, ij_dst, n1, n2);
            // summation version:
            kbd_stat_value_sum_w0 += bd_stat_value[0];
            kbd_stat_value_sum_w1 += bd_stat_value[1];
            // maximum K-1 version:
            bd_stat_w0_array[s] = bd_stat_value[0];
            bd_stat_w1_array[s] = bd_stat_value[1];
            s += 1;
            // maxsum version:
            bd_stat_w0_part_sum_array[i] += bd_stat_value[0];
            bd_stat_w1_part_sum_array[i] += bd_stat_value[1];
            free(ij_dst);
        }
    }
    // compute maximum K-1 version statistic:
    quick_sort(bd_stat_w0_array, bd_stat_number);
    quick_sort(bd_stat_w1_array, bd_stat_number);
    for (i = bd_stat_number - 1; i > (bd_stat_number - K); i--) {
        kbd_stat_value_max_w0 += bd_stat_w0_array[i];
        kbd_stat_value_max_w1 += bd_stat_w1_array[i];
    }
    free(bd_stat_w0_array);
    free(bd_stat_w1_array);
    free(cumulate_size);
    // compute maximum version statistic:
    quick_sort(bd_stat_w0_part_sum_array, K);
    quick_sort(bd_stat_w1_part_sum_array, K);
    kbd_stat_value_max1_w0 = bd_stat_w0_part_sum_array[K - 1];
    kbd_stat_value_max1_w1 = bd_stat_w1_part_sum_array[K - 1];
    free(bd_stat_w0_part_sum_array);
    free(bd_stat_w1_part_sum_array);
    //
    kbd_stat[0] = kbd_stat_value_sum_w0;
    kbd_stat[1] = kbd_stat_value_sum_w1;
    kbd_stat[2] = kbd_stat_value_max_w0;
    kbd_stat[3] = kbd_stat_value_max_w1;
    kbd_stat[4] = kbd_stat_value_max1_w0;
    kbd_stat[5] = kbd_stat_value_max1_w1;
}

/**
 * Compute K-sample ball divergence
 * @param kbd_stat : K-samples ball divergence statistic
 * @param xy : vectorized distance matrix
 * @param cumsum_size : the cumulative sums of size vector
 * @param size : an array contain sample size in each group
 * @param n : sample size
 * @param k : group number
 */
void K_Ball_Divergence(double *kbd_stat, double *xy, int *cumsum_size, int *size, int *n, int *k) {
    int i, j, s = 0;
    double *ij_dst, *bd_stat_w0_array, *bd_stat_w1_array, *bd_stat_w0_sum_array, *bd_stat_w1_sum_array;
    int two_group_size;
    double bd_stat_value[2];
    double kbd_stat_value_sum_w0 = 0.0, kbd_stat_value_sum_w1 = 0.0;
    double kbd_stat_value_max_w0 = 0.0, kbd_stat_value_max_w1 = 0.0;
    double kbd_stat_value_maxsum_w0, kbd_stat_value_maxsum_w1;
    int bd_stat_number = *k * (*k - 1) / 2;

    bd_stat_w0_array = (double *) malloc(bd_stat_number * sizeof(double));
    bd_stat_w1_array = (double *) malloc(bd_stat_number * sizeof(double));
    bd_stat_w0_sum_array = (double *) malloc(*k * sizeof(double));
    bd_stat_w1_sum_array = (double *) malloc(*k * sizeof(double));
    for (int l = 0; l < *k; ++l) {
        bd_stat_w0_sum_array[l] = 0;
        bd_stat_w1_sum_array[l] = 0;
    }

    // compute two type of KBD:
    for (i = 0; i < *k; i++) {
        for (j = (i + 1); j < *k; j++) {
            two_group_size = size[i] + size[j];
            two_group_size = two_group_size * two_group_size;
            ij_dst = (double *) malloc(two_group_size * sizeof(double));
            get_ij_dst(xy, ij_dst, cumsum_size, size, n, &i, &j);
            bd_value(bd_stat_value, ij_dst, &size[i], &size[j]);
            // summation version:
            kbd_stat_value_sum_w0 += bd_stat_value[0];
            kbd_stat_value_sum_w1 += bd_stat_value[1];
            // maximum K-1 version:
            bd_stat_w0_array[s] = bd_stat_value[0];
            bd_stat_w1_array[s] = bd_stat_value[1];
            s++;
            free(ij_dst);
        }
    }
    // compute maxsum version statistic:
    s = 0;
    for (int l = 0; l < *k; ++l) {
        for (int m = (l + 1); m < *k; ++m) {
            bd_stat_w0_sum_array[l] += bd_stat_w0_array[s];
            bd_stat_w0_sum_array[m] += bd_stat_w0_array[s];
            bd_stat_w1_sum_array[l] += bd_stat_w1_array[s];
            bd_stat_w1_sum_array[m] += bd_stat_w1_array[s];
            s++;
        }
    }
    quick_sort(bd_stat_w0_sum_array, *k);
    quick_sort(bd_stat_w1_sum_array, *k);
    kbd_stat_value_maxsum_w0 = bd_stat_w0_sum_array[*k - 1];
    kbd_stat_value_maxsum_w1 = bd_stat_w1_sum_array[*k - 1];
    free(bd_stat_w0_sum_array);
    free(bd_stat_w1_sum_array);

    // compute maximum K-1 version statistic:
    quick_sort(bd_stat_w0_array, bd_stat_number);
    quick_sort(bd_stat_w1_array, bd_stat_number);
    for (i = bd_stat_number - 1; i > (bd_stat_number - *k); i--) {
        kbd_stat_value_max_w0 += bd_stat_w0_array[i];
        kbd_stat_value_max_w1 += bd_stat_w1_array[i];
    }
    free(bd_stat_w0_array);
    free(bd_stat_w1_array);

    kbd_stat[0] = kbd_stat_value_sum_w0;
    kbd_stat[1] = kbd_stat_value_sum_w1;
    kbd_stat[2] = kbd_stat_value_max_w0;
    kbd_stat[3] = kbd_stat_value_max_w1;
    kbd_stat[4] = kbd_stat_value_maxsum_w0;
    kbd_stat[5] = kbd_stat_value_maxsum_w1;
}


/**
 * Ball divergence based k-sample test
 * @param kbd : K-samples ball divergence statistic or p-value
 * @param pvalue : p-value
 * @param xy : vectorized distance matrix
 * @param size : an array contain sample size in each group
 * @param n : sample size
 * @param k : group number
 * @param R
 */
void KBD(double *kbd, double *pvalue, double *xy, int *size, int *n, int *k, int *R, int *thread) {
    int *cumsum_size;
    cumsum_size = (int *) malloc(*k * sizeof(int));
    compute_cumulate_size(cumsum_size, size, k);
    K_Ball_Divergence(kbd, xy, cumsum_size, size, n, k);

    // permutation test:
    if ((*R) > 0) {
        int j;
        double *permuted_kbd_sum_w0, *permuted_kbd_sum_w1, *permuted_kbd_max_w0, *permuted_kbd_max_w1, *permuted_kbd_maxsum_w0, *permuted_kbd_maxsum_w1;
        permuted_kbd_sum_w0 = (double *) malloc(*R * sizeof(double));
        permuted_kbd_sum_w1 = (double *) malloc(*R * sizeof(double));
        permuted_kbd_max_w0 = (double *) malloc(*R * sizeof(double));
        permuted_kbd_max_w1 = (double *) malloc(*R * sizeof(double));
        permuted_kbd_maxsum_w0 = (double *) malloc(*R * sizeof(double));
        permuted_kbd_maxsum_w1 = (double *) malloc(*R * sizeof(double));

        int not_parallel = *thread == 1 ? 1 : 0;
#ifdef Ball_OMP_H_
        if (not_parallel != 1) {
            omp_set_dynamic(0);
            if (*thread < 0) {
                omp_set_num_threads(omp_get_num_procs());
            } else {
                omp_set_num_threads(*thread);
            }
        }
#endif

        int *index = (int *) malloc(*n * sizeof(int));
        for (int i = 0; i < *n; i++) {
            index[i] = i;
        }
        if (not_parallel) {
            double *new_xy, kbd_tmp[6];
            new_xy = (double *) malloc((*n * *n) * sizeof(double));

            for (j = 0; j < (*R); j++) {
                // stop permutation if user stop it manually:
                if (pending_interrupt()) {
                    print_stop_message();
                    break;
                }
                // permute data index:
                shuffle(index, n);
                // adjust vectorized distance matrix according to permuted index:
                permute_dst(xy, new_xy, index, n);
                // K-sample BD after permutation:
                K_Ball_Divergence(kbd_tmp, new_xy, cumsum_size, size, n, k);
                permuted_kbd_sum_w0[j] = kbd_tmp[0];
                permuted_kbd_sum_w1[j] = kbd_tmp[1];
                permuted_kbd_max_w0[j] = kbd_tmp[2];
                permuted_kbd_max_w1[j] = kbd_tmp[3];
                permuted_kbd_maxsum_w0[j] = kbd_tmp[4];
                permuted_kbd_maxsum_w1[j] = kbd_tmp[5];
            }
            free(new_xy);
        } else {
            int **index_matrix = alloc_int_matrix(*R, *n);
            shuffle_indicator_matrix(index_matrix, index, *R, *n);
#pragma omp parallel
            {
                int j_thread;
                double *new_xy, kbd_tmp[6];
                new_xy = (double *) malloc((*n * *n) * sizeof(double));
#pragma omp for
                for (j_thread = 0; j_thread < (*R); j_thread++) {
                    permute_dst(xy, new_xy, index_matrix[j_thread], n);
                    K_Ball_Divergence(kbd_tmp, new_xy, cumsum_size, size, n, k);
                    permuted_kbd_sum_w0[j_thread] = kbd_tmp[0];
                    permuted_kbd_sum_w1[j_thread] = kbd_tmp[1];
                    permuted_kbd_max_w0[j_thread] = kbd_tmp[2];
                    permuted_kbd_max_w1[j_thread] = kbd_tmp[3];
                    permuted_kbd_maxsum_w0[j_thread] = kbd_tmp[4];
                    permuted_kbd_maxsum_w1[j_thread] = kbd_tmp[5];
                }
                free(new_xy);
            };
            free_int_matrix(index_matrix, *R, *n);
            j = *R;
        }

        pvalue[0] = compute_pvalue(kbd[0], permuted_kbd_sum_w0, j);
        pvalue[1] = compute_pvalue(kbd[1], permuted_kbd_sum_w1, j);
        pvalue[2] = compute_pvalue(kbd[2], permuted_kbd_max_w0, j);
        pvalue[3] = compute_pvalue(kbd[3], permuted_kbd_max_w1, j);
        pvalue[4] = compute_pvalue(kbd[4], permuted_kbd_maxsum_w0, j);
        pvalue[5] = compute_pvalue(kbd[5], permuted_kbd_maxsum_w1, j);
        free(permuted_kbd_sum_w0);
        free(permuted_kbd_sum_w1);
        free(permuted_kbd_max_w0);
        free(permuted_kbd_max_w1);
        free(permuted_kbd_maxsum_w0);
        free(permuted_kbd_maxsum_w1);
        free(index);
    }
    free(cumsum_size);
}

/**
 * @deprecated
 */
void KBD_parallel(double *kbd, double *pvalue, double *xy, int *size, int *n, int *k, int *R, int *nthread) {
    int N = (*n);
    int dst_size = N * N;
    int i_permute;
    kbd_value(kbd, xy, size, n, k);

    // permutation test:
    if ((*R) > 0) {
        double *permuted_kbd_sum_w0, *permuted_kbd_sum_w1, *permuted_kbd_max_w0, *permuted_kbd_max_w1, *permuted_kbd_max1_w0, *permuted_kbd_max1_w1;
        permuted_kbd_sum_w0 = (double *) malloc(*R * sizeof(double));
        permuted_kbd_sum_w1 = (double *) malloc(*R * sizeof(double));
        permuted_kbd_max_w0 = (double *) malloc(*R * sizeof(double));
        permuted_kbd_max_w1 = (double *) malloc(*R * sizeof(double));
        permuted_kbd_max1_w0 = (double *) malloc(*R * sizeof(double));
        permuted_kbd_max1_w1 = (double *) malloc(*R * sizeof(double));

#ifdef Ball_OMP_H_
        omp_set_num_threads(*nthread);
#endif

#pragma omp parallel
        {
            int *index, j = 0;
            double *new_xy, kbd_tmp[6];
            new_xy = (double *) malloc(dst_size * sizeof(double));
            index = (int *) malloc(N * sizeof(int));
            for (i_permute = 0; i_permute < N; i_permute++) {
                index[i_permute] = i_permute;
            }

#pragma omp for
            for (j = 0; j < (*R); j++) {
                // permute data index:
#pragma omp critical
                { shuffle(index, n); };
                // adjust vectorized distance matrix according to permuted index:
                permute_dst(xy, new_xy, index, n);
                // K-sample BD after permutation:
                kbd_value(kbd_tmp, new_xy, size, n, k);
                permuted_kbd_sum_w0[j] = kbd_tmp[0];
                permuted_kbd_sum_w1[j] = kbd_tmp[1];
                permuted_kbd_max_w0[j] = kbd_tmp[2];
                permuted_kbd_max_w1[j] = kbd_tmp[3];
                permuted_kbd_max1_w0[j] = kbd_tmp[4];
                permuted_kbd_max1_w1[j] = kbd_tmp[5];
            }
            free(new_xy);
            free(index);
        };
        pvalue[0] = compute_pvalue(kbd[0], permuted_kbd_sum_w0, *R);
        pvalue[1] = compute_pvalue(kbd[1], permuted_kbd_sum_w1, *R);
        pvalue[2] = compute_pvalue(kbd[2], permuted_kbd_max_w0, *R);
        pvalue[3] = compute_pvalue(kbd[3], permuted_kbd_max_w1, *R);
        pvalue[4] = compute_pvalue(kbd[4], permuted_kbd_max1_w0, *R);
        pvalue[5] = compute_pvalue(kbd[5], permuted_kbd_max1_w1, *R);

        free(permuted_kbd_sum_w0);
        free(permuted_kbd_sum_w1);
        free(permuted_kbd_max_w0);
        free(permuted_kbd_max_w1);
        free(permuted_kbd_max1_w0);
        free(permuted_kbd_max1_w1);

    }
}

/*
* K-sample ball divergence test, computational complexity O(n^2)
* However, it seems that we can't control the type-I error rate
*/
void KBD2(double *kbd, double *pvalue, double *xy, int *size, int *n, int *k, int *R) {
    if (*R == 0) {
        kbd_value(kbd, xy, size, n, k);
    } else {
        int two_group_size;
        int i, j, r, s = 0;
        int **Ixy, **Rxy, **Rx, *i_perm, *i_perm_tmp;
        double **Dxy, **permuted_kbd_stat_value_w0_array, **permuted_kbd_stat_value_w1_array, *ij_dst;
        double kbd_stat_value_sum_w0, kbd_stat_value_sum_w1, kbd_stat_value_max_w0, kbd_stat_value_max_w1, kbd_stat_value_maxsum_w0, kbd_stat_value_maxsum_w1;
        double *bd_stat_w0_array, *bd_stat_w1_array, *bd_stat_w0_part_sum_array, *bd_stat_w1_part_sum_array;
        double bd_stat_value[2];

        kbd_stat_value_sum_w0 = 0.0;
        kbd_stat_value_sum_w1 = 0.0;
        kbd_stat_value_max_w0 = 0.0;
        kbd_stat_value_max_w1 = 0.0;

        int *cumulate_size;
        cumulate_size = (int *) malloc(*k * sizeof(int));
        compute_cumulate_size(cumulate_size, size, k);

        int bd_stat_number = (*k) * ((*k) - 1) / 2;
        bd_stat_w0_array = (double *) malloc(bd_stat_number * sizeof(double));
        bd_stat_w1_array = (double *) malloc(bd_stat_number * sizeof(double));
        bd_stat_w0_part_sum_array = (double *) malloc(*k * sizeof(double));
        bd_stat_w1_part_sum_array = (double *) malloc(*k * sizeof(double));
        double *permuted_kbd_stat_value_sum_w0, *permuted_kbd_stat_value_sum_w1,
                *permuted_kbd_stat_value_max_w0, *permuted_kbd_stat_value_max_w1,
                *permuted_kbd_stat_value_maxsum_w0, *permuted_kbd_stat_value_maxsum_w1;
        permuted_kbd_stat_value_sum_w0 = (double *) malloc(*R * sizeof(double));
        permuted_kbd_stat_value_sum_w1 = (double *) malloc(*R * sizeof(double));
        permuted_kbd_stat_value_max_w0 = (double *) malloc(*R * sizeof(double));
        permuted_kbd_stat_value_max_w1 = (double *) malloc(*R * sizeof(double));
        permuted_kbd_stat_value_maxsum_w0 = (double *) malloc(*R * sizeof(double));
        permuted_kbd_stat_value_maxsum_w1 = (double *) malloc(*R * sizeof(double));
        for (r = 0; r < *R; r++) {
            permuted_kbd_stat_value_sum_w0[r] = 0.0;
            permuted_kbd_stat_value_sum_w1[r] = 0.0;
            permuted_kbd_stat_value_max_w0[r] = 0.0;
            permuted_kbd_stat_value_max_w1[r] = 0.0;
        }

        permuted_kbd_stat_value_w0_array = alloc_matrix(*R, bd_stat_number);
        permuted_kbd_stat_value_w1_array = alloc_matrix(*R, bd_stat_number);

        for (int group_i = 0; group_i < *k; group_i++) {
            for (int group_j = (group_i + 1); group_j < *k; group_j++) {
                int *n1 = &size[group_i];
                int *n2 = &size[group_j];
                two_group_size = *n1 + *n2;
                ij_dst = (double *) malloc((two_group_size * two_group_size) * sizeof(double));

                // allocate space
                get_ij_dst(xy, ij_dst, cumulate_size, size, n, &group_i, &group_j);
                Dxy = alloc_matrix(two_group_size, two_group_size);
                Ixy = alloc_int_matrix(two_group_size, two_group_size);
                Rx = alloc_int_matrix(two_group_size, two_group_size);
                Rxy = alloc_int_matrix(two_group_size, two_group_size);
                i_perm = (int *) malloc(two_group_size * sizeof(int));
                i_perm_tmp = (int *) malloc(two_group_size * sizeof(int));

                // Init Dxy:
                vector2matrix(ij_dst, Dxy, two_group_size, two_group_size, 1);
                // Init Ixy:
                for (i = 0; i < two_group_size; i++) {
                    for (j = 0; j < two_group_size; j++) {
                        Ixy[i][j] = j;
                    }
                }
                // Init i_perm and i_perm_tmp:
                for (i = 0; i < two_group_size; i++) {
                    i_perm[i] = i < *n1 ? 1 : 0;
                    i_perm_tmp[i] = i;
                }
                // compute two sample ball divergence:
                for (i = 0; i < two_group_size; i++) {
                    quicksort(Dxy[i], Ixy[i], 0, two_group_size - 1);
                }
                ranksort2(two_group_size, Rxy, Dxy, Ixy);
                Findx(Rxy, Ixy, i_perm, n1, n2, Rx);
                Ball_Divergence(bd_stat_value, Rxy, Rx, i_perm_tmp, n1, n2);
                kbd_stat_value_sum_w0 += bd_stat_value[0];
                kbd_stat_value_sum_w1 += bd_stat_value[1];
                bd_stat_w0_array[s] = bd_stat_value[0];
                bd_stat_w1_array[s] = bd_stat_value[1];

                for (r = 0; r < *R; r++) {
                    resample_indicator_label(i_perm, i_perm_tmp, two_group_size, n1);
                    Findx(Rxy, Ixy, i_perm, n1, n2, Rx);
                    Ball_Divergence(bd_stat_value, Rxy, Rx, i_perm_tmp, n1, n2);

                    // permuted summation statistic:
                    permuted_kbd_stat_value_sum_w0[r] += bd_stat_value[0];
                    permuted_kbd_stat_value_sum_w1[r] += bd_stat_value[1];
                    // pre-process maximum statistic:
                    permuted_kbd_stat_value_w0_array[r][s] = bd_stat_value[0];
                    permuted_kbd_stat_value_w1_array[r][s] = bd_stat_value[1];
                }

                free_matrix(Dxy, two_group_size, two_group_size);
                free_int_matrix(Ixy, two_group_size, two_group_size);
                free_int_matrix(Rxy, two_group_size, two_group_size);
                free_int_matrix(Rx, two_group_size, two_group_size);
                free(i_perm);
                free(i_perm_tmp);
                free(ij_dst);
                s += 1;
            }
        }
        // compute the maxsum statistic:
        s = 0;
        for (int l = 0; l < *k; ++l) {
            for (int m = (l + 1); m < *k; ++m) {
                bd_stat_w0_part_sum_array[l] += bd_stat_w0_array[s];
                bd_stat_w0_part_sum_array[m] += bd_stat_w0_array[s];
                bd_stat_w1_part_sum_array[l] += bd_stat_w1_array[s];
                bd_stat_w1_part_sum_array[m] += bd_stat_w1_array[s];
                s++;
            }
        }
        quick_sort(bd_stat_w0_part_sum_array, *k);
        quick_sort(bd_stat_w1_part_sum_array, *k);
        kbd_stat_value_maxsum_w0 = bd_stat_w0_part_sum_array[*k - 1];
        kbd_stat_value_maxsum_w1 = bd_stat_w1_part_sum_array[*k - 1];

        // compute the maximum statistic:
        quick_sort(bd_stat_w0_array, bd_stat_number);
        quick_sort(bd_stat_w1_array, bd_stat_number);
        for (i = bd_stat_number - 1; i > (bd_stat_number - *k); i--) {
            kbd_stat_value_max_w0 += bd_stat_w0_array[i];
            kbd_stat_value_max_w1 += bd_stat_w1_array[i];
        }

        // compute the permuted maxsum statistic:
        double **bd_stat_w0_part_sum_array2, **bd_stat_w1_part_sum_array2;
        bd_stat_w0_part_sum_array2 = alloc_matrix(*R, *k);
        bd_stat_w1_part_sum_array2 = alloc_matrix(*R, *k);
        for (r = 0; r < *R; r++) {
            s = 0;
            for (int l = 0; l < *k; ++l) {
                for (int m = (l + 1); m < *k; ++m) {
                    bd_stat_w0_part_sum_array2[r][l] += permuted_kbd_stat_value_w0_array[r][s];
                    bd_stat_w0_part_sum_array2[r][m] += permuted_kbd_stat_value_w0_array[r][s];
                    bd_stat_w1_part_sum_array2[r][l] += permuted_kbd_stat_value_w1_array[r][s];
                    bd_stat_w1_part_sum_array2[r][m] += permuted_kbd_stat_value_w1_array[r][s];
                    s++;
                }
            }
            quick_sort(bd_stat_w0_part_sum_array2[r], *k);
            quick_sort(bd_stat_w1_part_sum_array2[r], *k);
            permuted_kbd_stat_value_maxsum_w0[r] = bd_stat_w0_part_sum_array2[r][*k - 1];
            permuted_kbd_stat_value_maxsum_w0[r] = bd_stat_w0_part_sum_array2[r][*k - 1];
        }

        // compute the permuted maximum statistic:
        for (r = 0; r < *R; r++) {
            for (s = 0; s < bd_stat_number; s++) {
                bd_stat_w0_array[s] = permuted_kbd_stat_value_w0_array[r][s];
                bd_stat_w1_array[s] = permuted_kbd_stat_value_w1_array[r][s];
            }
            quick_sort(bd_stat_w0_array, bd_stat_number);
            quick_sort(bd_stat_w1_array, bd_stat_number);
            for (i = bd_stat_number - 1; i > (bd_stat_number - *k); i--) {
                permuted_kbd_stat_value_max_w0[r] += bd_stat_w0_array[i];
                permuted_kbd_stat_value_max_w1[r] += bd_stat_w1_array[i];
            }
        }

        // compute p-value:
        kbd[0] = kbd_stat_value_sum_w0;
        kbd[1] = kbd_stat_value_sum_w1;
        kbd[2] = kbd_stat_value_max_w0;
        kbd[3] = kbd_stat_value_max_w1;
        kbd[4] = kbd_stat_value_maxsum_w0;
        kbd[5] = kbd_stat_value_maxsum_w1;
        pvalue[0] = compute_pvalue(kbd[0], permuted_kbd_stat_value_sum_w0, *R);
        pvalue[1] = compute_pvalue(kbd[1], permuted_kbd_stat_value_sum_w1, *R);
        pvalue[2] = compute_pvalue(kbd[2], permuted_kbd_stat_value_max_w0, *R);
        pvalue[3] = compute_pvalue(kbd[3], permuted_kbd_stat_value_max_w1, *R);
        pvalue[4] = compute_pvalue(kbd[4], permuted_kbd_stat_value_maxsum_w0, *R);
        pvalue[5] = compute_pvalue(kbd[5], permuted_kbd_stat_value_maxsum_w1, *R);

        free_matrix(permuted_kbd_stat_value_w0_array, *R, bd_stat_number);
        free_matrix(permuted_kbd_stat_value_w1_array, *R, bd_stat_number);
        free_matrix(bd_stat_w0_part_sum_array2, *R, *k);
        free_matrix(bd_stat_w1_part_sum_array2, *R, *k);
        free(permuted_kbd_stat_value_sum_w0);
        free(permuted_kbd_stat_value_sum_w1);
        free(permuted_kbd_stat_value_max_w0);
        free(permuted_kbd_stat_value_max_w1);
        free(permuted_kbd_stat_value_maxsum_w0);
        free(permuted_kbd_stat_value_maxsum_w1);
        free(bd_stat_w0_array);
        free(bd_stat_w1_array);
        free(bd_stat_w0_part_sum_array);
        free(bd_stat_w1_part_sum_array);
        free(cumulate_size);
    }
}


/**
 * return ball divergence value of two sample(univariate) with sample size n1, n2;
 * permutation procedure not consider at all
 * input:
 * xy: value of sample
 * n1, n2: sample sizes of group1 and group2
 */
void ubd_value(double *bd_stat, double *xy, int *n1, int *n2) {
    int i, j, n;
    int **Ixy, **Rxy, **Rx, *i_perm, *xy_idx;
    double p1, p2, p3, ans, TS_weight0 = 0, SS_weight0 = 0, TS_weight1 = 0.0, SS_weight1 = 0.0;

    n = *n1 + *n2;
    Ixy = alloc_int_matrix(n, n);
    Rx = alloc_int_matrix(n, n);
    Rxy = alloc_int_matrix(n, n);
    i_perm = (int *) malloc(n * sizeof(int));
    xy_idx = (int *) malloc(n * sizeof(int));
    for (i = 0; i < n; i++) {
        xy_idx[i] = i;
        for (j = 0; j < n; j++) {
            Ixy[i][j] = j;
        }
    }

    for (i = 0; i < n; i++) {
        i_perm[i] = i < (*n1) ? 1 : 0;
    }

    quicksort(xy, xy_idx, 0, n - 1);
    ranksort3(n, xy_idx, xy, Rxy, Ixy);
    Findx(Rxy, Ixy, i_perm, n1, n2, Rx);
    free_int_matrix(Ixy, n, n);
    free(xy_idx);
    free(i_perm);

    // Ball Divergence statistic:
    // Calculate A_{ij}^{X} and A_{ij}^{Y}:
    for (i = 0; i < *n1; i++) {
        for (j = 0; j < *n1; j++) {
            p1 = Rx[i][j] + 1;
            p2 = Rxy[i][j] - p1 + 1;
            p3 = (p1 + p2) / n;
            if (p3 * (1 - p3) == 0) { continue; }
            ans = p1 / (*n1) - p2 / (*n2);
            TS_weight0 += (ans * ans);
            // TODO: Weight Ball Divergence:
            TS_weight1 += (ans * ans) * 1.0;
        }
    }
    // Calculate C_{kl}^{X} and C_{kl}^{Y}:
    for (i = *n1; i < n; i++) {
        for (j = *n1; j < n; j++) {
            p1 = Rx[i][j] + 1;
            p2 = Rxy[i][j] - p1 + 1;
            p3 = (p1 + p2) / n;
            if (p3 * (1 - p3) == 0) { continue; }
            ans = p1 / (*n1) - p2 / (*n2);
            SS_weight0 += (ans * ans);
            // TODO: Weight Ball Divergence:
            SS_weight1 += (ans * ans) * 1.0;
        }
    }
    bd_stat[0] = TS_weight0 / (1.0 * (*n1) * (*n1)) + SS_weight0 / (1.0 * (*n2) * (*n2));
    bd_stat[1] = TS_weight1 / (1.0 * (*n1) * (*n1)) + SS_weight1 / (1.0 * (*n2) * (*n2));

    free_int_matrix(Rxy, n, n);
    free_int_matrix(Rx, n, n);
}

/*
 * Given input data xy, fetch the i,j group data and insert into ij_value
 * input:
 * xy: value of all sample
 * ij_value: contain sample value of sample i and j
 * size: sample size of each group
 * cumulate_size: the cumulative sums of size vector
 * (for example, size: [1, 2, 3] ==> cumulate_size: [0, 1, 3])
 * p point the i group
 * q point the j group
 */
void get_ij_value(double *xy, double *ij_value, int *cumulate_size, int *size, int *p, int *q) {
    int k1 = 0, k2 = 0;
    int i = *p;
    int j = *q;
    int n1 = size[i];
    int n2 = size[j];
    // for group1:
    int start_index_i = cumulate_size[i];
    int start_index_j = cumulate_size[j];
    for (k1 = 0; k1 < n1; k1++) {
        ij_value[k1] = xy[start_index_i + k1];
    }
    // for group2:
    for (k2 = 0; k2 < n2; k2++) {
        ij_value[k1] = xy[start_index_j + k2];
        k1 += 1;
    }
}


/*
 * Given input univariate data xy and sample size of each group, compute kbd;
 * input:
 * xy: value of all sample
 * size: sample size of each group
 * cumulate_size: the cumulative sums of size vector
 * (for example, size: [1, 2, 3] ==> cumulate_size: [0, 1, 3])
 * p point the i group
 * q point the j group
 */
void ukbd_value(double *kbd_stat, double *xy, int *size, int *k) {
    int K = *k;
    int i, j, s = 0, t = 0;
    int *cumusum_size;
    int two_group_size, tmp_value1 = 0, tmp_value2 = 0;
    int *n1 = &tmp_value1;
    int *n2 = &tmp_value2;
    double *ij_value, *bd_stat_w0_array, *bd_stat_w1_array, *bd_stat_w0_part_sum_array, *bd_stat_w1_part_sum_array;
    double bd_stat_value[2];
    double kbd_stat_value_sum_w0 = 0.0, kbd_stat_value_sum_w1 = 0.0;
    double kbd_stat_value_max_w0 = 0.0, kbd_stat_value_max_w1 = 0.0;
    double kbd_stat_value_maxsum_w0 = 0.0, kbd_stat_value_maxsum_w1 = 0.0;
    int bd_stat_number = K * (K - 1) / 2;

    bd_stat_w0_array = (double *) malloc(bd_stat_number * sizeof(double));
    bd_stat_w1_array = (double *) malloc(bd_stat_number * sizeof(double));

    bd_stat_w0_part_sum_array = (double *) malloc(K * sizeof(double));
    bd_stat_w1_part_sum_array = (double *) malloc(K * sizeof(double));

    cumusum_size = (int *) malloc(K * sizeof(int));
    compute_cumulate_size(cumusum_size, size, k);

    // KBD statistic:
    for (i = 0; i < K; i++) {
        bd_stat_w0_part_sum_array[i] = 0.0;
        bd_stat_w1_part_sum_array[i] = 0.0;
        for (t = 0; t < (i - 1); t++) {
            bd_stat_w0_part_sum_array[i] += bd_stat_w0_array[(t + 1) * i - t * (t + 1) / 2];
            bd_stat_w1_part_sum_array[i] += bd_stat_w1_array[(t + 1) * i - t * (t + 1) / 2];
        }
        for (j = (i + 1); j < K; j++) {
            tmp_value1 = size[i];
            tmp_value2 = size[j];
            two_group_size = tmp_value1 + tmp_value2;
            ij_value = (double *) malloc(two_group_size * sizeof(double));
            get_ij_value(xy, ij_value, cumusum_size, size, &i, &j);
            ubd_value(bd_stat_value, ij_value, n1, n2);
            // summation version:
            kbd_stat_value_sum_w0 += bd_stat_value[0];
            kbd_stat_value_sum_w1 += bd_stat_value[1];
            // maximum K-1 version:
            bd_stat_w0_array[s] = bd_stat_value[0];
            bd_stat_w1_array[s] = bd_stat_value[1];
            s += 1;
            // maximum version:
            bd_stat_w0_part_sum_array[i] += bd_stat_value[0];
            bd_stat_w1_part_sum_array[i] += bd_stat_value[1];
            free(ij_value);
        }
    }
    quick_sort(bd_stat_w0_array, bd_stat_number);
    quick_sort(bd_stat_w1_array, bd_stat_number);
    // compute maximum K-1 version statistic:
    for (i = bd_stat_number - 1; i > (bd_stat_number - K); i--) {
        kbd_stat_value_max_w0 += bd_stat_w0_array[i];
        kbd_stat_value_max_w1 += bd_stat_w1_array[i];
    }
    free(bd_stat_w0_array);
    free(bd_stat_w1_array);
    free(cumusum_size);
    // compute maximum version statistic:
    quick_sort(bd_stat_w0_part_sum_array, K);
    quick_sort(bd_stat_w1_part_sum_array, K);
    kbd_stat_value_maxsum_w0 = bd_stat_w0_part_sum_array[K - 1];
    kbd_stat_value_maxsum_w1 = bd_stat_w1_part_sum_array[K - 1];
    free(bd_stat_w0_part_sum_array);
    free(bd_stat_w1_part_sum_array);
    //
    kbd_stat[0] = kbd_stat_value_sum_w0;
    kbd_stat[1] = kbd_stat_value_sum_w1;
    kbd_stat[2] = kbd_stat_value_max_w0;
    kbd_stat[3] = kbd_stat_value_max_w1;
    kbd_stat[4] = kbd_stat_value_maxsum_w0;
    kbd_stat[5] = kbd_stat_value_maxsum_w1;
}

void U_K_Ball_Divergence(double *kbd_stat, double *xy, int *cumsum_size, int *size, int *k) {
    int i, j, s = 0;
    double *ij_value, *bd_stat_w0_array, *bd_stat_w1_array, *bd_stat_w0_sum_array, *bd_stat_w1_sum_array;
    double bd_stat_value[2];
    double kbd_stat_value_sum_w0 = 0.0, kbd_stat_value_sum_w1 = 0.0;
    double kbd_stat_value_max_w0 = 0.0, kbd_stat_value_max_w1 = 0.0;
    double kbd_stat_value_maxsum_w0, kbd_stat_value_maxsum_w1;
    int bd_stat_number = *k * (*k - 1) / 2;

    bd_stat_w0_array = (double *) malloc(bd_stat_number * sizeof(double));
    bd_stat_w1_array = (double *) malloc(bd_stat_number * sizeof(double));

    bd_stat_w0_sum_array = (double *) malloc(*k * sizeof(double));
    bd_stat_w1_sum_array = (double *) malloc(*k * sizeof(double));
    for (int l = 0; l < *k; ++l) {
        bd_stat_w0_sum_array[l] = 0.0;
        bd_stat_w1_sum_array[l] = 0.0;
    }

    // KBD statistic:
    for (i = 0; i < *k; i++) {
        for (j = (i + 1); j < *k; j++) {
            ij_value = (double *) malloc((size[i] + size[j]) * sizeof(double));
            get_ij_value(xy, ij_value, cumsum_size, size, &i, &j);
            ubd_value(bd_stat_value, ij_value, &size[i], &size[j]);
            // summation version:
            kbd_stat_value_sum_w0 += bd_stat_value[0];
            kbd_stat_value_sum_w1 += bd_stat_value[1];
            // maximum K-1 version:
            bd_stat_w0_array[s] = bd_stat_value[0];
            bd_stat_w1_array[s] = bd_stat_value[1];
            s += 1;
            free(ij_value);
        }
    }
    // compute maxsum version statistic:
    s = 0;
    for (int l = 0; l < *k; ++l) {
        for (int m = (l + 1); m < *k; ++m) {
            bd_stat_w0_sum_array[l] += bd_stat_w0_array[s];
            bd_stat_w0_sum_array[m] += bd_stat_w0_array[s];
            bd_stat_w1_sum_array[l] += bd_stat_w1_array[s];
            bd_stat_w1_sum_array[m] += bd_stat_w1_array[s];
            s++;
        }
    }
    quick_sort(bd_stat_w0_sum_array, *k);
    quick_sort(bd_stat_w1_sum_array, *k);
    kbd_stat_value_maxsum_w0 = bd_stat_w0_sum_array[*k - 1];
    kbd_stat_value_maxsum_w1 = bd_stat_w1_sum_array[*k - 1];
    free(bd_stat_w0_sum_array);
    free(bd_stat_w1_sum_array);

    // compute maximum K-1 version statistic:
    quick_sort(bd_stat_w0_array, bd_stat_number);
    quick_sort(bd_stat_w1_array, bd_stat_number);
    for (i = bd_stat_number - 1; i > (bd_stat_number - *k); i--) {
        kbd_stat_value_max_w0 += bd_stat_w0_array[i];
        kbd_stat_value_max_w1 += bd_stat_w1_array[i];
    }
    free(bd_stat_w0_array);
    free(bd_stat_w1_array);

    //
    kbd_stat[0] = kbd_stat_value_sum_w0;
    kbd_stat[1] = kbd_stat_value_sum_w1;
    kbd_stat[2] = kbd_stat_value_max_w0;
    kbd_stat[3] = kbd_stat_value_max_w1;
    kbd_stat[4] = kbd_stat_value_maxsum_w0;
    kbd_stat[5] = kbd_stat_value_maxsum_w1;
}

/**
 * implement ball divergence based k-sample test for univariate dataset
 * @param kbd : K-samples ball divergence statistic
 * @param pvalue
 * @param xy : dataset
 * @param size : an array contain sample size in each group
 * @param n : sample size
 * @param k : group number
 * @param R : permutation number
 */
void UKBD(double *kbd, double *pvalue, double *xy, int *size, int *n, int *k, int *R, int *thread) {
    int *cumsum_size;
    cumsum_size = (int *) malloc(*k * sizeof(int));
    compute_cumulate_size(cumsum_size, size, k);
    U_K_Ball_Divergence(kbd, xy, cumsum_size, size, k);

    if ((*R) > 0) {
        int j;
        double kbd_tmp[6];
        double *permuted_kbd_sum_w0, *permuted_kbd_sum_w1, *permuted_kbd_max_w0, *permuted_kbd_max_w1, *permuted_kbd_max1_w0, *permuted_kbd_max1_w1;
        permuted_kbd_sum_w0 = (double *) malloc(*R * sizeof(double));
        permuted_kbd_sum_w1 = (double *) malloc(*R * sizeof(double));
        permuted_kbd_max_w0 = (double *) malloc(*R * sizeof(double));
        permuted_kbd_max_w1 = (double *) malloc(*R * sizeof(double));
        permuted_kbd_max1_w0 = (double *) malloc(*R * sizeof(double));
        permuted_kbd_max1_w1 = (double *) malloc(*R * sizeof(double));

        int not_parallel = *thread == 1 ? 1 : 0;
#ifdef Ball_OMP_H_
        if (not_parallel != 1) {
            omp_set_dynamic(0);
            if (*thread < 0) {
                omp_set_num_threads(omp_get_num_procs());
            } else {
                omp_set_num_threads(*thread);
            }
        }
#endif

        if (not_parallel) {
            for (j = 0; j < (*R); j++) {
                if (pending_interrupt()) {
                    print_stop_message();
                    break;
                }
                shuffle_value(xy, n);
                U_K_Ball_Divergence(kbd_tmp, xy, cumsum_size, size, k);
                permuted_kbd_sum_w0[j] = kbd_tmp[0];
                permuted_kbd_sum_w1[j] = kbd_tmp[1];
                permuted_kbd_max_w0[j] = kbd_tmp[2];
                permuted_kbd_max_w1[j] = kbd_tmp[3];
                permuted_kbd_max1_w0[j] = kbd_tmp[4];
                permuted_kbd_max1_w1[j] = kbd_tmp[5];
            }
        } else {
            double **xy_new = alloc_matrix(*R, *n);
            shuffle_value_matrix(xy_new, xy, *R, *n);
#pragma omp parallel
            {
                int j_thread;
#pragma omp for
                for (j_thread = 0; j_thread < (*R); j_thread++) {
                    U_K_Ball_Divergence(kbd_tmp, xy_new[j_thread], cumsum_size, size, k);
                    permuted_kbd_sum_w0[j_thread] = kbd_tmp[0];
                    permuted_kbd_sum_w1[j_thread] = kbd_tmp[1];
                    permuted_kbd_max_w0[j_thread] = kbd_tmp[2];
                    permuted_kbd_max_w1[j_thread] = kbd_tmp[3];
                    permuted_kbd_max1_w0[j_thread] = kbd_tmp[4];
                    permuted_kbd_max1_w1[j_thread] = kbd_tmp[5];
                }
            }
            free_matrix(xy_new, *R, *n);
            j = *R;
        }
        pvalue[0] = compute_pvalue(kbd[0], permuted_kbd_sum_w0, j);
        pvalue[1] = compute_pvalue(kbd[1], permuted_kbd_sum_w1, j);
        pvalue[2] = compute_pvalue(kbd[2], permuted_kbd_max_w0, j);
        pvalue[3] = compute_pvalue(kbd[3], permuted_kbd_max_w1, j);
        pvalue[4] = compute_pvalue(kbd[4], permuted_kbd_max1_w0, j);
        pvalue[5] = compute_pvalue(kbd[5], permuted_kbd_max1_w1, j);
        free(permuted_kbd_sum_w0);
        free(permuted_kbd_sum_w1);
        free(permuted_kbd_max_w0);
        free(permuted_kbd_max_w1);
        free(permuted_kbd_max1_w0);
        free(permuted_kbd_max1_w1);
    }
    free(cumsum_size);
}

/**
 * @deprecated
 */
void UKBD_parallel(double *kbd, double *pvalue, double *xy, int *size, int *n, int *k, int *R) {
    ukbd_value(kbd, xy, size, k);
    if ((*R) > 0) {
#pragma omp parallel
        {
            int j;
            double kbd_tmp[6];
            double *permuted_kbd_sum_w0, *permuted_kbd_sum_w1, *permuted_kbd_max_w0, *permuted_kbd_max_w1, *permuted_kbd_max1_w0, *permuted_kbd_max1_w1, *xy_new;
            permuted_kbd_sum_w0 = (double *) malloc(*R * sizeof(double));
            permuted_kbd_sum_w1 = (double *) malloc(*R * sizeof(double));
            permuted_kbd_max_w0 = (double *) malloc(*R * sizeof(double));
            permuted_kbd_max_w1 = (double *) malloc(*R * sizeof(double));
            permuted_kbd_max1_w0 = (double *) malloc(*R * sizeof(double));
            permuted_kbd_max1_w1 = (double *) malloc(*R * sizeof(double));
            xy_new = (double *) malloc(*n * sizeof(double));
            for (int i = 0; i < *n; i++) { xy_new[i] = xy[i]; }

#pragma omp for
            for (j = 0; j < (*R); j++) {
#pragma omp critical
                {
                    shuffle_value(xy_new, n);
                };
                // K-sample UBD after permutation:
                ukbd_value(kbd_tmp, xy_new, size, k);
                permuted_kbd_sum_w0[j] = kbd_tmp[0];
                permuted_kbd_sum_w1[j] = kbd_tmp[1];
                permuted_kbd_max_w0[j] = kbd_tmp[2];
                permuted_kbd_max_w1[j] = kbd_tmp[3];
                permuted_kbd_max1_w0[j] = kbd_tmp[4];
                permuted_kbd_max1_w1[j] = kbd_tmp[5];
            }
            free(xy_new);
            // compute p-value
            pvalue[0] = compute_pvalue(kbd[0], permuted_kbd_sum_w0, *R);
            pvalue[1] = compute_pvalue(kbd[1], permuted_kbd_sum_w1, *R);
            pvalue[2] = compute_pvalue(kbd[2], permuted_kbd_max_w0, *R);
            pvalue[3] = compute_pvalue(kbd[3], permuted_kbd_max_w1, *R);
            pvalue[4] = compute_pvalue(kbd[4], permuted_kbd_max1_w0, *R);
            pvalue[5] = compute_pvalue(kbd[5], permuted_kbd_max1_w1, *R);
            // free
            free(permuted_kbd_sum_w0);
            free(permuted_kbd_sum_w1);
            free(permuted_kbd_max_w0);
            free(permuted_kbd_max_w1);
            free(permuted_kbd_max1_w0);
            free(permuted_kbd_max1_w1);
        };
    }
}

/////////////////////////////////////////////////////////////////
//////////// four functions below wrap key function /////////////
// bd_stat, kbd_stat return ball divergence statistic 
// bd_stat, kbd_stat execute ball divergence statistic based test
/////////////////////////////////////////////////////////////////

/**
 * R function call this function to execute Two sample ball divergence test
 */
void bd_test(double *bd, double *pvalue, double *xy, int *size, int *n, int *k, int *dst, int *R, int *nthread) {
    int single_thread = 1;
    if ((*R) <= 200) {
        *nthread = single_thread;
    }
    if ((*k) == 2) {
        if (*dst) {
            BD(bd, pvalue, xy, &size[0], &size[1], R, nthread);
        } else {
            UBD(bd, pvalue, xy, &size[0], &size[1], R, nthread);
        }
    } else {
        if (*dst) {
            if (*n > 1000) {
                KBD2(bd, pvalue, xy, size, n, k, R);
            } else {
                KBD(bd, pvalue, xy, size, n, k, R, nthread);
            }
        } else {
            UKBD(bd, pvalue, xy, size, n, k, R, nthread);
        }
    }
}