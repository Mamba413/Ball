//
// Created by JinZhu on 2019/3/13.
//

#include "gtest/gtest.h"

static double ABSOLUTE_ERROR = 0.000001;

extern "C" {
#include "test_setting.h"
#include "BD.h"
#include "BI.h"
#include "bcor.h"
#include "kbcov.h"
#include "utilities.h"
#include "Ball_omp.h"
#include "stats_crude.h"
#include "kbd.h"
}

TEST(BD, univariate_bd_value) {
    double ball_stat_value[2], ball_stat_value_golden[2], p_value[2];
    double *xy, *xy_dst;
    xy = (double *) malloc(20 * sizeof(double));
    xy_dst = (double *) malloc(400 * sizeof(double));
    memcpy(xy, X1_X2_CONTINUOUS, 20 * sizeof(double));
    memcpy(xy_dst, X1_X2_CONTINUOUS_DST, 190 * sizeof(double));

    int size[2] = {10, 10};
    int n = 20, R = 0, K = 2, dst = 0, nth = 1;

    // Software of Chengfeng Liu output result (golden standard): 0.0614
    double **dx = alloc_matrix(20, 20);
    distance2matrix(xy_dst, dx, 20);
    Ball_Divergence_Crude(ball_stat_value_golden, dx, size[0], size[1]);
    EXPECT_NEAR(ball_stat_value_golden[0], 0.0614, ABSOLUTE_ERROR);

    bd_test(ball_stat_value, p_value, xy, size, &n, &K, &dst, &R, &nth);
    EXPECT_NEAR(ball_stat_value[0], ball_stat_value_golden[0], ABSOLUTE_ERROR);
}

TEST(BD, multivariate_bd_value) {
    double ball_stat_value[2], ball_stat_value_golden[2], p_value[2];

    int size[2] = {10, 10};
    int n = 20, R = 0, K = 2, dst = 1, nth = 1;

    double **dx = alloc_matrix(20, 20);
    distance2matrix(X1_X2_CONTINUOUS_DST, dx, 20);
    Ball_Divergence_Crude(ball_stat_value_golden, dx, size[0], size[1]);

    bd_test(ball_stat_value, p_value, X1_X2_CONTINUOUS_DST, size, &n, &K, &dst, &R, &nth);
    EXPECT_NEAR(ball_stat_value[0], ball_stat_value_golden[0], ABSOLUTE_ERROR);
}

TEST(BCOV, univariate_bcov_value) {
    double ball_stat_value[3], ball_stat_value_golden[3], p_value[3];
    int n = 10, R = 0, dst = 0, nth = 1;

    double **dx, **dy;
    dx = alloc_matrix(n, n);
    dy = alloc_matrix(n, n);

    // Discrete case:
    distance2matrix(X1_DISCRETE_DST, dx, n);
    distance2matrix(X2_DISCRETE_DST, dy, n);
    Ball_Covariance_Crude(ball_stat_value_golden, dx, dy, n);
    // Software of Chengfeng Liu output result (golden standard): 0.034214
    EXPECT_NEAR(ball_stat_value_golden[0], 0.034214, ABSOLUTE_ERROR);

    bcov_test(ball_stat_value, p_value, X1_DISCRETE, X2_DISCRETE, &n, &R, &dst, &nth);
    EXPECT_NEAR(ball_stat_value[0], ball_stat_value_golden[0], ABSOLUTE_ERROR);

    // Continuous case:
    distance2matrix(X1_CONTINUOUS_DST, dx, n);
    distance2matrix(X2_CONTINUOUS_DST, dy, n);
    Ball_Covariance_Crude(ball_stat_value_golden, dx, dy, n);
    // Software of Chengfeng Liu output result (golden standard): 0.00671
    EXPECT_NEAR(ball_stat_value_golden[0], 0.00671, ABSOLUTE_ERROR);

    bcov_test(ball_stat_value, p_value, X1_CONTINUOUS, X2_CONTINUOUS, &n, &R, &dst, &nth);
    EXPECT_NEAR(ball_stat_value[0], ball_stat_value_golden[0], ABSOLUTE_ERROR);
}

TEST(BCOV, multivariate_bcov_value) {
    double ball_stat_value[3], ball_stat_value_golden[3], p_value[3];
    int n = 10, R = 0, dst = 1, nth = 1;

    double **dx, **dy;
    dx = alloc_matrix(n, n);
    dy = alloc_matrix(n, n);
    distance2matrix(X1_DISCRETE_DST, dx, n);
    distance2matrix(X2_DISCRETE_DST, dy, n);
    Ball_Covariance_Crude(ball_stat_value_golden, dx, dy, n);

    bcov_test(ball_stat_value, p_value, X1_DISCRETE_DST, X2_DISCRETE_DST, &n, &R, &dst, &nth);
    EXPECT_NEAR(ball_stat_value[0], ball_stat_value_golden[0], ABSOLUTE_ERROR);

    // Continuous case:
    distance2matrix(X1_CONTINUOUS_DST, dx, n);
    distance2matrix(X2_CONTINUOUS_DST, dy, n);
    Ball_Covariance_Crude(ball_stat_value_golden, dx, dy, n);

    bcov_test(ball_stat_value, p_value, X1_CONTINUOUS_DST, X2_CONTINUOUS_DST, &n, &R, &dst, &nth);
    EXPECT_NEAR(ball_stat_value[0], ball_stat_value_golden[0], ABSOLUTE_ERROR);
}

TEST(KBD, multivariate_kbd_value) {
    double kbd_stat[6], pvalue[6];
    int n = 30, k = 3, R = 199, nth = 1, size[3] = {10, 10, 10};
    KBD3(kbd_stat, pvalue, X1_X2_X3_CONTINUOUS_DST, size, &n, &k, &R, &nth);
    EXPECT_NEAR(kbd_stat[0], 0.1372, ABSOLUTE_ERROR);
    EXPECT_NEAR(kbd_stat[2], 0.1042, ABSOLUTE_ERROR);
    EXPECT_NEAR(kbd_stat[4], 0.1042, ABSOLUTE_ERROR);
    nth = 2;
    KBD3(kbd_stat, pvalue, X1_X2_X3_CONTINUOUS_DST, size, &n, &k, &R, &nth);
    EXPECT_NEAR(kbd_stat[0], 0.1372, ABSOLUTE_ERROR);
    EXPECT_NEAR(kbd_stat[2], 0.1042, ABSOLUTE_ERROR);
    EXPECT_NEAR(kbd_stat[4], 0.1042, ABSOLUTE_ERROR);
}

TEST(KBCOV, multivariate_kbcov_value) {
    double ball_stat_value[3], ball_stat_value_golden[3], p_value[3];
    int n = 10, k = 3, R = 0, dst = 1, nth = 1, dst_num = n * (n - 1)>>1;;
    double ***dxyz = alloc_3d_matrix(n, n, k);
    double **dx = alloc_matrix(n, n);
    double **dy = alloc_matrix(n, n);
    double **dz = alloc_matrix(n, n);
    double *dxyz_vector;
    dxyz_vector = (double *) malloc((k * ((n * (n - 1))>>1)) * sizeof(double));

    // Discrete case:
    distance2matrix(X1_DISCRETE_DST, dx, n);
    distance2matrix(X2_DISCRETE_DST, dy, n);
    distance2matrix(X3_DISCRETE_DST, dz, n);
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            dxyz[i][j][0] = dx[i][j];
            dxyz[i][j][1] = dy[i][j];
            dxyz[i][j][2] = dz[i][j];
        }
    }
    K_Ball_Covariance_Crude(ball_stat_value_golden, dxyz, n, k);
    // Software of R output result (Ball_1.0.0.tar.gz, bcov function: 0.0789831
    EXPECT_NEAR(ball_stat_value_golden[0], 0.0789831, ABSOLUTE_ERROR);
    for (int l = 0; l < ((n * (n - 1))>>1); ++l) {
        dxyz_vector[l] = X1_DISCRETE_DST[l];
        dxyz_vector[l + dst_num] = X2_DISCRETE_DST[l];
        dxyz_vector[l + (dst_num << 1)] = X3_DISCRETE_DST[l];
    }
    kbcov_test(ball_stat_value, p_value, dxyz_vector, &k, &n, &R, &dst, &nth);
    EXPECT_NEAR(ball_stat_value[0], ball_stat_value_golden[0], ABSOLUTE_ERROR);

    // Continuous case:
    distance2matrix(X1_CONTINUOUS_DST, dx, n);
    distance2matrix(X2_CONTINUOUS_DST, dy, n);
    distance2matrix(X3_CONTINUOUS_DST, dz, n);
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            dxyz[i][j][0] = dx[i][j];
            dxyz[i][j][1] = dy[i][j];
            dxyz[i][j][2] = dz[i][j];
        }
    }
    K_Ball_Covariance_Crude(ball_stat_value_golden, dxyz, n, k);
    // Software of R output result (Ball_1.0.0.tar.gz, bcov function: 0.01611825
    EXPECT_NEAR(ball_stat_value_golden[0], 0.01611825, ABSOLUTE_ERROR);
    for (int l = 0; l < ((n * (n - 1))>>1); ++l) {
        dxyz_vector[l] = X1_CONTINUOUS_DST[l];
        dxyz_vector[l + dst_num] = X2_CONTINUOUS_DST[l];
        dxyz_vector[l + (dst_num << 1)] = X3_CONTINUOUS_DST[l];
    }
    kbcov_test(ball_stat_value, p_value, dxyz_vector, &k, &n, &R, &dst, &nth);
    EXPECT_NEAR(ball_stat_value[0], ball_stat_value_golden[0], ABSOLUTE_ERROR);
}

TEST(BCov, bcor_value) {
    double ball_stat_value[6];
    double x[20] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1};
    double y[10] = {11, 12, 13, 14, 15, 16, 17, 18, 19, 20};
    int x_number[2] = {1, 1};
    int f_number = 2, n = 10, dst_y = 0, dst_x = 0, k = 0, p = 1, nth = 2, size = 2;

    // univariate case:
    bcor_test(ball_stat_value, y, x, x_number, &f_number, &size, &n, &p, &k, &dst_y, &dst_x, &nth);
    EXPECT_NEAR(ball_stat_value[0], 1, ABSOLUTE_ERROR);;
    EXPECT_NEAR(ball_stat_value[1], 1, ABSOLUTE_ERROR);;
    EXPECT_NEAR(ball_stat_value[2], 1, ABSOLUTE_ERROR);;

    EXPECT_NEAR(ball_stat_value[3], 1, ABSOLUTE_ERROR);;
    EXPECT_NEAR(ball_stat_value[4], 1, ABSOLUTE_ERROR);;
    EXPECT_NEAR(ball_stat_value[5], 1, ABSOLUTE_ERROR);;

    // multivariate case:
    double y_dst[100] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 2, 1, 0, 1, 2, 3, 4, 5, 6, 7, 3, 2,
                         1, 0, 1, 2, 3, 4, 5, 6, 4, 3, 2, 1, 0, 1, 2, 3, 4, 5, 5, 4, 3, 2, 1, 0, 1, 2, 3, 4, 6, 5, 4, 3,
                         2, 1, 0, 1, 2, 3, 7, 6, 5, 4, 3, 2, 1, 0, 1, 2, 8, 7, 6, 5, 4, 3, 2, 1, 0, 1, 9, 8, 7, 6, 5, 4,
                         3, 2, 1, 0};
    dst_y = 1;
    p = 0;
    bcor_test(ball_stat_value, y_dst, x, x_number, &f_number, &size, &n, &p, &k, &dst_y, &dst_x, &nth);
    EXPECT_NEAR(ball_stat_value[0], 1, ABSOLUTE_ERROR);;
    EXPECT_NEAR(ball_stat_value[1], 1, ABSOLUTE_ERROR);;
    EXPECT_NEAR(ball_stat_value[2], 1, ABSOLUTE_ERROR);;

    EXPECT_NEAR(ball_stat_value[3], 1, ABSOLUTE_ERROR);;
    EXPECT_NEAR(ball_stat_value[4], 1, ABSOLUTE_ERROR);;
    EXPECT_NEAR(ball_stat_value[5], 1, ABSOLUTE_ERROR);;

    // distance matrix case:
    double x_dst[100] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 2, 1, 0, 1, 2, 3, 4, 5, 6, 7, 3, 2,
                         1, 0, 1, 2, 3, 4, 5, 6, 4, 3, 2, 1, 0, 1, 2, 3, 4, 5, 5, 4, 3, 2, 1, 0, 1, 2, 3, 4, 6, 5, 4, 3,
                         2, 1, 0, 1, 2, 3, 7, 6, 5, 4, 3, 2, 1, 0, 1, 2, 8, 7, 6, 5, 4, 3, 2, 1, 0, 1, 9, 8, 7, 6, 5, 4,
                         3, 2, 1, 0};
    dst_x = 1;
    bcor_test(ball_stat_value, y_dst, x_dst, x_number, &f_number, &size, &n, &p, &k, &dst_y, &dst_x, &nth);
    EXPECT_NEAR(ball_stat_value[0], 1, ABSOLUTE_ERROR);;
    EXPECT_NEAR(ball_stat_value[1], 1, ABSOLUTE_ERROR);;
    EXPECT_NEAR(ball_stat_value[2], 1, ABSOLUTE_ERROR);;
}