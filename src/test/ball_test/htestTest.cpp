//
// Created by JinZhu on 2019/3/16.
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

TEST(BDTEST, two_sample_test_univariate) {
    double ball_stat_value[2], p_value[2];
    int size[2] = {10, 10};
    int nth, k = 2, n = 20, dst = 0, R = 299;
    double *xy;
    xy = (double *) malloc(20 * sizeof(double));
    memmove(xy, X1_X2_CONTINUOUS, 20 * sizeof(double));
    nth = 1;
    bd_test(ball_stat_value, p_value, xy, size, &n, &k, &dst, &R, &nth);
    EXPECT_GE(p_value[0], 0.05);
    free(xy);

    xy = (double *) malloc(20 * sizeof(double));
    memmove(xy, X1_X2_CONTINUOUS, 20 * sizeof(double));
    nth = 2;
    bd_test(ball_stat_value, p_value, xy, size, &n, &k, &dst, &R, &nth);
    EXPECT_GE(p_value[0], 0.05);
    free(xy);

    xy = (double *) malloc(20 * sizeof(double));
    memmove(xy, X1_X2_DISCRETE, 20 * sizeof(double));
    nth = 1;
    bd_test(ball_stat_value, p_value, xy, size, &n, &k, &dst, &R, &nth);
    EXPECT_LE(p_value[0], 0.05);
    free(xy);

    xy = (double *) malloc(20 * sizeof(double));
    memmove(xy, X1_X2_DISCRETE, 20 * sizeof(double));
    nth = 2;
    bd_test(ball_stat_value, p_value, xy, size, &n, &k, &dst, &R, &nth);
    EXPECT_LE(p_value[0], 0.05);
    free(xy);
}

TEST(BDTEST, two_sample_test_multivariate) {
    double ball_stat_value[2], p_value[2];
    int size[2] = {10, 10};
    int nth, k = 2, n = 20, dst = 1, R = 299;
    nth = 1;
    bd_test(ball_stat_value, p_value, X1_X2_CONTINUOUS_DST, size, &n, &k, &dst, &R, &nth);
    EXPECT_GE(p_value[0], 0.05);

    nth = 2;
    bd_test(ball_stat_value, p_value, X1_X2_CONTINUOUS_DST, size, &n, &k, &dst, &R, &nth);
    EXPECT_GE(p_value[0], 0.05);

    nth = 1;
    bd_test(ball_stat_value, p_value, X1_X2_DISCRETE_DST, size, &n, &k, &dst, &R, &nth);
    EXPECT_LE(p_value[0], 0.05);

    nth = 2;
    bd_test(ball_stat_value, p_value, X1_X2_DISCRETE_DST, size, &n, &k, &dst, &R, &nth);
    EXPECT_LE(p_value[0], 0.05);
}

TEST(BDTEST, k_sample_test_univariate) {
    // TODO:
}

TEST(BDTEST, k_sample_test_multivariate) {
    // TODO:
}

TEST(KBDTEST, k_sample_test_multivariate) {
    double kbd_stat[6], pvalue[6];
    int n = 30, k = 3, R = 199, nth = 1, size[3] = {10, 10, 10};
    KBD3(kbd_stat, pvalue, X1_X2_X3_CONTINUOUS_DST, size, &n, &k, &R, &nth);
    EXPECT_GE(pvalue[0], 0.05);
    EXPECT_GE(pvalue[2], 0.05);
    EXPECT_GE(pvalue[4], 0.05);
    nth = 2;
    KBD3(kbd_stat, pvalue, X1_X2_X3_CONTINUOUS_DST, size, &n, &k, &R, &nth);
    EXPECT_GE(pvalue[0], 0.05);
    EXPECT_GE(pvalue[2], 0.05);
    EXPECT_GE(pvalue[4], 0.05);
}

TEST(BCOVTEST, independence_test_univariate) {
    double ball_stat_value[3], p_value[3];
    double *x, *y;
    int nth, n = 10, R = 299, dst = 0;

    x = (double *) malloc(n * sizeof(double));
    y = (double *) malloc(n * sizeof(double));
    memmove(x, X1_CONTINUOUS, n * sizeof(double));
    memmove(y, X2_CONTINUOUS, n * sizeof(double));
    nth = 1;
    bcov_test(ball_stat_value, p_value, x, y, &n, &R, &dst, &nth);
    EXPECT_GE(p_value[0], 0.05);
    free(x);
    free(y);

    x = (double *) malloc(n * sizeof(double));
    y = (double *) malloc(n * sizeof(double));
    memmove(x, X1_CONTINUOUS, n * sizeof(double));
    memmove(y, X2_CONTINUOUS, n * sizeof(double));
    nth = 2;
    bcov_test(ball_stat_value, p_value, x, y, &n, &R, &dst, &nth);
    EXPECT_GE(p_value[0], 0.05);
    free(x);
    free(y);

    x = (double *) malloc(n * sizeof(double));
    y = (double *) malloc(n * sizeof(double));
    memmove(x, X1_DISCRETE, n * sizeof(double));
    memmove(y, X2_DISCRETE, n * sizeof(double));
    nth = 1;
    bcov_test(ball_stat_value, p_value, x, y, &n, &R, &dst, &nth);
    EXPECT_LE(p_value[0], 0.05);
    free(x);
    free(y);

    x = (double *) malloc(n * sizeof(double));
    y = (double *) malloc(n * sizeof(double));
    memmove(x, X1_DISCRETE, n * sizeof(double));
    memmove(y, X2_DISCRETE, n * sizeof(double));
    nth = 2;
    bcov_test(ball_stat_value, p_value, x, y, &n, &R, &dst, &nth);
    EXPECT_LE(p_value[0], 0.05);
    free(x);
    free(y);
}

TEST(BCOVTEST, independence_test_multivariate) {
    double ball_stat_value[3], p_value[3];
    int nth, n = 10, R = 299, dst = 1;

    // Note: there exists some test dependence, and hence, the error sometimes may occurs.
    nth = 1;
    bcov_test(ball_stat_value, p_value, X1_CONTINUOUS_DST, X2_CONTINUOUS_DST, &n, &R, &dst, &nth);
    EXPECT_GE(p_value[0], 0.05);
    nth = 2;
    bcov_test(ball_stat_value, p_value, X1_CONTINUOUS_DST, X2_CONTINUOUS_DST, &n, &R, &dst, &nth);
    EXPECT_GE(p_value[0], 0.05);

    nth = 1;
    bcov_test(ball_stat_value, p_value, X1_DISCRETE_DST, X2_DISCRETE_DST, &n, &R, &dst, &nth);
    EXPECT_LE(p_value[0], 0.05);
    nth = 2;
    bcov_test(ball_stat_value, p_value, X1_DISCRETE_DST, X2_DISCRETE_DST, &n, &R, &dst, &nth);
    EXPECT_LE(p_value[0], 0.05);
}

TEST(KBCOVTEST, mutual_independence_test_multivariate) {
    double ball_stat_value[3], p_value[3];
    int nth, n = 10, k = 3, R = 299, dst = 1, dst_num = n * (n - 1)>>1;;
    double *dxyz_vector;
    dxyz_vector = (double *) malloc((k * ((n * (n - 1))>>1)) * sizeof(double));

    // Discrete case:
    nth = 1;
    for (int l = 0; l < ((n * (n - 1))>>1); ++l) {
        dxyz_vector[l] = X1_DISCRETE_DST[l];
        dxyz_vector[l + dst_num] = X2_DISCRETE_DST[l];
        dxyz_vector[l + (dst_num << 1)] = X3_DISCRETE_DST[l];
    }
    kbcov_test(ball_stat_value, p_value, dxyz_vector, &k, &n, &R, &dst, &nth);
    EXPECT_LE(p_value[0], 0.05);
    EXPECT_LE(p_value[1], 0.05);
    EXPECT_LE(p_value[2], 0.05);
    nth = 2;
    kbcov_test(ball_stat_value, p_value, dxyz_vector, &k, &n, &R, &dst, &nth);
    EXPECT_LE(p_value[0], 0.05);
    EXPECT_LE(p_value[1], 0.05);
    EXPECT_LE(p_value[2], 0.05);

    // Continuous case:
    nth = 1;
    for (int l = 0; l < ((n * (n - 1))>>1); ++l) {
        dxyz_vector[l] = X1_CONTINUOUS_DST[l];
        dxyz_vector[l + dst_num] = X2_CONTINUOUS_DST[l];
        dxyz_vector[l + (dst_num << 1)] = X3_CONTINUOUS_DST[l];
    }
    kbcov_test(ball_stat_value, p_value, dxyz_vector, &k, &n, &R, &dst, &nth);
    EXPECT_NEAR(p_value[0], 0.05, 0.05);
    EXPECT_NEAR(p_value[1], 0.05, 0.05);
    EXPECT_NEAR(p_value[2], 0.05, 0.05);
    nth = 2;
    kbcov_test(ball_stat_value, p_value, dxyz_vector, &k, &n, &R, &dst, &nth);
    EXPECT_NEAR(p_value[0], 0.05, 0.05);
    EXPECT_NEAR(p_value[1], 0.05, 0.05);
    EXPECT_NEAR(p_value[2], 0.05, 0.05);
}

TEST(BDTEST, bd_gwas_test) {
    double ball_stat_value[4], permuted_stat_value[398], p_value[4];
    int snp[40] = {0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1,
                   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
    int nth, p = 2, n = 20, R = 199;
    int unique_k_num = 1;
    int each_k_num[10] = {3, 3};
    nth = 1;
    bd_gwas_screening(ball_stat_value, permuted_stat_value, p_value, X1_X2_CONTINUOUS_DST, snp, &n, &p,
                      &unique_k_num, each_k_num, &R, &nth);
    EXPECT_GE(p_value[1], 0.05);
    EXPECT_GE(p_value[3], 0.05);

    nth = 2;
    bd_gwas_screening(ball_stat_value, permuted_stat_value, p_value, X1_X2_CONTINUOUS_DST, snp, &n, &p,
                      &unique_k_num, each_k_num, &R, &nth);
    EXPECT_GE(p_value[1], 0.05);
    EXPECT_GE(p_value[3], 0.05);

    int snp2[40] = {0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1,
                    0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2};
    nth = 1;
    bd_gwas_screening(ball_stat_value, permuted_stat_value, p_value, X1_X2_CONTINUOUS_DST, snp2, &n, &p,
                      &unique_k_num, each_k_num, &R, &nth);
    EXPECT_GE(p_value[1], 0.05);
    EXPECT_GE(p_value[3], 0.05);

    nth = 2;
    bd_gwas_screening(ball_stat_value, permuted_stat_value, p_value, X1_X2_CONTINUOUS_DST, snp2, &n, &p,
                      &unique_k_num, each_k_num, &R, &nth);
    EXPECT_GE(p_value[1], 0.05);
    EXPECT_GE(p_value[3], 0.05);

    R = 49999;
    double permuted_stat_value2[99999];
    bd_gwas_screening(ball_stat_value, permuted_stat_value2, p_value, X1_X2_CONTINUOUS_DST, snp2, &n, &p,
                      &unique_k_num, each_k_num, &R, &nth);
    EXPECT_GE(p_value[1], 0.05);
    EXPECT_GE(p_value[3], 0.05);
}