//
// Created by JinZhu on 2019/3/16.
//

#include "gtest/gtest.h"

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

TEST(BD, two_sample_test_univariate) {
    double ball_stat_value[2], p_value[2];
    int size[2] = {10, 10};
    int nth, k = 2, n = 20, dst = 0, R = 299;
    double *xy;
    xy = (double *) malloc(20 * sizeof(double));

    memcpy(xy, X1_X2_CONTINUOUS, 20 * sizeof(double));
    nth = 1;
    bd_test(ball_stat_value, p_value, xy, size, &n, &k, &dst, &R, &nth);
    EXPECT_GE(p_value[0], 0.05);

    memcpy(xy, X1_X2_CONTINUOUS, 20 * sizeof(double));
    nth = 2;
    bd_test(ball_stat_value, p_value, xy, size, &n, &k, &dst, &R, &nth);
    EXPECT_GE(p_value[0], 0.05);

    memcpy(xy, X1_X2_DISCRETE, 20 * sizeof(double));
    nth = 1;
    bd_test(ball_stat_value, p_value, xy, size, &n, &k, &dst, &R, &nth);
    EXPECT_LE(p_value[0], 0.05);

    memcpy(xy, X1_X2_DISCRETE, 20 * sizeof(double));
    nth = 2;
    bd_test(ball_stat_value, p_value, xy, size, &n, &k, &dst, &R, &nth);
    EXPECT_LE(p_value[0], 0.05);
}

TEST(BD, two_sample_test_multivariate) {
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

TEST(BD, k_sample_test_univariate) {
    // TODO:
}

TEST(BD, k_sample_test_multivariate) {
    // TODO:
}

TEST(KBD, k_sample_test_multivariate) {
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

TEST(BCov, independence_test_univariate) {
    double ball_stat_value[3], p_value[3];
    double *x, *y;
    int nth, n = 10, R = 299, dst = 0;
    x = (double *) malloc(n * sizeof(double));
    y = (double *) malloc(n * sizeof(double));

    memcpy(x, X1_CONTINUOUS, n * sizeof(double));
    memcpy(y, X2_CONTINUOUS, n * sizeof(double));
    nth = 1;
    bcov_test(ball_stat_value, p_value, x, y, &n, &R, &dst, &nth);
    EXPECT_GE(p_value[0], 0.05);
    memcpy(x, X1_CONTINUOUS, n * sizeof(double));
    memcpy(y, X2_CONTINUOUS, n * sizeof(double));
    nth = 2;
    bcov_test(ball_stat_value, p_value, x, y, &n, &R, &dst, &nth);
    EXPECT_GE(p_value[0], 0.05);

    memcpy(x, X1_DISCRETE, n * sizeof(double));
    memcpy(y, X2_DISCRETE, n * sizeof(double));
    nth = 1;
    bcov_test(ball_stat_value, p_value, x, y, &n, &R, &dst, &nth);
    EXPECT_LE(p_value[0], 0.05);
    memcpy(x, X1_DISCRETE, n * sizeof(double));
    memcpy(y, X2_DISCRETE, n * sizeof(double));
    nth = 2;
    bcov_test(ball_stat_value, p_value, x, y, &n, &R, &dst, &nth);
    EXPECT_LE(p_value[0], 0.05);
}

TEST(BCOV, independence_test_multivariate) {
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
