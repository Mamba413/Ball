//
// Created by JinZhu on 2019/3/13.
//

#include "gtest/gtest.h"

static double ABSOLUATE_ERROR = 0.000001;

extern "C" {
#include "test_setting.h"
#include "BD.h"
#include "BI.h"
#include "bcor.h"
#include "kbcov.h"
#include "utilities.h"
#include "Ball_omp.h"
}

TEST(BD, bd_value) {
    double ball_stat_value[2], p_value[2];
    double xy[20] = {-0.26162, -1.17788, 0.07105, 0.31527, 0.26728, 1.02388, -1.33132, 1.5708, 2.05537, -0.60046,
                     0.32773, 0.74089, 0.40954, 1.95479, -1.24204, 0.80392, -0.87819, -0.87891, -0.25194, -0.28557};
    double xy_dst[400] = {0, 0.91626, 0.33267, 0.57689, 0.5289, 1.2855, 1.0697, 1.83242, 2.31699, 0.33884, 0.58935,
                          1.00251, 0.67116, 2.21641, 0.98042, 1.06554, 0.61657, 0.61729, 0.00968000000000002, 0.02395,
                          0.91626, 0, 1.24893, 1.49315, 1.44516, 2.20176, 0.15344, 2.74868, 3.23325, 0.57742, 1.50561,
                          1.91877, 1.58742, 3.13267, 0.06416, 1.9818, 0.29969, 0.29897, 0.92594, 0.89231, 0.33267,
                          1.24893, 0, 0.24422, 0.19623, 0.95283, 1.40237, 1.49975, 1.98432, 0.67151, 0.25668, 0.66984,
                          0.33849, 1.88374, 1.31309, 0.73287, 0.94924, 0.94996, 0.32299, 0.35662, 0.57689, 1.49315,
                          0.24422, 0, 0.04799, 0.70861, 1.64659, 1.25553, 1.7401, 0.91573, 0.01246, 0.42562, 0.09427,
                          1.63952, 1.55731, 0.48865, 1.19346, 1.19418, 0.56721, 0.60084, 0.5289, 1.44516, 0.19623,
                          0.04799, 0, 0.7566, 1.5986, 1.30352, 1.78809, 0.86774, 0.06045, 0.47361, 0.14226, 1.68751,
                          1.50932, 0.53664, 1.14547, 1.14619, 0.51922, 0.55285, 1.2855, 2.20176, 0.95283, 0.70861,
                          0.7566, 0, 2.3552, 0.54692, 1.03149, 1.62434, 0.69615, 0.28299, 0.61434, 0.93091, 2.26592,
                          0.21996, 1.90207, 1.90279, 1.27582, 1.30945, 1.0697, 0.15344, 1.40237, 1.64659, 1.5986,
                          2.3552, 0, 2.90212, 3.38669, 0.73086, 1.65905, 2.07221, 1.74086, 3.28611, 0.08928, 2.13524,
                          0.45313, 0.45241, 1.07938, 1.04575, 1.83242, 2.74868, 1.49975, 1.25553, 1.30352, 0.54692,
                          2.90212, 0, 0.48457, 2.17126, 1.24307, 0.82991, 1.16126, 0.38399, 2.81284, 0.76688, 2.44899,
                          2.44971, 1.82274, 1.85637, 2.31699, 3.23325, 1.98432, 1.7401, 1.78809, 1.03149, 3.38669,
                          0.48457, 0, 2.65583, 1.72764, 1.31448, 1.64583, 0.10058, 3.29741, 1.25145, 2.93356, 2.93428,
                          2.30731, 2.34094, 0.33884, 0.57742, 0.67151, 0.91573, 0.86774, 1.62434, 0.73086, 2.17126,
                          2.65583, 0, 0.92819, 1.34135, 1.01, 2.55525, 0.64158, 1.40438, 0.27773, 0.27845, 0.34852,
                          0.31489, 0.58935, 1.50561, 0.25668, 0.01246, 0.06045, 0.69615, 1.65905, 1.24307, 1.72764,
                          0.92819, 0, 0.41316, 0.08181, 1.62706, 1.56977, 0.47619, 1.20592, 1.20664, 0.57967, 0.6133,
                          1.00251, 1.91877, 0.66984, 0.42562, 0.47361, 0.28299, 2.07221, 0.82991, 1.31448, 1.34135,
                          0.41316, 0, 0.33135, 1.2139, 1.98293, 0.0630299999999999, 1.61908, 1.6198, 0.99283, 1.02646,
                          0.67116, 1.58742, 0.33849, 0.09427, 0.14226, 0.61434, 1.74086, 1.16126, 1.64583, 1.01,
                          0.08181, 0.33135, 0, 1.54525, 1.65158, 0.39438, 1.28773, 1.28845, 0.66148, 0.69511, 2.21641,
                          3.13267, 1.88374, 1.63952, 1.68751, 0.93091, 3.28611, 0.38399, 0.10058, 2.55525, 1.62706,
                          1.2139, 1.54525, 0, 3.19683, 1.15087, 2.83298, 2.8337, 2.20673, 2.24036, 0.98042, 0.06416,
                          1.31309, 1.55731, 1.50932, 2.26592, 0.08928, 2.81284, 3.29741, 0.64158, 1.56977, 1.98293,
                          1.65158, 3.19683, 0, 2.04596, 0.36385, 0.36313, 0.9901, 0.95647, 1.06554, 1.9818, 0.73287,
                          0.48865, 0.53664, 0.21996, 2.13524, 0.76688, 1.25145, 1.40438, 0.47619, 0.0630299999999999,
                          0.39438, 1.15087, 2.04596, 0, 1.68211, 1.68283, 1.05586, 1.08949, 0.61657, 0.29969, 0.94924,
                          1.19346, 1.14547, 1.90207, 0.45313, 2.44899, 2.93356, 0.27773, 1.20592, 1.61908, 1.28773,
                          2.83298, 0.36385, 1.68211, 0, 0.000719999999999943, 0.62625, 0.59262, 0.61729, 0.29897,
                          0.94996, 1.19418, 1.14619, 1.90279, 0.45241, 2.44971, 2.93428, 0.27845, 1.20664, 1.6198,
                          1.28845, 2.8337, 0.36313, 1.68283, 0.000719999999999943, 0, 0.62697, 0.59334,
                          0.00968000000000002, 0.92594, 0.32299, 0.56721, 0.51922, 1.27582, 1.07938, 1.82274, 2.30731,
                          0.34852, 0.57967, 0.99283, 0.66148, 2.20673, 0.9901, 1.05586, 0.62625, 0.62697, 0, 0.03363,
                          0.02395, 0.89231, 0.35662, 0.60084, 0.55285, 1.30945, 1.04575, 1.85637, 2.34094, 0.31489,
                          0.6133, 1.02646, 0.69511, 2.24036, 0.95647, 1.08949, 0.59262, 0.59334, 0.03363, 0};
    int size[2] = {10, 10};
    int n = 20, R = 0, K = 2, dst = 0, nth = 1;

    // Software of Chengfeng Liu output result (golden standard): 0.0274
    // univariate case:
    bd_test(ball_stat_value, p_value, xy, size, &n, &K, &dst, &R, &nth);
    printf("Univariate Ball Divergence: %f; \n", ball_stat_value[0]);
    EXPECT_NEAR(ball_stat_value[0], 0.0274, ABSOLUATE_ERROR);

    // multivariate case:
    dst = 1;
    bd_test(ball_stat_value, p_value, xy_dst, size, &n, &K, &dst, &R, &nth);
    printf("Multivariate Ball Divergence: %f; \n", ball_stat_value[0]);
    EXPECT_NEAR(ball_stat_value[0], 0.0274, ABSOLUATE_ERROR);
}

TEST(BCov, bcov_value) {
    double ball_stat_value[3], p_value[3];
    double x[10] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
    double y[10] = {11, 12, 13, 14, 15, 16, 17, 18, 19, 20};
    double x_dst[100] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 2, 1, 0, 1, 2, 3, 4, 5, 6, 7, 3, 2,
                         1, 0, 1, 2, 3, 4, 5, 6, 4, 3, 2, 1, 0, 1, 2, 3, 4, 5, 5, 4, 3, 2, 1, 0, 1, 2, 3, 4, 6, 5, 4, 3,
                         2, 1, 0, 1, 2, 3, 7, 6, 5, 4, 3, 2, 1, 0, 1, 2, 8, 7, 6, 5, 4, 3, 2, 1, 0, 1, 9, 8, 7, 6, 5, 4,
                         3, 2, 1, 0};
    double y_dst[100] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 2, 1, 0, 1, 2, 3, 4, 5, 6, 7, 3, 2,
                         1, 0, 1, 2, 3, 4, 5, 6, 4, 3, 2, 1, 0, 1, 2, 3, 4, 5, 5, 4, 3, 2, 1, 0, 1, 2, 3, 4, 6, 5, 4, 3,
                         2, 1, 0, 1, 2, 3, 7, 6, 5, 4, 3, 2, 1, 0, 1, 2, 8, 7, 6, 5, 4, 3, 2, 1, 0, 1, 9, 8, 7, 6, 5, 4,
                         3, 2, 1, 0};
    int n = 10;
    int R = 0, dst = 0, nth = 1;

    // Software of ChengFeng Liu output result (golden standard): 0.034214
    // univariate case:
    bcov_test(ball_stat_value, p_value, x, y, &n, &R, &dst, &nth);
    printf("Univariate Ball Covariance: %f; \n", ball_stat_value[0]);
    EXPECT_NEAR(ball_stat_value[0], 0.034214, ABSOLUATE_ERROR);

    // multivariate case:
    dst = 1;
    bcov_test(ball_stat_value, p_value, x_dst, y_dst, &n, &R, &dst, &nth);
    printf("Multivariate Ball Covariance: %f; \n", ball_stat_value[0]);
    EXPECT_NEAR(ball_stat_value[0], 0.034214, ABSOLUATE_ERROR);
}

TEST(BCov, bcor_value) {
    printf("\n");
    double ball_stat_value[6];
    double x[20] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1};
    double y[10] = {11, 12, 13, 14, 15, 16, 17, 18, 19, 20};
    int x_number[2] = {1, 1};
    int f_number = 2, n = 10, dst_y = 0, dst_x = 0, k = 0, p = 1, nth = 2, size = 2;

    // univariate case:
    bcor_test(ball_stat_value, y, x, x_number, &f_number, &size, &n, &p, &k, &dst_y, &dst_x, &nth);
    EXPECT_NEAR(ball_stat_value[0], 1, ABSOLUATE_ERROR);;
    EXPECT_NEAR(ball_stat_value[1], 1, ABSOLUATE_ERROR);;
    EXPECT_NEAR(ball_stat_value[2], 1, ABSOLUATE_ERROR);;

    EXPECT_NEAR(ball_stat_value[3], 1, ABSOLUATE_ERROR);;
    EXPECT_NEAR(ball_stat_value[4], 1, ABSOLUATE_ERROR);;
    EXPECT_NEAR(ball_stat_value[5], 1, ABSOLUATE_ERROR);;

    // multivariate case:
    double y_dst[100] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 2, 1, 0, 1, 2, 3, 4, 5, 6, 7, 3, 2,
                         1, 0, 1, 2, 3, 4, 5, 6, 4, 3, 2, 1, 0, 1, 2, 3, 4, 5, 5, 4, 3, 2, 1, 0, 1, 2, 3, 4, 6, 5, 4, 3,
                         2, 1, 0, 1, 2, 3, 7, 6, 5, 4, 3, 2, 1, 0, 1, 2, 8, 7, 6, 5, 4, 3, 2, 1, 0, 1, 9, 8, 7, 6, 5, 4,
                         3, 2, 1, 0};
    dst_y = 1;
    p = 0;
    bcor_test(ball_stat_value, y_dst, x, x_number, &f_number, &size, &n, &p, &k, &dst_y, &dst_x, &nth);
    EXPECT_NEAR(ball_stat_value[0], 1, ABSOLUATE_ERROR);;
    EXPECT_NEAR(ball_stat_value[1], 1, ABSOLUATE_ERROR);;
    EXPECT_NEAR(ball_stat_value[2], 1, ABSOLUATE_ERROR);;

    EXPECT_NEAR(ball_stat_value[3], 1, ABSOLUATE_ERROR);;
    EXPECT_NEAR(ball_stat_value[4], 1, ABSOLUATE_ERROR);;
    EXPECT_NEAR(ball_stat_value[5], 1, ABSOLUATE_ERROR);;

    // distance matrix case:
    double x_dst[100] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 2, 1, 0, 1, 2, 3, 4, 5, 6, 7, 3, 2,
                         1, 0, 1, 2, 3, 4, 5, 6, 4, 3, 2, 1, 0, 1, 2, 3, 4, 5, 5, 4, 3, 2, 1, 0, 1, 2, 3, 4, 6, 5, 4, 3,
                         2, 1, 0, 1, 2, 3, 7, 6, 5, 4, 3, 2, 1, 0, 1, 2, 8, 7, 6, 5, 4, 3, 2, 1, 0, 1, 9, 8, 7, 6, 5, 4,
                         3, 2, 1, 0};
    dst_x = 1;
    bcor_test(ball_stat_value, y_dst, x_dst, x_number, &f_number, &size, &n, &p, &k, &dst_y, &dst_x, &nth);
    EXPECT_NEAR(ball_stat_value[0], 1, ABSOLUATE_ERROR);;
    EXPECT_NEAR(ball_stat_value[1], 1, ABSOLUATE_ERROR);;
    EXPECT_NEAR(ball_stat_value[2], 1, ABSOLUATE_ERROR);;
}

TEST(BCov, kbcov_value) {
    double ball_stat_value[3], p_value[3];
    double dst_x[300] = {0, 0.50667, 0.39184, 0.28863, 0.55886, 0.48608, 0.18851, 0.05645, 1.74563, 0.65341, 0.50667, 0,
                         0.89851, 0.21805, 0.05218, 0.99276, 0.69519, 0.45022, 2.2523, 1.16008, 0.39184, 0.89851, 0,
                         0.68047, 0.9507, 0.09424, 0.20333, 0.44829, 1.35379, 0.26157, 0.28863, 0.21805, 0.68047, 0,
                         0.27023, 0.77471, 0.47714, 0.23217, 2.03426, 0.94204, 0.55886, 0.05218, 0.9507, 0.27023, 0,
                         1.04494, 0.74737, 0.5024, 2.30449, 1.21227, 0.48608, 0.99276, 0.09424, 0.77471, 1.04494, 0,
                         0.29757, 0.54254, 1.25955, 0.16733, 0.18851, 0.69519, 0.20333, 0.47714, 0.74737, 0.29757, 0,
                         0.24497, 1.55712, 0.4649, 0.05645, 0.45022, 0.44829, 0.23217, 0.5024, 0.54254, 0.24497, 0,
                         1.80208, 0.70987, 1.74563, 2.2523, 1.35379, 2.03426, 2.30449, 1.25955, 1.55712, 1.80208, 0,
                         1.09222, 0.65341, 1.16008, 0.26157, 0.94204, 1.21227, 0.16733, 0.4649, 0.70987, 1.09222, 0, 0,
                         2.70912, 1.79127, 0.19072, 1.49571, 1.79843, 2.58357, 1.09798, 0.68462, 1.03527, 2.70912, 0,
                         0.91784, 2.5184, 1.2134, 0.91068, 0.12555, 1.61114, 2.0245, 1.67385, 1.79127, 0.91784, 0,
                         1.60055, 0.29556, 0.00716, 0.7923, 0.69329, 1.10666, 0.75601, 0.19072, 2.5184, 1.60055, 0,
                         1.30499, 1.60771, 2.39285, 0.90726, 0.4939, 0.84455, 1.49571, 1.2134, 0.29556, 1.30499, 0,
                         0.30272, 1.08786, 0.39773, 0.8111, 0.46045, 1.79843, 0.91068, 0.00716, 1.60771, 0.30272, 0,
                         0.78514, 0.70045, 1.11382, 0.76317, 2.58357, 0.12555, 0.7923, 2.39285, 1.08786, 0.78514, 0,
                         1.48559, 1.89895, 1.5483, 1.09798, 1.61114, 0.69329, 0.90726, 0.39773, 0.70045, 1.48559, 0,
                         0.41336, 0.06271, 0.68462, 2.0245, 1.10666, 0.4939, 0.8111, 1.11382, 1.89895, 0.41336, 0,
                         0.35065, 1.03527, 1.67385, 0.75601, 0.84455, 0.46045, 0.76317, 1.5483, 0.06271, 0.35065, 0, 0,
                         1.56849, 0.49364, 1.21639, 0.8943, 0.31705, 0.23985, 0.45008, 0.94375, 0.52013, 1.56849, 0,
                         2.06213, 2.78487, 0.67419, 1.25144, 1.80833, 2.01857, 2.51224, 1.04835, 0.49364, 2.06213, 0,
                         0.72274, 1.38794, 0.81069, 0.2538, 0.04357, 0.45011, 1.01378, 1.21639, 2.78487, 0.72274, 0,
                         2.11068, 1.53344, 0.97654, 0.76631, 0.27263, 1.73652, 0.8943, 0.67419, 1.38794, 2.11068, 0,
                         0.57725, 1.13415, 1.34438, 1.83805, 0.37416, 0.31705, 1.25144, 0.81069, 1.53344, 0.57725, 0,
                         0.5569, 0.76713, 1.26081, 0.20308, 0.23985, 1.80833, 0.2538, 0.97654, 1.13415, 0.5569, 0,
                         0.21023, 0.70391, 0.75998, 0.45008, 2.01857, 0.04357, 0.76631, 1.34438, 0.76713, 0.21023, 0,
                         0.49368, 0.97021, 0.94375, 2.51224, 0.45011, 0.27263, 1.83805, 1.26081, 0.70391, 0.49368, 0,
                         1.46389, 0.52013, 1.04835, 1.01378, 1.73652, 0.37416, 0.20308, 0.75998, 0.97021, 1.46389, 0};
    int n = 10, k = 3;
    int R = 99, dst = 1, nth = 1;

//    TODO: check once again:
    // Implement of R output result (golden standard): 0.010842
    // Single-thread:
    kbcov_test(ball_stat_value, p_value, dst_x, &k, &n, &R, &dst, &nth);
    EXPECT_NEAR(ball_stat_value[0], 0.010842, ABSOLUATE_ERROR);
    EXPECT_NEAR(ball_stat_value[1], 1.14240, ABSOLUATE_ERROR);

    // Multi-thread:
    nth = 2;
    kbcov_test(ball_stat_value, p_value, dst_x, &k, &n, &R, &dst, &nth);
    EXPECT_NEAR(ball_stat_value[0], 0.010842, ABSOLUATE_ERROR);
    EXPECT_NEAR(ball_stat_value[1], 1.14240, ABSOLUATE_ERROR);
}

TEST(BD, two_sample_univariate_bd_test) {
    double ball_stat_value[2], p_value[2];
    double x[500] = {2.25175, 1.55615, 0.54524, 2.85418, 1.10683, 0.72145, 1.97587, 0.84344, -1.70603, -0.78955,
                     2.32068, 0.97631, 0.77115, 1.82555, -0.40834, 0.70672, 1.35385, 1.65415, 0.54051, 1.76525, 0.64433,
                     1.7401, 0.16161, 1.713, 1.08914, 0.39928, 1.47068, 1.04822, 0.74917, 0.2367, 1.22338, 1.93922,
                     -0.84661, 0.0847, 2.92419, 0.57285, 0.84396, -1.94022, 0.14856, 1.23254, 1.07918, 1.79337, 2.72291,
                     0.63817, 1.01389, 0.95096, 1.5309, 1.28475, 0.56053, 0.50629, 0.18289, 1.03813, 1.24354, 1.53103,
                     0.14486, 0.97716, 0.91053, 0.52578, 0.72429, 0.43408, 1.67109, 1.27183, 1.69888, 0.76495, 0.67365,
                     0.35637, 0.01117, 1.18014, -0.3262, 0.64187, -0.02459, 0.31699, 2.20765, -0.61115, 0.16935,
                     0.96341, -0.0671, 1.43816, 0.43653, 2.5803, 1.36074, 0.10524, -0.4532, 0.908, 1.71979, 2.05462,
                     0.28815, 0.63903, 0.75906, 2.10066, 1.64896, 2.14244, 1.01404, 2.49235, 2.23282, 1.34103, 0.66155,
                     2.235, -0.35812, 1.51349, 0.90529, 1.14303, 0.36791, 1.37934, -1.43776, 1.07899, 0.80347, 2.65376,
                     1.75877, 1.57104, 0.31, 2.09237, 0.38631, 1.6475, 1.99678, 0.81316, 0.57989, -0.02215, 0.99435,
                     0.40954, 2.20811, 1.479, 1.49441, 0.0701, 0.87551, 0.39109, 1.29121, 0.54331, 1.66132, 1.1505,
                     3.23168, 1.58024, 2.08576, 2.66405, 1.79259, -0.953, 1.70187, -0.11666, 0.59768, -0.06712, 1.8215,
                     0.48039, -0.75134, 0.13491, 1.53026, 0.38776, 0.54211, 2.80549, 2.10628, 1.96142, 1.22193, 3.45196,
                     0.52233, 1.14369, 2.39459, 0.12202, -0.04252, 1.64233, 1.86858, 1.62409, 1.74682, -1.05572,
                     1.42309, 0.70916, -0.75358, 0.99209, -0.88119, -0.0581, -0.76142, 1.67469, 1.88033, 1.38972,
                     1.00137, 1.3897, 1.98954, 2.12393, 0.49024, 1.65855, 0.49964, 1.92137, 2.68253, 1.55927, -0.07893,
                     0.5343, -0.91064, 1.17404, 2.10121, 0.43845, 0.73953, 0.66175, 1.7428, 0.53864, 1.02051, 2.66725,
                     0.92197, 0.22764, -0.25545, -0.17605, -0.4281, -0.24533, 0.62213, 0.91474, -0.58603, 0.36518,
                     1.19489, 0.18535, 0.08253, 0.77608, 2.37213, 0.58334, 0.3253, 1.99076, 0.49968, 1.79846, -0.45763,
                     1.63836, 1.30264, 1.92182, 2.48508, 1.13147, -0.19058, 0.88212, -1.01598, 0.77181, -0.31026,
                     0.62621, 0.16373, 0.2404, 1.06689, -0.71961, 0.45225, 2.05092, 1.7834, 0.76702, -0.41365, 0.98748,
                     -1.33853, 0.24561, 0.06651, -0.27744, -0.76812, -0.41024, -0.38618, 0.32252, 1.63112, -0.06034,
                     0.25776, -1.54519, 0.24131, 1.54197, 0.51826, 1.97117, 1.91757, -0.67011, 0.41305, -0.05486,
                     -0.35572, 2.02861, 0.23772, 1.64604, -0.29029, -0.18933, 1.12945, -0.29203, -0.12381, 0.38669,
                     -0.48716, 0.18902, 1.0314, -0.64156, -0.99785, -0.95451, -1.00394, -1.07203, 0.79471, 3.1304,
                     0.87868, 0.62735, -0.00596, 0.56395, -0.16823, 0.581, 0.6356, 0.94813, -0.92096, 1.66221, 0.16091,
                     0.48368, -0.16963, 0.30523, -1.48444, -1.43197, 0.00041, -0.96845, 0.22431, 0.42396, -1.2707,
                     0.30458, 2.58439, 1.01013, 0.02986, -0.16052, 0.73194, -0.83472, 0.03374, -0.4647, -0.49565,
                     0.73078, 1.3263, -1.43352, 1.14522, 0.87447, 0.59477, -2.91983, 0.9821, 0.04228, -0.7579, 0.72066,
                     0.63426, -1.98498, -1.84243, 0.23127, 1.36689, 0.65751, -1.32793, 0.76658, 0.91408, -0.36984,
                     0.24676, -1.32251, 1.11424, 0.28732, -1.71797, 0.68414, 0.35187, -0.38769, 0.23127, -1.53572,
                     1.3105, 0.84261, 0.84039, -0.39682, -0.59393, -0.4193, 0.78651, -1.87045, -0.05073, 1.12713,
                     -1.13818, 0.41322, -0.85635, 0.85415, 0.09091, -1.02064, -0.86351, 0.97995, -0.44872, 1.56491,
                     1.01164, -0.31884, 0.50394, 1.38312, -0.29272, 1.17322, 1.4049, -0.33293, 2.18284, -0.1389, 1.5242,
                     1.25374, -0.55504, -1.91363, 0.45471, -1.07961, -0.69629, -1.45287, -0.48276, -0.29733, -0.1064,
                     1.7312, 0.64748, -1.12057, 0.6634, -0.48684, -0.43941, 0.40728, -0.07195, 0.61082, 0.68196,
                     -2.13277, 0.30082, -0.58186, 0.29594, -0.93074, 0.95474, -0.80438, 2.12899, 0.48902, 0.81091,
                     0.4143, -0.12889, 1.50433, 1.88676, 0.89654, 1.01675, 0.00424, 1.42727, -0.72508, 1.05808,
                     -0.02322, 0.87603, 1.03567, 0.21079, -0.09184, -0.90076, -0.4944, -1.22989, 0.78523, -0.14292,
                     0.1735, -0.87049, 1.01657, -0.07105, -0.77789, -1.67587, -0.4794, -0.50866, 1.81928, -0.78227,
                     0.93179, 1.11011, 0.53856, 1.00276, -1.01718, -0.72113, -1.22834, 0.16478, -1.49958, -0.01627,
                     0.61458, 0.30104, -0.15816, -1.32077, -2.45525, 0.43569, 0.37217, -0.9515, 0.30982, -1.27753,
                     0.48857, -0.33971, 0.01417, -1.51227, -1.08933, 1.11125, 1.00337, -0.55228, -2.49924, -0.56945,
                     1.857, -0.76726, -0.5512, 0.33967, 0.45384, -0.42878, 0.70092, 0.84725, -0.08575, 0.96277,
                     -0.88731, -0.05063, 0.56627, -0.18277, 0.00336, 0.60146, -1.4978, 0.99308, -0.59198, 0.00807,
                     0.33094, 0.67669, -1.58132, -0.05644, 0.9292, -0.94089, -0.62233, -0.14153, -1.2477, -0.13975,
                     -1.28383, 0.59095, 0.42241, -1.60643, 0.10822, -0.50155, -0.3301, 1.86959, 0.30328, -0.15598,
                     -2.63078};
    int size[2] = {250, 250};
    int k = 2;
    int n = 500;
    int dst = 0, R = 199, nth = 2;
    bd_test(ball_stat_value, p_value, x, size, &n, &k, &dst, &R, &nth);

    // Software of ChengFeng Liu output result (golden standard): 0.059821
    EXPECT_NEAR(ball_stat_value[0], 0.059821, ABSOLUATE_ERROR);
    EXPECT_NEAR(ball_stat_value[1], 0.059821, ABSOLUATE_ERROR);
    EXPECT_LE(p_value[0], 0.05);
    EXPECT_LE(p_value[1], 0.05);
}

TEST(BD, two_sample_multivariate_bd_test) {
    double ball_stat_value[2], p_value[2];
    double x[400] = {0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0,
                     13.0, 1.0, 0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0,
                     11.0, 12.0, 2.0, 1.0, 0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0,
                     9.0, 10.0, 11.0, 3.0, 2.0, 1.0, 0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0,
                     7.0, 8.0, 9.0, 10.0, 4.0, 3.0, 2.0, 1.0, 0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 0.0, 1.0, 2.0, 3.0, 4.0,
                     5.0, 6.0, 7.0, 8.0, 9.0, 5.0, 4.0, 3.0, 2.0, 1.0, 0.0, 1.0, 2.0, 3.0, 4.0, 1.0, 0.0, 1.0, 2.0, 3.0,
                     4.0, 5.0, 6.0, 7.0, 8.0, 6.0, 5.0, 4.0, 3.0, 2.0, 1.0, 0.0, 1.0, 2.0, 3.0, 2.0, 1.0, 0.0, 1.0, 2.0,
                     3.0, 4.0, 5.0, 6.0, 7.0, 7.0, 6.0, 5.0, 4.0, 3.0, 2.0, 1.0, 0.0, 1.0, 2.0, 3.0, 2.0, 1.0, 0.0, 1.0,
                     2.0, 3.0, 4.0, 5.0, 6.0, 8.0, 7.0, 6.0, 5.0, 4.0, 3.0, 2.0, 1.0, 0.0, 1.0, 4.0, 3.0, 2.0, 1.0, 0.0,
                     1.0, 2.0, 3.0, 4.0, 5.0, 9.0, 8.0, 7.0, 6.0, 5.0, 4.0, 3.0, 2.0, 1.0, 0.0, 5.0, 4.0, 3.0, 2.0, 1.0,
                     0.0, 1.0, 2.0, 3.0, 4.0, 4.0, 3.0, 2.0, 1.0, 0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 0.0, 1.0, 2.0, 3.0, 4.0,
                     5.0, 6.0, 7.0, 8.0, 9.0, 5.0, 4.0, 3.0, 2.0, 1.0, 0.0, 1.0, 2.0, 3.0, 4.0, 1.0, 0.0, 1.0, 2.0, 3.0,
                     4.0, 5.0, 6.0, 7.0, 8.0, 6.0, 5.0, 4.0, 3.0, 2.0, 1.0, 0.0, 1.0, 2.0, 3.0, 2.0, 1.0, 0.0, 1.0, 2.0,
                     3.0, 4.0, 5.0, 6.0, 7.0, 7.0, 6.0, 5.0, 4.0, 3.0, 2.0, 1.0, 0.0, 1.0, 2.0, 3.0, 2.0, 1.0, 0.0, 1.0,
                     2.0, 3.0, 4.0, 5.0, 6.0, 8.0, 7.0, 6.0, 5.0, 4.0, 3.0, 2.0, 1.0, 0.0, 1.0, 4.0, 3.0, 2.0, 1.0, 0.0,
                     1.0, 2.0, 3.0, 4.0, 5.0, 9.0, 8.0, 7.0, 6.0, 5.0, 4.0, 3.0, 2.0, 1.0, 0.0, 5.0, 4.0, 3.0, 2.0, 1.0,
                     0.0, 1.0, 2.0, 3.0, 4.0, 10.0, 9.0, 8.0, 7.0, 6.0, 5.0, 4.0, 3.0, 2.0, 1.0, 6.0, 5.0, 4.0, 3.0,
                     2.0, 1.0, 0.0, 1.0, 2.0, 3.0, 11.0, 10.0, 9.0, 8.0, 7.0, 6.0, 5.0, 4.0, 3.0, 2.0, 7.0, 6.0, 5.0,
                     4.0, 3.0, 2.0, 1.0, 0.0, 1.0, 2.0, 12.0, 11.0, 10.0, 9.0, 8.0, 7.0, 6.0, 5.0, 4.0, 3.0, 8.0, 7.0,
                     6.0, 5.0, 4.0, 3.0, 2.0, 1.0, 0.0, 1.0, 13.0, 12.0, 11.0, 10.0, 9.0, 8.0, 7.0, 6.0, 5.0, 4.0, 9.0,
                     8.0, 7.0, 6.0, 5.0, 4.0, 3.0, 2.0, 1.0, 0.0};
    int size[2] = {10, 10};
    int k = 2;
    int n = 20;
    int dst = 1, R = 300, nth = 1;
    bd_test(ball_stat_value, p_value, x, size, &n, &k, &dst, &R, &nth);
    printf("Ball statistics: %f; ", ball_stat_value[0]);
    printf("p-value: %f \n", p_value[0]);
    printf("Ball statistics: %f; ", ball_stat_value[1]);
    printf("p-value: %f \n", p_value[1]);
    // Software of ChengFeng Liu output result (golden standard): 0.059821
    EXPECT_NEAR(ball_stat_value[0], 0.141, ABSOLUATE_ERROR);
    EXPECT_NEAR(ball_stat_value[1], 0.141, ABSOLUATE_ERROR);
    EXPECT_GE(p_value[0], 0.05);
    EXPECT_GE(p_value[1], 0.05);
}

TEST(BD, k_sample_univariate_permutation) {
    // TODO:
}

TEST(BD, k_sample_multivariate_permutation) {
    // TODO:
}

TEST(BCov, bcov_test) {
    double ball_stat_value[3], p_value[3];
    double *x, *y, *x_dst, *y_dst;
    int n = 30;
    x = (double *) malloc(n * sizeof(double));
    y = (double *) malloc(n * sizeof(double));
    x_dst = (double *) malloc((n * n) * sizeof(double));
    y_dst = (double *) malloc((n * n) * sizeof(double));
    memcpy(x, BCOV_X, n * sizeof(double));
    memcpy(y, BCOV_Y, n * sizeof(double));
    memcpy(x_dst, BCOV_X_DST, (n * n) * sizeof(double));
    memcpy(y_dst, BCOV_Y_DST, (n * n) * sizeof(double));

    int R = 1000, dst = 0, nth = 2;

    // univariate case:
    bcov_test(ball_stat_value, p_value, x, y, &n, &R, &dst, &nth);
    EXPECT_NEAR(ball_stat_value[0], 0.001256868, ABSOLUATE_ERROR);

    // multivariate case:
    dst = 1;
    bcov_test(ball_stat_value, p_value, x_dst, y_dst, &n, &R, &dst, &nth);
    EXPECT_NEAR(ball_stat_value[0], 0.001256868, ABSOLUATE_ERROR);

    memcpy(x, BCOV_X, n * sizeof(double));
    memcpy(y, BCOV_Y, n * sizeof(double));
    memcpy(x_dst, BCOV_X_DST, (n * n) * sizeof(double));
    memcpy(y_dst, BCOV_Y_DST, (n * n) * sizeof(double));

    nth = 1;
    dst = 0;
    // univariate case:
    bcov_test(ball_stat_value, p_value, x, y, &n, &R, &dst, &nth);
    EXPECT_NEAR(ball_stat_value[0], 0.001256868, ABSOLUATE_ERROR);

    // multivariate case:
    dst = 1;
    bcov_test(ball_stat_value, p_value, x_dst, y_dst, &n, &R, &dst, &nth);
    EXPECT_NEAR(ball_stat_value[0], 0.001256868, ABSOLUATE_ERROR);
}