//
// Created by JinZhu on 2019/3/13.
//


#include "gtest/gtest.h"

extern "C" {
#include "utilities.h"
}

//TEST(utilities, quick_sort_macro) {
//    double x[5] = {5.0, 4.0, 3.0, 2.0, 1.0};
//    int x_ind[5] = {0, 1, 2, 3, 4};
//    int num = 5;
//    quick_sort_macro(x, x_ind, num);
//    double true_x[5] = {1.0, 2.0, 3.0, 4.0, 5.0};
//    int true_x_ind[5] = {4, 3, 2, 1, 0};
//    for (int i = 0; i < 5; ++i) {
//        EXPECT_EQ(true_x[i], x[i]);
//        EXPECT_EQ(true_x_ind[i], x_ind[i]);
//    }
//}

TEST(utilities, quicksort) {
    double x[9] = {5.0, 4.0, 3.0, 2.0, 1.0, 6.0, 7.0, 10.0, 8.0};
    int x_ind[9] = {0, 1, 2, 3, 4, 5, 6, 7, 8};
    quicksort(x, x_ind, 0, 8);
    double true_x[9] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 10.0};
    int true_x_ind[9] = {4, 3, 2, 1, 0, 5, 6, 8, 7};
    for (int i = 0; i < 9; ++i) {
        EXPECT_EQ(true_x[i], x[i]);
        EXPECT_EQ(true_x_ind[i], x_ind[i]);
    }
}

