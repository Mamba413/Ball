//
// Created by JinZhu on 2019/3/13.
//


#include "gtest/gtest.h"

static double ABSOLUTE_ERROR = 0.000001;

extern "C" {
#include "utilities.h"
#include "test_setting.h"
#include "kbd.h"
}

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

TEST(utilities, count_of_number_before_self) {
    const int num = 7;
    int yrank[num] = {5, 3, 6, 2, 4, 1, 7};
    int isource[num] = {0, 1, 2, 3, 4, 5, 6};
    int icount[num] = {0, 0, 0, 0, 0, 0, 0};
    Inversions(yrank, isource, icount, num - 1, num);
    int true_icount[num - 1] = {0, 1, 0, 3, 2, 5};
    for (int i = 0; i < (num - 1); ++i) {
        EXPECT_EQ(true_icount[i], icount[i]);
    }
}

TEST(utilities, count_smaller_number_after_self_solution) {
    const int num = 7;
    double y_vector[num] = {5, 3, 6, 2, 4, 1, 7};
    int inv_count[num] = {0, 0, 0, 0, 0, 0, 0};
    count_smaller_number_after_self_solution(y_vector, inv_count, num);
    int true_icount[num] = {4, 2, 3, 1, 1, 0, 0};
    for (int i = 0; i < num; ++i) {
        EXPECT_EQ(true_icount[i], inv_count[i]);
    }

    for (int i = 0; i < num; ++i) {
        y_vector[i] = 1.0;
        inv_count[i] = 0;
        true_icount[i] = num - 1 - i;
    }
    count_smaller_number_after_self_solution(y_vector, inv_count, num);
    for (int i = 0; i < num; ++i) {
        EXPECT_EQ(true_icount[i], inv_count[i]);
    }
}

TEST(KBD, sub_rank_finder) {
    int k = 3, n = 30, size[3] = {10, 10, 10}, cumsum_size[3] = {0, 10, 20}, label[30];
    double *x123_continuous, **distance_matrix;
    x123_continuous = (double *) malloc(n * sizeof(double));
    distance_matrix = alloc_matrix(n, n);
    for (int i = 0; i < 10; ++i) {
        x123_continuous[i] = X1_CONTINUOUS[i];
        x123_continuous[i + 10] = X2_CONTINUOUS[i];
        x123_continuous[i + 20] = X3_CONTINUOUS[i];
        label[i] = 0;
        label[i + 10] = 1;
        label[i + 20] = 2;
    }
    Euclidean_distance(x123_continuous, distance_matrix, n, 1);

    int *group_relative_location = (int *) malloc(n * sizeof(int));
    find_group_relative_location(group_relative_location, label, cumsum_size, n, k);

    int **index_matrix = alloc_int_matrix(n, n);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            index_matrix[i][j] = j;
        }
    }
    double *distance_matrix_copy;
    distance_matrix_copy = (double *) malloc(n * sizeof(double));
    for (int i = 0; i < n; i++) {
        memcpy(distance_matrix_copy, distance_matrix[i], n * sizeof(double));
        quicksort(distance_matrix_copy, index_matrix[i], 0, n - 1);
    }
    free(distance_matrix_copy);
    int ***sub_rank = alloc_int_square_matrix_list(size, k);
    sub_rank_finder(sub_rank, distance_matrix, index_matrix, label, group_relative_location, cumsum_size, n, k - 1);

    int true_rank[300] = {1, 5, 3, 10, 6, 2, 7, 9, 8, 4, 7, 1, 9, 10, 2, 8, 3, 6, 4, 5, 3, 5, 1, 10, 6, 2, 7, 9, 8, 4,
                          8, 6, 10, 1, 5, 9, 4, 2, 3, 7, 7, 2, 9, 10, 1, 8, 3, 5, 4, 6, 3, 5, 2, 10, 6, 1, 7, 9, 8, 4,
                          8, 5, 10, 7, 3, 9, 1, 4, 2, 6, 8, 5, 10, 6, 4, 9, 3, 1, 2, 7, 8, 5, 10, 7, 4, 9, 2, 3, 1, 6,
                          2, 3, 5, 10, 6, 4, 7, 9, 8, 1, 1, 6, 9, 10, 2, 8, 7, 3, 4, 5, 9, 1, 8, 10, 7, 5, 3, 6, 4, 2,
                          10, 4, 1, 8, 9, 2, 3, 7, 6, 5, 10, 5, 2, 1, 9, 3, 4, 8, 7, 6, 4, 6, 9, 10, 1, 8, 7, 2, 3, 5,
                          9, 3, 4, 10, 8, 1, 2, 7, 6, 5, 9, 3, 4, 10, 8, 2, 1, 7, 6, 5, 6, 5, 9, 10, 3, 8, 7, 1, 2, 4,
                          6, 5, 9, 10, 4, 8, 7, 2, 1, 3, 8, 2, 9, 10, 5, 7, 6, 4, 3, 1, 1, 2, 5, 10, 3, 6, 7, 9, 8, 4,
                          2, 1, 5, 10, 3, 6, 7, 9, 8, 4, 8, 7, 1, 10, 5, 2, 3, 9, 6, 4, 10, 9, 6, 1, 8, 5, 4, 2, 3, 7,
                          4, 2, 5, 10, 1, 6, 7, 9, 8, 3, 8, 7, 3, 10, 6, 1, 2, 9, 4, 5, 8, 7, 3, 10, 6, 2, 1, 9, 4, 5,
                          10, 9, 6, 2, 8, 5, 4, 1, 3, 7, 9, 8, 4, 10, 7, 3, 2, 6, 1, 5, 6, 4, 3, 10, 2, 5, 7, 9, 8, 1};
    int s = 0;
    for (int l = 0; l < 3; ++l) {
        for (int i = 0; i < 10; ++i) {
            for (int j = 0; j < 10; ++j) {
                EXPECT_EQ(true_rank[s++], sub_rank[l][i][j]);
            }
        }
    }

    int label_new[30] = {0, 1, 1, 2, 0, 2, 2, 1, 1, 0, 0, 0, 1, 0, 1, 2, 2, 2, 1, 0, 0, 2, 2, 0, 2, 1, 2, 1, 0, 1};
    find_group_relative_location(group_relative_location, label_new, cumsum_size, n, k);
    sub_rank_finder(sub_rank, distance_matrix, index_matrix, label_new, group_relative_location, cumsum_size, n, k - 1);

    int true_rank_new[300] = {1, 4, 3, 10, 5, 9, 6, 8, 7, 2, 7, 1, 5, 8, 2, 10, 3, 4, 9, 6, 3, 4, 1, 9, 5, 10, 6, 7, 8,
                              2, 8, 5, 6, 1, 4, 10, 3, 2, 9, 7, 7, 2, 5, 8, 1, 10, 3, 4, 9, 6, 3, 6, 5, 10, 7, 1, 8, 9,
                              2, 4, 8, 3, 5, 6, 2, 10, 1, 4, 9, 7, 8, 4, 6, 5, 3, 10, 2, 1, 9, 7, 3, 6, 5, 10, 7, 2, 8,
                              9, 1, 4, 2, 4, 3, 10, 5, 9, 6, 7, 8, 1, 1, 9, 5, 4, 7, 8, 6, 3, 10, 2, 5, 1, 8, 7, 2, 10,
                              9, 4, 3, 6, 6, 9, 1, 3, 8, 5, 2, 7, 10, 4, 5, 9, 3, 1, 8, 6, 4, 7, 10, 2, 4, 2, 8, 7, 1,
                              10, 9, 3, 5, 6, 6, 9, 3, 4, 8, 1, 2, 7, 10, 5, 6, 9, 2, 3, 8, 4, 1, 7, 10, 5, 2, 6, 7, 5,
                              4, 9, 8, 1, 10, 3, 5, 2, 8, 7, 3, 10, 9, 4, 1, 6, 3, 9, 4, 2, 8, 7, 5, 6, 10, 1, 1, 10, 5,
                              8, 7, 2, 3, 6, 4, 9, 10, 1, 6, 3, 4, 9, 8, 5, 7, 2, 9, 10, 1, 7, 6, 5, 3, 4, 2, 8, 10, 7,
                              5, 1, 2, 9, 8, 4, 6, 3, 10, 8, 5, 2, 1, 9, 7, 3, 6, 4, 5, 10, 4, 8, 7, 1, 2, 6, 3, 9, 7,
                              10, 4, 8, 6, 2, 1, 5, 3, 9, 10, 9, 5, 3, 2, 8, 7, 1, 6, 4, 9, 10, 2, 7, 6, 4, 3, 5, 1, 8,
                              10, 6, 5, 2, 3, 9, 8, 4, 7, 1};
    s = 0;
    for (int l = 0; l < 3; ++l) {
        for (int i = 0; i < 10; ++i) {
            for (int j = 0; j < 10; ++j) {
                EXPECT_EQ(true_rank_new[s++], sub_rank[l][i][j]);
            }
        }
    }

    free(x123_continuous);
    free_int_square_matrix_list(sub_rank, size, k);
    free_int_matrix(index_matrix, n, n);
    free_matrix(distance_matrix, n, n);
    free(group_relative_location);
}

TEST(KBD, sub_rank_finder_imbalanced_multiple_group) {
    int s, k = 5, n = 30, size[5] = {6, 5, 5, 6, 8}, cumsum_size[5] = {0, 6, 11, 16, 22}, label[30];
    double **distance_matrix;
    distance_matrix = alloc_matrix(n, n);
    Euclidean_distance(X1_X2_X3_CONTINUOUS, distance_matrix, n, 1);

    s = 0;
    for (int i = 0; i < k; ++i) {
        for (int j = 0; j < size[i]; ++j) {
            label[s++] = i;
        }
    }

    int *group_relative_location = (int *) malloc(n * sizeof(int));
    find_group_relative_location(group_relative_location, label, cumsum_size, n, k);

    int **index_matrix = alloc_int_matrix(n, n);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            index_matrix[i][j] = j;
        }
    }
    double *distance_matrix_copy;
    distance_matrix_copy = (double *) malloc(n * sizeof(double));
    for (int i = 0; i < n; i++) {
        memcpy(distance_matrix_copy, distance_matrix[i], n * sizeof(double));
        quicksort(distance_matrix_copy, index_matrix[i], 0, n - 1);
    }
    free(distance_matrix_copy);

    int ***sub_rank = alloc_int_square_matrix_list(size, k);
    sub_rank_finder(sub_rank, distance_matrix, index_matrix, label, group_relative_location, cumsum_size, n, k - 1);

    int true_rank[] = {1, 4, 3, 6, 5, 2, 3, 1, 5, 6, 2, 4, 3, 4, 1, 6, 5, 2, 4, 3, 6, 1, 2, 5, 3, 2, 5, 6, 1, 4, 3, 4,
                       2, 6, 5, 1, 1, 3, 2, 4, 5, 3, 1, 2, 5, 4, 2, 3, 1, 4, 5, 2, 4, 3, 1, 5, 4, 2, 3, 5, 1, 1, 4, 5,
                       3, 2, 3, 1, 4, 5, 2, 4, 2, 1, 5, 3, 2, 4, 5, 1, 3, 2, 3, 5, 4, 1, 1, 6, 4, 2, 5, 3, 6, 1, 3, 5,
                       2, 4, 6, 4, 1, 5, 3, 2, 6, 5, 3, 1, 4, 2, 6, 2, 3, 5, 1, 4, 6, 4, 2, 5, 3, 1, 1, 8, 5, 2, 3, 7,
                       6, 4, 6, 1, 8, 5, 4, 2, 3, 7, 3, 8, 1, 4, 5, 7, 6, 2, 3, 8, 6, 1, 2, 7, 4, 5, 3, 8, 6, 2, 1, 7,
                       4, 5, 6, 2, 8, 5, 4, 1, 3, 7, 4, 8, 7, 3, 2, 6, 1, 5, 3, 8, 2, 4, 5, 7, 6, 1};
    s = 0;
    for (int l = 0; l < 5; ++l) {
        for (int i = 0; i < size[l]; ++i) {
            for (int j = 0; j < size[l]; ++j) {
                EXPECT_EQ(true_rank[s++], sub_rank[l][i][j]);
            }
        }
    }
}

TEST(KBD, full_rank_finder) {
    int s, k = 3, n = 30, size[3] = {10, 10, 10}, cumsum_size[3] = {0, 10, 20}, label[30];
    int bd_stat_number = (((k - 1) * (k)) >> 1);
    double **distance_matrix;
    distance_matrix = alloc_matrix(n, n);
    Euclidean_distance(X1_X2_X3_CONTINUOUS, distance_matrix, n, 1);

    s = 0;
    for (int i = 0; i < k; ++i) {
        for (int j = 0; j < size[i]; ++j) {
            label[s++] = i;
        }
    }

    int *group_relative_location = (int *) malloc(n * sizeof(int));
    find_group_relative_location(group_relative_location, label, cumsum_size, n, k);

    int **index_matrix = alloc_int_matrix(n, n);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            index_matrix[i][j] = j;
        }
    }
    double *distance_matrix_copy;
    distance_matrix_copy = (double *) malloc(n * sizeof(double));
    for (int i = 0; i < n; i++) {
        memcpy(distance_matrix_copy, distance_matrix[i], n * sizeof(double));
        quicksort(distance_matrix_copy, index_matrix[i], 0, n - 1);
    }
    free(distance_matrix_copy);

    int *pairwise_size = (int *) malloc(bd_stat_number * sizeof(int));
    s = 0;
    for (int i = 0; i < (k - 1); ++i) {
        for (int j = (i + 1); j <= (k - 1); ++j) {
            pairwise_size[s] = size[i] + size[j];
            s++;
        }
    }

    int ***full_rank = alloc_int_square_matrix_list(pairwise_size, bd_stat_number);
    full_rank_finder(full_rank, distance_matrix, index_matrix, label, group_relative_location,
                     cumsum_size, size, n, k - 1);
    int full_rank_vector1[400] = {1, 8, 4, 20, 9, 3, 11, 14, 12, 5, 19, 10, 2, 17, 18, 6, 7, 16, 15, 13, 14, 1, 17, 19,
                                  2, 16, 6, 10, 7, 9, 18, 4, 13, 20, 15, 5, 3, 12, 11, 8, 3, 8, 1, 20, 9, 2, 11, 15, 13,
                                  5, 19, 10, 4, 12, 18, 6, 7, 17, 16, 14, 17, 12, 19, 1, 11, 18, 9, 6, 8, 15, 2, 10, 16,
                                  20, 3, 14, 13, 4, 5, 7, 15, 3, 17, 19, 1, 16, 4, 9, 5, 12, 18, 2, 14, 20, 13, 8, 7,
                                  11, 10, 6, 3, 8, 2, 20, 9, 1, 11, 15, 13, 5, 19, 10, 4, 12, 18, 6, 7, 17, 16, 14, 17,
                                  7, 19, 15, 5, 18, 1, 6, 2, 13, 14, 3, 16, 20, 12, 11, 10, 9, 8, 4, 17, 10, 19, 14, 9,
                                  18, 6, 1, 4, 15, 12, 7, 16, 20, 8, 13, 11, 5, 2, 3, 17, 9, 19, 15, 7, 18, 3, 4, 1, 13,
                                  14, 5, 16, 20, 10, 12, 11, 8, 6, 2, 5, 6, 8, 19, 9, 7, 11, 14, 12, 1, 18, 10, 4, 20,
                                  17, 2, 3, 16, 15, 13, 17, 12, 19, 2, 11, 18, 9, 6, 8, 15, 1, 10, 16, 20, 3, 14, 13, 4,
                                  5, 7, 15, 6, 19, 17, 2, 18, 3, 7, 4, 12, 16, 1, 14, 20, 13, 10, 8, 11, 9, 5, 2, 8, 4,
                                  20, 9, 3, 11, 14, 12, 5, 19, 10, 1, 17, 18, 6, 7, 16, 15, 13, 4, 9, 2, 20, 10, 3, 12,
                                  15, 13, 6, 19, 11, 5, 1, 18, 7, 8, 17, 16, 14, 17, 12, 19, 6, 11, 18, 9, 4, 8, 15, 5,
                                  10, 16, 20, 1, 14, 13, 2, 3, 7, 9, 3, 14, 19, 5, 12, 7, 13, 10, 4, 18, 6, 8, 20, 17,
                                  1, 2, 16, 15, 11, 11, 3, 14, 19, 5, 13, 7, 12, 8, 4, 18, 6, 9, 20, 17, 2, 1, 16, 15,
                                  10, 17, 12, 19, 11, 10, 18, 7, 4, 6, 15, 9, 8, 16, 20, 3, 14, 13, 1, 2, 5, 17, 10, 19,
                                  12, 9, 18, 7, 2, 5, 15, 11, 8, 16, 20, 6, 14, 13, 3, 1, 4, 17, 9, 19, 15, 7, 18, 3, 4,
                                  2, 13, 14, 5, 16, 20, 10, 12, 11, 8, 6, 1};
    int full_rank_vector2[400] = {1, 9, 4, 20, 11, 3, 13, 17, 14, 5, 19, 18, 8, 16, 15, 7, 6, 10, 2, 12, 15, 1, 17, 18,
                                  3, 16, 6, 11, 8, 10, 14, 12, 2, 20, 9, 5, 7, 19, 13, 4, 3, 10, 1, 20, 12, 2, 14, 17,
                                  15, 5, 19, 18, 9, 11, 16, 8, 7, 6, 4, 13, 16, 10, 18, 1, 9, 17, 7, 4, 6, 14, 2, 3, 11,
                                  20, 5, 12, 13, 19, 15, 8, 15, 3, 17, 18, 1, 16, 4, 9, 5, 13, 12, 10, 6, 20, 7, 8, 11,
                                  19, 14, 2, 3, 10, 2, 20, 11, 1, 14, 17, 15, 5, 19, 18, 9, 12, 16, 8, 7, 6, 4, 13, 16,
                                  8, 18, 15, 5, 17, 1, 6, 3, 13, 10, 7, 9, 20, 4, 11, 12, 19, 14, 2, 16, 9, 18, 12, 8,
                                  17, 6, 1, 4, 14, 5, 2, 10, 20, 3, 11, 13, 19, 15, 7, 16, 9, 18, 14, 7, 17, 3, 5, 1,
                                  13, 8, 6, 10, 20, 2, 11, 12, 19, 15, 4, 5, 7, 9, 20, 10, 8, 12, 15, 13, 1, 18, 16, 6,
                                  19, 14, 4, 2, 17, 3, 11, 16, 10, 18, 9, 8, 17, 6, 3, 5, 14, 1, 2, 11, 20, 4, 12, 13,
                                  19, 15, 7, 16, 9, 18, 11, 8, 17, 6, 2, 5, 14, 3, 1, 10, 20, 4, 12, 13, 19, 15, 7, 13,
                                  2, 17, 18, 5, 16, 8, 12, 9, 7, 15, 14, 1, 20, 10, 3, 4, 19, 11, 6, 5, 11, 3, 20, 12,
                                  4, 14, 17, 15, 7, 19, 18, 10, 1, 16, 9, 8, 2, 6, 13, 16, 9, 18, 14, 7, 17, 4, 3, 2,
                                  13, 8, 5, 10, 20, 1, 11, 12, 19, 15, 6, 10, 4, 14, 19, 6, 13, 9, 15, 11, 5, 17, 16, 3,
                                  20, 12, 1, 2, 18, 7, 8, 7, 6, 12, 19, 8, 11, 10, 15, 13, 3, 17, 16, 4, 20, 14, 2, 1,
                                  18, 5, 9, 5, 11, 3, 20, 12, 4, 14, 17, 15, 7, 19, 18, 10, 2, 16, 9, 8, 1, 6, 13, 2, 9,
                                  6, 20, 10, 5, 12, 16, 14, 3, 18, 17, 8, 19, 15, 7, 4, 13, 1, 11, 15, 6, 18, 16, 3, 17,
                                  2, 7, 4, 13, 11, 9, 8, 20, 5, 10, 12, 19, 14, 1};
    int full_rank_vector3[400] = {1, 10, 17, 20, 2, 13, 12, 3, 5, 8, 4, 6, 11, 19, 7, 14, 15, 18, 16, 9, 17, 1, 16, 20,
                                  14, 9, 7, 13, 8, 3, 11, 6, 5, 19, 4, 10, 12, 18, 15, 2, 20, 9, 1, 18, 19, 5, 6, 17,
                                  15, 11, 16, 14, 7, 13, 12, 4, 3, 8, 2, 10, 20, 11, 4, 1, 19, 8, 9, 18, 16, 13, 17, 15,
                                  10, 2, 14, 7, 6, 3, 5, 12, 6, 10, 17, 20, 1, 13, 12, 2, 4, 8, 3, 5, 11, 19, 7, 14, 15,
                                  18, 16, 9, 18, 7, 9, 20, 16, 1, 3, 15, 13, 10, 14, 12, 5, 19, 11, 2, 4, 17, 6, 8, 18,
                                  6, 9, 20, 16, 2, 1, 15, 13, 10, 14, 12, 4, 19, 11, 3, 5, 17, 8, 7, 10, 9, 17, 20, 5,
                                  13, 12, 1, 3, 7, 2, 4, 11, 19, 6, 14, 15, 18, 16, 8, 10, 9, 17, 20, 7, 13, 12, 4, 1,
                                  6, 3, 2, 11, 19, 5, 14, 15, 18, 16, 8, 15, 5, 17, 20, 10, 12, 11, 8, 6, 1, 7, 4, 9,
                                  19, 2, 13, 14, 18, 16, 3, 10, 9, 17, 20, 5, 13, 12, 2, 3, 7, 1, 4, 11, 19, 6, 14, 15,
                                  18, 16, 8, 11, 9, 17, 20, 7, 13, 12, 4, 2, 6, 3, 1, 10, 19, 5, 14, 15, 18, 16, 8, 17,
                                  6, 11, 20, 16, 3, 2, 15, 13, 8, 14, 12, 1, 19, 9, 4, 5, 18, 10, 7, 20, 11, 4, 2, 19,
                                  8, 9, 18, 16, 13, 17, 15, 10, 1, 14, 7, 6, 3, 5, 12, 15, 6, 17, 20, 9, 12, 11, 8, 4,
                                  2, 7, 3, 10, 19, 1, 13, 14, 18, 16, 5, 18, 7, 9, 20, 16, 2, 3, 15, 13, 10, 14, 12, 5,
                                  19, 11, 1, 4, 17, 6, 8, 18, 8, 7, 20, 16, 3, 4, 15, 13, 10, 14, 12, 5, 19, 11, 2, 1,
                                  17, 6, 9, 20, 11, 4, 3, 19, 8, 9, 18, 16, 13, 17, 15, 10, 2, 14, 7, 6, 1, 5, 12, 20,
                                  8, 2, 19, 18, 5, 6, 16, 14, 11, 15, 13, 7, 17, 12, 4, 3, 10, 1, 9, 17, 2, 16, 20, 14,
                                  9, 8, 12, 7, 3, 11, 6, 5, 19, 4, 10, 13, 18, 15, 1};

    s = 0;
    for (int i = 0; i < 20; ++i) {
        for (int j = 0; j < 20; ++j) {
            EXPECT_EQ(full_rank_vector1[s++], full_rank[0][i][j]);
        }
    }
    s = 0;
    for (int i = 0; i < 20; ++i) {
        for (int j = 0; j < 20; ++j) {
            EXPECT_EQ(full_rank_vector2[s++], full_rank[1][i][j]);
        }
    }
    s = 0;
    for (int i = 0; i < 20; ++i) {
        for (int j = 0; j < 20; ++j) {
            EXPECT_EQ(full_rank_vector3[s++], full_rank[2][i][j]);
        }
    }

    int label_new[30] = {0, 1, 1, 2, 0, 2, 2, 1, 1, 0, 0, 0, 1, 0, 1, 2, 2, 2, 1, 0, 0, 2, 2, 0, 2, 1, 2, 1, 0, 1};
    find_group_relative_location(group_relative_location, label_new, cumsum_size, n, k);
    full_rank_finder(full_rank, distance_matrix, index_matrix, label_new,
                     group_relative_location, cumsum_size, size, n, k - 1);

    int full_rank_vector_new1[400] = {1, 9, 5, 20, 10, 18, 13, 17, 14, 3, 7, 4, 15, 12, 2, 19, 16, 6, 8, 11, 15, 1, 11,
                                      17, 2, 20, 6, 10, 19, 13, 4, 16, 8, 5, 14, 12, 9, 7, 18, 3, 5, 8, 1, 19, 9, 20,
                                      12, 16, 18, 2, 6, 7, 13, 11, 4, 17, 14, 3, 15, 10, 16, 10, 13, 1, 9, 20, 6, 3, 19,
                                      14, 11, 17, 5, 7, 15, 2, 4, 12, 18, 8, 15, 3, 11, 16, 1, 20, 5, 10, 19, 13, 6, 17,
                                      7, 4, 14, 12, 8, 9, 18, 2, 5, 11, 8, 20, 12, 1, 15, 18, 2, 7, 10, 4, 16, 14, 6,
                                      19, 17, 9, 3, 13, 16, 7, 12, 13, 5, 20, 1, 8, 19, 14, 9, 17, 3, 2, 15, 10, 6, 11,
                                      18, 4, 16, 9, 13, 10, 8, 20, 5, 1, 19, 14, 11, 17, 3, 6, 15, 4, 2, 12, 18, 7, 5,
                                      11, 8, 20, 12, 2, 15, 18, 1, 7, 10, 4, 16, 14, 6, 19, 17, 9, 3, 13, 3, 8, 4, 20,
                                      9, 19, 13, 16, 17, 1, 7, 5, 14, 12, 2, 18, 15, 6, 11, 10, 14, 2, 8, 17, 3, 20, 7,
                                      12, 19, 11, 1, 16, 9, 6, 13, 15, 10, 5, 18, 4, 2, 10, 5, 20, 11, 13, 15, 18, 9, 4,
                                      8, 1, 16, 14, 3, 19, 17, 7, 6, 12, 16, 9, 13, 11, 7, 20, 3, 5, 19, 14, 10, 17, 1,
                                      4, 15, 8, 2, 12, 18, 6, 16, 7, 12, 13, 5, 20, 2, 8, 19, 14, 9, 17, 4, 1, 15, 10,
                                      6, 11, 18, 3, 2, 9, 5, 20, 10, 18, 13, 17, 15, 3, 7, 4, 14, 12, 1, 19, 16, 6, 8,
                                      11, 16, 10, 13, 5, 9, 20, 6, 2, 19, 14, 11, 17, 4, 7, 15, 1, 3, 12, 18, 8, 16, 9,
                                      13, 11, 8, 20, 4, 3, 19, 14, 10, 17, 2, 5, 15, 6, 1, 12, 18, 7, 9, 4, 3, 18, 6,
                                      20, 11, 15, 19, 5, 2, 12, 13, 10, 8, 16, 14, 1, 17, 7, 5, 11, 8, 20, 12, 4, 15,
                                      18, 2, 7, 10, 3, 16, 14, 6, 19, 17, 9, 1, 13, 15, 3, 12, 16, 2, 20, 5, 10, 19, 13,
                                      6, 17, 7, 4, 14, 11, 8, 9, 18, 1};

    int full_rank_vector_new2[400] = {1, 9, 4, 19, 10, 18, 12, 16, 14, 2, 20, 3, 11, 6, 7, 17, 15, 8, 13, 5, 15, 1, 13,
                                      17, 2, 20, 5, 11, 19, 14, 18, 16, 3, 8, 7, 12, 9, 4, 6, 10, 6, 9, 1, 18, 10, 20,
                                      12, 15, 17, 3, 19, 8, 11, 4, 5, 16, 14, 7, 13, 2, 17, 10, 15, 1, 9, 20, 7, 4, 19,
                                      16, 2, 18, 8, 13, 12, 3, 5, 11, 6, 14, 15, 2, 13, 16, 1, 20, 4, 10, 19, 14, 17,
                                      18, 3, 9, 8, 12, 7, 6, 5, 11, 4, 11, 6, 19, 12, 1, 14, 17, 2, 5, 20, 3, 13, 8, 9,
                                      18, 16, 10, 15, 7, 17, 6, 13, 14, 5, 20, 1, 7, 19, 16, 15, 18, 3, 11, 10, 8, 4, 9,
                                      2, 12, 17, 8, 15, 9, 7, 20, 5, 1, 19, 16, 10, 18, 6, 13, 12, 2, 3, 11, 4, 14, 4,
                                      11, 6, 19, 12, 2, 14, 17, 1, 5, 20, 3, 13, 8, 9, 18, 16, 10, 15, 7, 2, 9, 3, 19,
                                      10, 18, 12, 15, 17, 1, 20, 5, 11, 6, 7, 16, 14, 8, 13, 4, 17, 10, 15, 2, 9, 20, 7,
                                      4, 19, 16, 1, 18, 8, 13, 12, 3, 5, 11, 6, 14, 2, 9, 4, 19, 11, 13, 14, 17, 10, 3,
                                      20, 1, 12, 6, 7, 18, 16, 8, 15, 5, 17, 5, 13, 15, 2, 20, 3, 8, 19, 14, 16, 18, 1,
                                      11, 10, 9, 6, 7, 4, 12, 10, 6, 5, 17, 8, 20, 11, 15, 19, 7, 18, 13, 9, 1, 2, 16,
                                      14, 4, 12, 3, 11, 6, 5, 17, 7, 20, 10, 15, 19, 8, 18, 14, 9, 2, 1, 16, 13, 3, 12,
                                      4, 17, 9, 15, 8, 7, 20, 5, 2, 19, 16, 10, 18, 6, 13, 12, 1, 3, 11, 4, 14, 17, 8,
                                      15, 10, 7, 20, 5, 2, 19, 16, 12, 18, 6, 13, 11, 3, 1, 9, 4, 14, 12, 5, 7, 17, 6,
                                      20, 9, 14, 19, 11, 18, 16, 8, 3, 2, 15, 13, 1, 10, 4, 17, 6, 14, 13, 5, 20, 2, 7,
                                      19, 16, 15, 18, 3, 11, 10, 8, 4, 9, 1, 12, 7, 8, 4, 17, 9, 20, 12, 15, 19, 6, 18,
                                      11, 10, 2, 3, 16, 14, 5, 13, 1};
    int full_rank_vector_new3[400] = {1, 18, 11, 9, 15, 16, 13, 6, 20, 5, 19, 17, 7, 4, 3, 14, 12, 2, 10, 8, 10, 1, 15,
                                      13, 3, 19, 17, 6, 4, 11, 20, 2, 12, 7, 8, 18, 16, 9, 14, 5, 10, 19, 1, 5, 17, 9,
                                      3, 14, 20, 8, 15, 18, 7, 13, 12, 6, 2, 11, 4, 16, 9, 19, 5, 1, 17, 11, 7, 14, 20,
                                      4, 16, 18, 3, 13, 12, 8, 6, 10, 2, 15, 9, 3, 15, 13, 1, 19, 17, 5, 10, 11, 20, 2,
                                      12, 6, 7, 18, 16, 8, 14, 4, 11, 19, 5, 8, 17, 1, 3, 15, 20, 10, 6, 18, 9, 14, 13,
                                      2, 4, 12, 7, 16, 10, 19, 3, 6, 17, 7, 1, 15, 20, 9, 12, 18, 8, 14, 13, 4, 2, 11,
                                      5, 16, 6, 13, 14, 10, 9, 18, 16, 1, 19, 7, 20, 12, 8, 2, 3, 17, 15, 5, 11, 4, 10,
                                      2, 15, 13, 4, 19, 17, 6, 1, 11, 20, 3, 12, 7, 8, 18, 16, 9, 14, 5, 5, 19, 6, 3,
                                      16, 15, 9, 12, 20, 1, 17, 18, 2, 11, 10, 13, 8, 7, 4, 14, 11, 19, 6, 8, 17, 2, 4,
                                      15, 20, 10, 1, 18, 9, 14, 13, 3, 5, 12, 7, 16, 10, 2, 15, 13, 3, 19, 17, 6, 4, 11,
                                      20, 1, 12, 7, 8, 18, 16, 9, 14, 5, 7, 19, 5, 3, 17, 14, 8, 13, 20, 2, 16, 18, 1,
                                      12, 11, 10, 6, 9, 4, 15, 6, 14, 13, 10, 9, 18, 16, 2, 19, 7, 20, 12, 8, 1, 3, 17,
                                      15, 5, 11, 4, 6, 15, 12, 9, 10, 18, 16, 3, 19, 7, 20, 14, 8, 2, 1, 17, 13, 4, 11,
                                      5, 11, 19, 5, 7, 17, 4, 2, 15, 20, 9, 10, 18, 8, 14, 13, 1, 3, 12, 6, 16, 10, 19,
                                      3, 6, 17, 8, 2, 15, 20, 9, 13, 18, 7, 14, 12, 4, 1, 11, 5, 16, 3, 17, 11, 9, 12,
                                      18, 14, 5, 20, 7, 19, 16, 8, 4, 2, 15, 13, 1, 10, 6, 9, 19, 3, 2, 17, 10, 6, 14,
                                      20, 7, 16, 18, 4, 13, 12, 8, 5, 11, 1, 15, 6, 11, 14, 12, 7, 18, 16, 2, 19, 8, 20,
                                      10, 9, 3, 4, 17, 15, 5, 13, 1};

    s = 0;
    for (int i = 0; i < 20; ++i) {
        for (int j = 0; j < 20; ++j) {
            EXPECT_EQ(full_rank_vector_new1[s++], full_rank[0][i][j]);
        }
    }
    s = 0;
    for (int i = 0; i < 20; ++i) {
        for (int j = 0; j < 20; ++j) {
            EXPECT_EQ(full_rank_vector_new2[s++], full_rank[1][i][j]);
        }
    }
    s = 0;
    for (int i = 0; i < 20; ++i) {
        for (int j = 0; j < 20; ++j) {
            EXPECT_EQ(full_rank_vector_new3[s++], full_rank[2][i][j]);
        }
    }
}

TEST(KBD, full_rank_finder_imbalanced_multiple_group) {
    int s, k = 5, n = 30, size[5] = {6, 5, 5, 6, 8}, cumsum_size[5] = {0, 6, 11, 16, 22}, label[30];
    int bd_stat_number = (((k - 1) * (k)) >> 1);
    double **distance_matrix;
    distance_matrix = alloc_matrix(n, n);
    Euclidean_distance(X1_X2_X3_CONTINUOUS, distance_matrix, n, 1);

    s = 0;
    for (int i = 0; i < k; ++i) {
        for (int j = 0; j < size[i]; ++j) {
            label[s++] = i;
        }
    }
    int *group_relative_location = (int *) malloc(n * sizeof(int));
    find_group_relative_location(group_relative_location, label, cumsum_size, n, k);
    int **index_matrix = alloc_int_matrix(n, n);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            index_matrix[i][j] = j;
        }
    }
    double *distance_matrix_copy;
    distance_matrix_copy = (double *) malloc(n * sizeof(double));
    for (int i = 0; i < n; i++) {
        memcpy(distance_matrix_copy, distance_matrix[i], n * sizeof(double));
        quicksort(distance_matrix_copy, index_matrix[i], 0, n - 1);
    }
    free(distance_matrix_copy);

    int *pairwise_size = (int *) malloc(bd_stat_number * sizeof(int));
    s = 0;
    for (int i = 0; i < (k - 1); ++i) {
        for (int j = (i + 1); j <= (k - 1); ++j) {
            pairwise_size[s] = size[i] + size[j];
            s++;
        }
    }

    int ***full_rank = alloc_int_square_matrix_list(pairwise_size, bd_stat_number);
    full_rank_finder(full_rank, distance_matrix, index_matrix, label, group_relative_location,
                     cumsum_size, size, n, k - 1);

    int true_rank[] = {1, 5, 3, 11, 6, 2, 7, 9, 8, 4, 10, 7, 1, 9, 11, 2, 8, 3, 6, 4, 5, 10, 3, 5, 1, 11, 6, 2, 7, 9, 8,
                       4, 10, 9, 7, 11, 1, 6, 10, 5, 3, 4, 8, 2, 7, 2, 9, 11, 1, 8, 3, 5, 4, 6, 10, 3, 5, 2, 11, 6, 1,
                       7, 9, 8, 4, 10, 9, 5, 11, 8, 3, 10, 1, 4, 2, 6, 7, 9, 5, 11, 7, 4, 10, 3, 1, 2, 8, 6, 9, 5, 11,
                       8, 4, 10, 2, 3, 1, 6, 7, 2, 3, 5, 11, 6, 4, 7, 9, 8, 1, 10, 9, 7, 11, 2, 6, 10, 5, 3, 4, 8, 1, 1,
                       6, 4, 11, 7, 3, 8, 2, 9, 10, 5, 6, 1, 9, 10, 2, 8, 3, 5, 11, 7, 4, 3, 6, 1, 11, 7, 2, 8, 4, 9,
                       10, 5, 8, 5, 10, 1, 4, 9, 3, 7, 11, 2, 6, 7, 3, 9, 10, 1, 8, 2, 6, 11, 5, 4, 3, 6, 2, 11, 7, 1,
                       8, 4, 9, 10, 5, 7, 3, 10, 8, 2, 9, 1, 6, 11, 5, 4, 2, 6, 4, 11, 7, 3, 8, 1, 9, 10, 5, 4, 7, 2,
                       11, 8, 3, 9, 5, 1, 10, 6, 8, 5, 10, 2, 4, 9, 3, 7, 11, 1, 6, 6, 2, 8, 10, 3, 7, 4, 5, 11, 9, 1,
                       1, 5, 3, 12, 6, 2, 4, 11, 9, 7, 10, 8, 9, 1, 11, 12, 2, 10, 3, 8, 6, 4, 7, 5, 3, 5, 1, 12, 6, 2,
                       4, 11, 9, 7, 10, 8, 10, 8, 12, 1, 7, 11, 9, 2, 4, 6, 3, 5, 9, 2, 11, 12, 1, 10, 4, 8, 6, 3, 7, 5,
                       3, 5, 2, 12, 6, 1, 4, 11, 9, 7, 10, 8, 5, 2, 8, 12, 3, 7, 1, 11, 9, 4, 10, 6, 10, 8, 12, 7, 6,
                       11, 9, 1, 3, 5, 2, 4, 10, 7, 12, 8, 6, 11, 9, 4, 1, 5, 3, 2, 10, 7, 12, 9, 4, 11, 8, 6, 3, 1, 5,
                       2, 10, 8, 12, 7, 6, 11, 9, 2, 3, 5, 1, 4, 10, 7, 12, 9, 6, 11, 8, 4, 2, 5, 3, 1, 1, 8, 4, 14, 10,
                       3, 7, 13, 12, 6, 5, 9, 2, 11, 9, 1, 11, 12, 3, 10, 2, 14, 7, 5, 6, 13, 8, 4, 3, 9, 1, 14, 11, 2,
                       8, 10, 13, 7, 6, 5, 4, 12, 10, 5, 12, 1, 4, 11, 6, 14, 2, 7, 8, 13, 9, 3, 9, 3, 11, 12, 1, 10, 4,
                       14, 5, 6, 7, 13, 8, 2, 3, 9, 2, 14, 10, 1, 8, 11, 13, 7, 6, 5, 4, 12, 9, 2, 11, 12, 5, 10, 1, 14,
                       7, 3, 4, 13, 8, 6, 5, 10, 3, 14, 11, 4, 9, 1, 13, 8, 7, 2, 6, 12, 10, 4, 12, 8, 3, 11, 5, 14, 1,
                       6, 7, 13, 9, 2, 8, 4, 11, 13, 5, 10, 3, 14, 9, 1, 2, 12, 6, 7, 6, 5, 10, 13, 7, 9, 3, 14, 11, 2,
                       1, 12, 4, 8, 5, 10, 3, 14, 11, 4, 9, 2, 13, 8, 7, 1, 6, 12, 2, 8, 5, 14, 9, 4, 7, 13, 12, 6, 3,
                       11, 1, 10, 9, 4, 12, 10, 2, 11, 5, 14, 3, 6, 7, 13, 8, 1, 1, 4, 2, 7, 8, 3, 9, 10, 6, 5, 3, 1, 2,
                       8, 6, 4, 9, 10, 5, 7, 2, 3, 1, 7, 8, 4, 9, 10, 5, 6, 5, 7, 6, 1, 9, 4, 3, 10, 8, 2, 5, 3, 4, 8,
                       1, 6, 9, 10, 2, 7, 2, 4, 3, 6, 9, 1, 8, 10, 7, 5, 5, 7, 6, 2, 10, 4, 1, 8, 9, 3, 6, 8, 7, 3, 10,
                       5, 2, 1, 9, 4, 5, 2, 4, 8, 3, 6, 9, 10, 1, 7, 4, 7, 6, 2, 9, 3, 5, 10, 8, 1, 1, 4, 2, 10, 11, 9,
                       8, 6, 3, 7, 5, 8, 1, 5, 11, 10, 9, 7, 3, 4, 6, 2, 3, 4, 1, 10, 11, 9, 8, 6, 2, 7, 5, 3, 6, 4, 1,
                       11, 2, 10, 8, 5, 9, 7, 9, 6, 8, 11, 1, 10, 2, 4, 7, 3, 5, 3, 6, 4, 2, 11, 1, 10, 8, 5, 9, 7, 8,
                       5, 7, 11, 9, 10, 1, 3, 6, 2, 4, 8, 3, 7, 11, 9, 10, 5, 1, 6, 4, 2, 3, 4, 2, 10, 11, 9, 8, 6, 1,
                       7, 5, 8, 5, 7, 11, 9, 10, 2, 3, 6, 1, 4, 8, 3, 7, 11, 9, 10, 5, 2, 6, 4, 1, 1, 5, 3, 9, 11, 6,
                       13, 4, 7, 8, 12, 10, 2, 4, 1, 3, 10, 7, 6, 13, 2, 8, 9, 12, 11, 5, 3, 5, 1, 9, 10, 6, 13, 2, 7,
                       8, 12, 11, 4, 7, 10, 8, 1, 13, 5, 12, 9, 4, 2, 11, 3, 6, 5, 2, 4, 10, 1, 7, 13, 3, 8, 9, 12, 11,
                       6, 6, 10, 7, 5, 11, 1, 13, 8, 2, 3, 12, 9, 4, 9, 12, 10, 4, 13, 7, 1, 11, 6, 5, 2, 3, 8, 4, 3, 2,
                       10, 9, 6, 13, 1, 7, 8, 12, 11, 5, 7, 10, 8, 4, 12, 3, 13, 9, 1, 2, 11, 5, 6, 7, 10, 8, 3, 12, 4,
                       13, 9, 2, 1, 11, 5, 6, 9, 12, 10, 4, 13, 7, 2, 11, 6, 5, 1, 3, 8, 7, 11, 9, 2, 13, 5, 12, 10, 4,
                       3, 8, 1, 6, 2, 5, 3, 9, 11, 6, 13, 4, 7, 8, 12, 10, 1, 1, 10, 11, 9, 6, 4, 8, 5, 2, 7, 3, 4, 1,
                       10, 11, 2, 3, 9, 7, 5, 8, 6, 5, 2, 1, 11, 3, 4, 10, 8, 6, 9, 7, 7, 10, 11, 1, 9, 8, 2, 4, 6, 3,
                       5, 3, 4, 11, 10, 1, 2, 9, 7, 5, 8, 6, 3, 4, 11, 10, 2, 1, 9, 7, 5, 8, 6, 7, 10, 11, 5, 9, 8, 1,
                       3, 6, 2, 4, 7, 10, 11, 6, 9, 8, 4, 1, 5, 3, 2, 3, 10, 11, 7, 9, 8, 6, 4, 1, 5, 2, 7, 10, 11, 5,
                       9, 8, 2, 3, 6, 1, 4, 7, 10, 11, 6, 9, 8, 4, 2, 5, 3, 1, 1, 10, 13, 8, 5, 4, 12, 3, 6, 7, 11, 9,
                       2, 8, 1, 12, 13, 5, 6, 11, 10, 4, 3, 7, 2, 9, 10, 4, 1, 13, 8, 9, 2, 12, 7, 6, 3, 5, 11, 4, 10,
                       13, 1, 6, 5, 12, 2, 7, 8, 11, 9, 3, 6, 8, 13, 10, 1, 4, 12, 9, 2, 3, 11, 5, 7, 5, 9, 13, 10, 2,
                       1, 12, 7, 3, 4, 11, 8, 6, 10, 4, 2, 13, 8, 9, 1, 12, 7, 6, 3, 5, 11, 3, 10, 13, 4, 6, 5, 12, 1,
                       7, 8, 11, 9, 2, 6, 8, 13, 10, 2, 4, 12, 9, 1, 3, 11, 5, 7, 7, 6, 13, 10, 3, 4, 12, 9, 2, 1, 11,
                       5, 8, 10, 4, 3, 13, 8, 9, 2, 12, 7, 6, 1, 5, 11, 7, 2, 13, 12, 5, 6, 11, 10, 4, 3, 9, 1, 8, 2,
                       10, 13, 8, 5, 4, 12, 3, 6, 7, 11, 9, 1, 1, 12, 10, 7, 11, 9, 3, 14, 8, 2, 4, 13, 6, 5, 9, 1, 3,
                       6, 2, 4, 8, 14, 5, 10, 11, 13, 12, 7, 9, 4, 1, 6, 3, 2, 8, 14, 5, 10, 11, 13, 12, 7, 9, 7, 5, 1,
                       6, 4, 8, 14, 2, 10, 11, 13, 12, 3, 9, 2, 3, 6, 1, 4, 8, 14, 5, 10, 11, 13, 12, 7, 9, 4, 2, 6, 3,
                       1, 8, 14, 5, 10, 11, 13, 12, 7, 2, 12, 10, 6, 11, 9, 1, 14, 7, 3, 4, 13, 8, 5, 6, 14, 12, 9, 13,
                       11, 7, 1, 10, 5, 4, 2, 3, 8, 9, 7, 4, 2, 6, 3, 8, 14, 1, 10, 11, 13, 12, 5, 2, 12, 10, 7, 11, 9,
                       4, 14, 8, 1, 3, 13, 5, 6, 3, 12, 10, 7, 11, 9, 4, 14, 8, 2, 1, 13, 5, 6, 6, 14, 12, 9, 13, 11, 7,
                       2, 10, 5, 4, 1, 3, 8, 4, 13, 11, 8, 12, 10, 5, 14, 9, 3, 2, 7, 1, 6, 7, 10, 6, 2, 9, 5, 4, 14, 3,
                       8, 11, 13, 12, 1};
    s = 0;
    for (int l = 0; l < bd_stat_number; ++l) {
        for (int i = 0; i < pairwise_size[l]; ++i) {
            for (int j = 0; j < pairwise_size[l]; ++j) {
                EXPECT_EQ(true_rank[s++], full_rank[l][i][j]);
            }
        }
    }
}

TEST(KBD, permute_bd_value) {
    int s, k = 3, n = 30, size[3] = {10, 10, 10}, cumsum_size[3] = {0, 10, 20};
    int bd_stat_number = (((k - 1) * (k)) >> 1);
    double *x123_continuous, **distance_matrix;
    x123_continuous = (double *) malloc(n * sizeof(double));
    distance_matrix = alloc_matrix(n, n);
    for (int i = 0; i < 10; ++i) {
        x123_continuous[i] = X1_CONTINUOUS[i];
        x123_continuous[i + 10] = X2_CONTINUOUS[i];
        x123_continuous[i + 20] = X3_CONTINUOUS[i];
    }
    Euclidean_distance(x123_continuous, distance_matrix, n, 1);

    int label[30] = {0, 1, 1, 2, 0, 2, 2, 1, 1, 0, 0, 0, 1, 0, 1, 2, 2, 2, 1, 0, 0, 2, 2, 0, 2, 1, 2, 1, 0, 1};
    int *group_relative_location = (int *) malloc(n * sizeof(int));
    find_group_relative_location(group_relative_location, label, cumsum_size, n, k);

    int **index_matrix = alloc_int_matrix(n, n);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            index_matrix[i][j] = j;
        }
    }
    double *distance_matrix_copy;
    distance_matrix_copy = (double *) malloc(n * sizeof(double));
    for (int i = 0; i < n; i++) {
        memcpy(distance_matrix_copy, distance_matrix[i], n * sizeof(double));
        quicksort(distance_matrix_copy, index_matrix[i], 0, n - 1);
    }
    free(distance_matrix_copy);

    int *pairwise_size = (int *) malloc(bd_stat_number * sizeof(int));
    s = 0;
    for (int i = 0; i < (k - 1); ++i) {
        for (int j = (i + 1); j <= (k - 1); ++j) {
            pairwise_size[s] = size[i] + size[j];
            s++;
        }
    }

    int ***sub_rank = alloc_int_square_matrix_list(size, k);
    sub_rank_finder(sub_rank, distance_matrix, index_matrix, label, group_relative_location, cumsum_size, n, k - 1);
    int ***full_rank = alloc_int_square_matrix_list(pairwise_size, bd_stat_number);
    full_rank_finder(full_rank, distance_matrix, index_matrix, label, group_relative_location,
                     cumsum_size, size, n, k - 1);
    double **bd_stat_array = alloc_matrix(bd_stat_number, 2);
    ball_divergence_array(bd_stat_array, full_rank, sub_rank, size, k);
    double kbd_stat_tmp[6];
    k_ball_divergence_from_by_sample_ball_divergence(kbd_stat_tmp, bd_stat_array, bd_stat_number, k);
    EXPECT_NEAR(kbd_stat_tmp[0], 0.1284, ABSOLUTE_ERROR);
}

TEST(KBD, sort_ints) {
    int array[10] = {0, 1, 2, 0, 1, 2, 0, 1, 2, 0};
    sort_ints(array, 10);
    int sorted_array[10] = {0, 0, 0, 0, 1, 1, 1, 2, 2, 2};
    for (int i = 0; i < 10; ++i) {
        EXPECT_EQ(sorted_array[i], array[i]);
    }
}