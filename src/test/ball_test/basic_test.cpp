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

TEST(utilities, batch_pvalue) {
    double stat[100], p_value[100], true_p_value[100];
    double permuted_stat[1000];
    for (int i = 0; i < 100; ++i) {
        stat[i] = i / 100.0;
    }
    for (int i = 0; i < 1000; ++i) {
        permuted_stat[i] = i / 1000.0;
    }
    compute_batch_pvalue(stat, permuted_stat, p_value, 100, 1000);
    for (int i = 0; i < 100; ++i) {
        true_p_value[i] = compute_pvalue(stat[i], permuted_stat, 1000);
    }
    for (int i = 0; i < 100; ++i) {
        EXPECT_NEAR(true_p_value[i], p_value[i], ABSOLUTE_ERROR);
    }
}

TEST(utilities, quick_rank_min) {
    double x[6] = {1.2, 1.3, 1.3, 1.3, 1.3, 0.9};
    int r[6] = {0, 0, 0, 0, 0, 0}, true_r[6] = {2, 3, 3, 3, 3, 1};
    quick_rank_min(x, r, 6);
    for (int i = 0; i < 6; ++i) {
        EXPECT_EQ(r[i], true_r[i]);
    }
}

TEST(utilities, count_of_larger_number_before_self) {
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

TEST(utilities, beautify_time) {
    char result1[200] = "";
    int second = 100000;
    beautify_time(result1, second);
    char true_result1[200] = "1 day, 3 hours, 46 minutes, 40 seconds";
    for (int i = 0; i < strlen(result1); ++i) {
        EXPECT_EQ(true_result1[i], result1[i]);
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
    sub_rank_finder(sub_rank, distance_matrix, index_matrix, label, group_relative_location, cumsum_size, size, n,
                    k - 1);
    int s = 0;
    for (int l = 0; l < 3; ++l) {
        for (int i = 0; i < 10; ++i) {
            for (int j = 0; j < 10; ++j) {
                EXPECT_EQ(SUB_RANK_BALANCE_CONTINUOUS1[s++], sub_rank[l][i][j]);
            }
        }
    }
    sub_rank_finder_tie(sub_rank, distance_matrix, index_matrix, label, group_relative_location, cumsum_size, size, n,
                        k - 1);
    s = 0;
    for (int l = 0; l < 3; ++l) {
        for (int i = 0; i < 10; ++i) {
            for (int j = 0; j < 10; ++j) {
                EXPECT_EQ(SUB_RANK_BALANCE_CONTINUOUS1[s++], sub_rank[l][i][j]);
            }
        }
    }

    int label_new[30] = {0, 1, 1, 2, 0, 2, 2, 1, 1, 0, 0, 0, 1, 0, 1, 2, 2, 2, 1, 0, 0, 2, 2, 0, 2, 1, 2, 1, 0, 1};
    find_group_relative_location(group_relative_location, label_new, cumsum_size, n, k);
    sub_rank_finder(sub_rank, distance_matrix, index_matrix, label_new, group_relative_location, cumsum_size, size,
                    n, k - 1);
    s = 0;
    for (int l = 0; l < 3; ++l) {
        for (int i = 0; i < 10; ++i) {
            for (int j = 0; j < 10; ++j) {
                EXPECT_EQ(SUB_RANK_BALANCE_CONTINUOUS2[s++], sub_rank[l][i][j]);
            }
        }
    }
    sub_rank_finder_tie(sub_rank, distance_matrix, index_matrix, label_new, group_relative_location, cumsum_size, size,
                        n, k - 1);
    s = 0;
    for (int l = 0; l < 3; ++l) {
        for (int i = 0; i < 10; ++i) {
            for (int j = 0; j < 10; ++j) {
                EXPECT_EQ(SUB_RANK_BALANCE_CONTINUOUS2[s++], sub_rank[l][i][j]);
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
    sub_rank_finder(sub_rank, distance_matrix, index_matrix, label, group_relative_location, cumsum_size, size,
                    n, k - 1);

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
    s = 0;
    for (int l = 0; l < 3; ++l) {
        for (int i = 0; i < 20; ++i) {
            for (int j = 0; j < 20; ++j) {
                EXPECT_EQ(FULL_RANK_BALANCE_CONTINUOUS1[s++], full_rank[l][i][j]);
            }
        }
    }
    full_rank_finder_tie(full_rank, distance_matrix, index_matrix, label, group_relative_location,
                         cumsum_size, size, n, k - 1);
    s = 0;
    for (int l = 0; l < 3; ++l) {
        for (int i = 0; i < 20; ++i) {
            for (int j = 0; j < 20; ++j) {
                EXPECT_EQ(FULL_RANK_BALANCE_CONTINUOUS1[s++], full_rank[l][i][j]);
            }
        }
    }

    int label_new[30] = {0, 1, 1, 2, 0, 2, 2, 1, 1, 0, 0, 0, 1, 0, 1, 2, 2, 2, 1, 0, 0, 2, 2, 0, 2, 1, 2, 1, 0, 1};
    find_group_relative_location(group_relative_location, label_new, cumsum_size, n, k);
    full_rank_finder(full_rank, distance_matrix, index_matrix, label_new,
                     group_relative_location, cumsum_size, size, n, k - 1);
    s = 0;
    for (int l = 0; l < 3; ++l) {
        for (int i = 0; i < 20; ++i) {
            for (int j = 0; j < 20; ++j) {
                EXPECT_EQ(FULL_RANK_BALANCE_CONTINUOUS2[s++], full_rank[l][i][j]);
            }
        }
    }
    full_rank_finder_tie(full_rank, distance_matrix, index_matrix, label_new,
                         group_relative_location, cumsum_size, size, n, k - 1);
    s = 0;
    for (int l = 0; l < 3; ++l) {
        for (int i = 0; i < 20; ++i) {
            for (int j = 0; j < 20; ++j) {
                EXPECT_EQ(FULL_RANK_BALANCE_CONTINUOUS2[s++], full_rank[l][i][j]);
            }
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
    s = 0;
    for (int l = 0; l < bd_stat_number; ++l) {
        for (int i = 0; i < pairwise_size[l]; ++i) {
            for (int j = 0; j < pairwise_size[l]; ++j) {
                EXPECT_EQ(FULL_RANK_IMBALANCE_CONTINUOUS1[s++], full_rank[l][i][j]);
            }
        }
    }
}

TEST(KBD, sub_rank_finder_tie) {
    int k = 3, n = 30, size[3] = {10, 10, 10}, cumsum_size[3] = {0, 10, 20}, label[30];
    double *x123_discrete, **distance_matrix;
    x123_discrete = (double *) malloc(n * sizeof(double));
    distance_matrix = alloc_matrix(n, n);
    for (int i = 0; i < 10; ++i) {
        x123_discrete[i] = X1_DISCRETE[i];
        x123_discrete[i + 10] = X2_DISCRETE[i];
        x123_discrete[i + 20] = X3_DISCRETE[i];
        label[i] = 0;
        label[i + 10] = 1;
        label[i + 20] = 2;
    }
    Euclidean_distance(x123_discrete, distance_matrix, n, 1);

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
    sub_rank_finder_tie(sub_rank, distance_matrix, index_matrix, label, group_relative_location, cumsum_size, size, n,
                        k - 1);

    int true_rank[300] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 3, 1, 3, 4, 5, 6, 7, 8, 9, 10, 5, 3, 1, 3, 5, 6, 7, 8, 9, 10,
                          7, 5, 3, 1, 3, 5, 7, 8, 9, 10, 9, 7, 5, 3, 1, 3, 5, 7, 9, 10, 10, 9, 7, 5, 3, 1, 3, 5, 7, 9,
                          10, 9, 8, 7, 5, 3, 1, 3, 5, 7, 10, 9, 8, 7, 6, 5, 3, 1, 3, 5, 10, 9, 8, 7, 6, 5, 4, 3, 1, 3,
                          10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 3, 1, 3, 4, 5, 6, 7, 8, 9, 10,
                          5, 3, 1, 3, 5, 6, 7, 8, 9, 10, 7, 5, 3, 1, 3, 5, 7, 8, 9, 10, 9, 7, 5, 3, 1, 3, 5, 7, 9, 10,
                          10, 9, 7, 5, 3, 1, 3, 5, 7, 9, 10, 9, 8, 7, 5, 3, 1, 3, 5, 7, 10, 9, 8, 7, 6, 5, 3, 1, 3, 5,
                          10, 9, 8, 7, 6, 5, 4, 3, 1, 3, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10,
                          3, 1, 3, 4, 5, 6, 7, 8, 9, 10, 5, 3, 1, 3, 5, 6, 7, 8, 9, 10, 7, 5, 3, 1, 3, 5, 7, 8, 9, 10,
                          9, 7, 5, 3, 1, 3, 5, 7, 9, 10, 10, 9, 7, 5, 3, 1, 3, 5, 7, 9, 10, 9, 8, 7, 5, 3, 1, 3, 5, 7,
                          10, 9, 8, 7, 6, 5, 3, 1, 3, 5, 10, 9, 8, 7, 6, 5, 4, 3, 1, 3, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1};
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
    sub_rank_finder_tie(sub_rank, distance_matrix, index_matrix, label_new, group_relative_location, cumsum_size, size,
                        n, k - 1);

    int true_rank_new[300] = {1, 3, 6, 7, 8, 9, 10, 6, 4, 2, 4, 1, 6, 7, 8, 9, 10, 6, 2, 3, 9, 7, 2, 3, 4, 6, 10, 2, 5,
                              8, 10, 7, 4, 1, 4, 5, 9, 4, 6, 9, 10, 7, 5, 2, 1, 5, 8, 5, 6, 9, 10, 8, 5, 3, 2, 1, 6, 5,
                              7, 9, 10, 8, 6, 4, 3, 2, 1, 6, 7, 9, 9, 7, 2, 3, 4, 6, 10, 2, 5, 8, 8, 2, 4, 5, 7, 9, 10,
                              4, 1, 7, 2, 3, 6, 7, 8, 9, 10, 6, 4, 1, 1, 4, 6, 7, 8, 9, 10, 5, 4, 4, 3, 2, 6, 7, 8, 9,
                              10, 5, 2, 5, 7, 6, 1, 2, 6, 9, 10, 3, 6, 9, 8, 7, 2, 1, 4, 7, 10, 4, 7, 9, 9, 8, 4, 3, 1,
                              2, 5, 6, 8, 10, 9, 8, 5, 4, 2, 1, 3, 6, 8, 10, 9, 8, 5, 4, 3, 2, 1, 6, 8, 10, 5, 3, 5, 7,
                              8, 9, 10, 1, 3, 7, 3, 2, 6, 7, 8, 9, 10, 5, 2, 5, 2, 4, 6, 7, 8, 9, 10, 5, 4, 1, 2, 4, 5,
                              8, 9, 10, 7, 6, 4, 2, 6, 2, 3, 8, 9, 10, 7, 6, 2, 6, 7, 4, 1, 8, 9, 10, 5, 4, 4, 7, 10, 8,
                              6, 1, 2, 3, 4, 5, 8, 10, 10, 8, 6, 3, 1, 3, 4, 5, 8, 10, 10, 8, 6, 3, 2, 1, 4, 5, 8, 10,
                              7, 5, 3, 8, 9, 10, 1, 2, 5, 7, 7, 5, 3, 8, 9, 10, 3, 1, 5, 7, 6, 2, 3, 8, 9, 10, 7, 6, 2,
                              6, 2, 4, 5, 8, 9, 10, 7, 6, 4, 2};
    s = 0;
    for (int l = 0; l < 3; ++l) {
        for (int i = 0; i < 10; ++i) {
            for (int j = 0; j < 10; ++j) {
                EXPECT_EQ(true_rank_new[s++], sub_rank[l][i][j]);
            }
        }
    }

    free(x123_discrete);
    free_int_square_matrix_list(sub_rank, size, k);
    free_int_matrix(index_matrix, n, n);
    free_matrix(distance_matrix, n, n);
    free(group_relative_location);
}

TEST(KBD, sub_rank_finder_imbalanced_multiple_group_tie) {
    int s, k = 5, n = 30, size[5] = {6, 5, 5, 6, 8}, cumsum_size[5] = {0, 6, 11, 16, 22}, label[30];
    double **distance_matrix;
    distance_matrix = alloc_matrix(n, n);
    double X1_X2_X3_DISCRETE[30] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20,
                                    10, 9, 8, 7, 6, 5, 4, 3, 2, 1};
    Euclidean_distance(X1_X2_X3_DISCRETE, distance_matrix, n, 1);

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
    sub_rank_finder_tie(sub_rank, distance_matrix, index_matrix, label, group_relative_location, cumsum_size, size,
                        n, k - 1);

    int true_rank[] = {1, 2, 3, 4, 5, 6, 3, 1, 3, 4, 5, 6, 5, 3, 1, 3, 5, 6, 6, 5, 3, 1, 3, 5, 6, 5, 4, 3, 1, 3, 6, 5,
                       4, 3, 2, 1, 1, 2, 3, 4, 5, 3, 1, 3, 4, 5, 5, 3, 1, 3, 5, 5, 4, 3, 1, 3, 5, 4, 3, 2, 1, 1, 2, 3,
                       4, 5, 3, 1, 3, 4, 5, 5, 3, 1, 3, 5, 5, 4, 3, 1, 3, 5, 4, 3, 2, 1, 1, 2, 3, 4, 5, 6, 3, 1, 3, 4,
                       5, 6, 4, 3, 1, 3, 5, 6, 4, 3, 2, 1, 5, 6, 3, 4, 5, 6, 1, 2, 3, 4, 5, 6, 2, 1, 1, 2, 3, 4, 5, 6,
                       7, 8, 3, 1, 3, 4, 5, 6, 7, 8, 5, 3, 1, 3, 5, 6, 7, 8, 7, 5, 3, 1, 3, 5, 7, 8, 8, 7, 5, 3, 1, 3,
                       5, 7, 8, 7, 6, 5, 3, 1, 3, 5, 8, 7, 6, 5, 4, 3, 1, 3, 8, 7, 6, 5, 4, 3, 2, 1};
    s = 0;
    for (int l = 0; l < 5; ++l) {
        for (int i = 0; i < size[l]; ++i) {
            for (int j = 0; j < size[l]; ++j) {
                EXPECT_EQ(true_rank[s++], sub_rank[l][i][j]);
            }
        }
    }
}

TEST(KBD, full_rank_finder_tie) {
    int s, k = 3, n = 30, size[3] = {10, 10, 10}, cumsum_size[3] = {0, 10, 20}, label[30];
    int bd_stat_number = (((k - 1) * (k)) >> 1);
    double **distance_matrix;
    distance_matrix = alloc_matrix(n, n);
    double X1_X2_X3_DISCRETE[30] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20,
                                    10, 9, 8, 7, 6, 5, 4, 3, 2, 1};
    Euclidean_distance(X1_X2_X3_DISCRETE, distance_matrix, n, 1);

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
    full_rank_finder_tie(full_rank, distance_matrix, index_matrix, label, group_relative_location,
                         cumsum_size, size, n, k - 1);
    s = 0;
    for (int l = 0; l < 3; ++l) {
        for (int i = 0; i < 20; ++i) {
            for (int j = 0; j < 20; ++j) {
                EXPECT_EQ(FULL_RANK_BALANCE_DISCRETE1[s++], full_rank[l][i][j]);
            }
        }
    }


    int label_new[30] = {0, 1, 1, 2, 0, 2, 2, 1, 1, 0, 0, 0, 1, 0, 1, 2, 2, 2, 1, 0, 0, 2, 2, 0, 2, 1, 2, 1, 0, 1};
    find_group_relative_location(group_relative_location, label_new, cumsum_size, n, k);
    full_rank_finder_tie(full_rank, distance_matrix, index_matrix, label_new,
                         group_relative_location, cumsum_size, size, n, k - 1);

    s = 0;
    for (int l = 0; l < 3; ++l) {
        for (int i = 0; i < 20; ++i) {
            for (int j = 0; j < 20; ++j) {
                EXPECT_EQ(FULL_RANK_BALANCE_DISCRETE2[s++], full_rank[l][i][j]);
            }
        }
    }
}

TEST(KBD, full_rank_finder_imbalanced_multiple_group_tie) {
    int s, k = 5, n = 30, size[5] = {6, 5, 5, 6, 8}, cumsum_size[5] = {0, 6, 11, 16, 22}, label[30];
    int bd_stat_number = (((k - 1) * (k)) >> 1);
    double **distance_matrix;
    distance_matrix = alloc_matrix(n, n);
    double X1_X2_X3_DISCRETE[30] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20,
                                    10, 9, 8, 7, 6, 5, 4, 3, 2, 1};
    Euclidean_distance(X1_X2_X3_DISCRETE, distance_matrix, n, 1);

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
    full_rank_finder_tie(full_rank, distance_matrix, index_matrix, label, group_relative_location,
                         cumsum_size, size, n, k - 1);
    s = 0;
    for (int l = 0; l < bd_stat_number; ++l) {
        for (int i = 0; i < pairwise_size[l]; ++i) {
            for (int j = 0; j < pairwise_size[l]; ++j) {
                EXPECT_EQ(FULL_RANK_IMBALANCE_DISCRETE1[s++], full_rank[l][i][j]);
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
    sub_rank_finder(sub_rank, distance_matrix, index_matrix, label, group_relative_location, cumsum_size, size,
                    n, k - 1);
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

TEST(KBD, bd_gwas_refining) {
    int R = 10;
    double bd_stat[4], pvalue[4], permuted_bd_stat[40];
    int refine_size[5] = {6, 7, 7, 10, 10};
    int refine_k_num[2] = {3, 2};
    int nth = 1, n = 20, refine_num = 2, verbose = 0;
    bd_gwas_refining(bd_stat, permuted_bd_stat, pvalue, X1_X2_CONTINUOUS_DST, &n,
                     &refine_num, refine_size, refine_k_num, &R, &nth, &verbose);
    for (int i = 0; i < refine_num; ++i) {
        EXPECT_EQ(pvalue[2 * i], 1.0);
        EXPECT_EQ(pvalue[2 * i + 1], 1.0);
    }

    nth = 2;
    bd_gwas_refining(bd_stat, permuted_bd_stat, pvalue, X1_X2_CONTINUOUS_DST, &n,
                     &refine_num, refine_size, refine_k_num, &R, &nth, &verbose);
    for (int i = 0; i < refine_num; ++i) {
        EXPECT_EQ(pvalue[2 * i], 1.0);
        EXPECT_EQ(pvalue[2 * i + 1], 1.0);
    }
}