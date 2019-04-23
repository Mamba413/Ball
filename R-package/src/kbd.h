//
// Created by JinZhu on 2019/3/18.
//

#ifndef BALL_KBD_H
#define BALL_KBD_H

#include "stdlib.h"

void sub_rank_finder(int ***sub_rank, double **distance_matrix, int **index_matrix, int *label,
                     int *group_relative_location, int *cumsum_size, int *size, int num, int max_k);

void sub_rank_finder_tie(int ***sub_rank, double **distance_matrix, int **index_matrix, const int *label,
                         int *group_relative_location, int *cumsum_size, int *size, int num, int k_max);

void full_rank_finder(int ***full_rank, double **distance_matrix, int **index_matrix, int *label,
                      int *group_relative_location, int *cumsum_size, int *size, int num, int k_max);

void full_rank_finder_tie(int ***full_rank, double **distance_matrix, int **index_matrix, int *label,
                          int *group_relative_location, int *cumsum_size, int *size, int num, int k_max);

void ball_divergence_array(double **bd_stat_array, int ***full_rank, int ***sub_rank, int *size, int K);

void k_ball_divergence_from_by_sample_ball_divergence(double *kbd_stat, double **bd_stat_array,
                                                      int bd_stat_number, int k);

void sort_ints(int *a, size_t n);

void KBD3(double *, double *, double *, int *, int *, int *, int *, int *);

void bd_gwas_test(double *bd_stat, double *permuted_bd_stat, double *pvalue, double *xy, const int *snp,
                  const int *n, const int *p, int *R, int *nthread);

#endif //BALL_KBD_H
