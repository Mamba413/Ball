//
// Created by JinZhu on 2019/3/18.
//

#ifndef BALL_KBD_H
#define BALL_KBD_H

void sub_rank_finder(int ***sub_rank, double **distance_matrix, int **index_matrix, int *label,
                     int *group_relative_location, int *cumsum_size, int num, int max_k);

void full_rank_finder(int ***full_rank, double **distance_matrix, int **index_matrix, int *label,
                      int *group_relative_location, int *cumsum_size, int *size, int num, int k_max);

void ball_divergence_array(double **bd_stat_array, int ***full_rank, int ***sub_rank, int *size, int K);

void k_ball_divergence_from_by_sample_ball_divergence(double *kbd_stat, double **bd_stat_array,
                                                      int bd_stat_number, int k);

void KBD3(double *kbd_stat, double *pvalue, double *xy, int *size, int *n, int *k, int *R, int *thread);

#endif //BALL_KBD_H
