//
// Created by JinZhu on 2019/10/13.
//

#ifndef BALL_BDD_MATRIX_H
#define BALL_BDD_MATRIX_H

void bdd_matrix_bias_two_group(double *b_dd, double *x, int *n1_num, int *n2_num, int *nthread);
void bdd_matrix_bias(double *b_dd, double *x, int *n, int *nthread, int *weight_type);

#endif //BALL_BDD_MATRIX_H
