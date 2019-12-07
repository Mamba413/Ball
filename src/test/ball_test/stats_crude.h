//
// Created by JinZhu on 2019/3/15.
//

#ifndef BALL_STATS_CRUDE_H
#define BALL_STATS_CRUDE_H

void Ball_Divergence_Crude(double *bd_stat, double **Dx, int n1, int n2);
void Ball_Covariance_Crude(double *bcov_stat, double **Dx, double **Dy, int n);
void Ball_Correlation_Crude(double *bcor_stat, double **Dx, double **Dy, int n);
void K_Ball_Covariance_Crude(double *kbcov_stat, double ***Dx, int n, int k);

#endif //BALL_STATS_CRUDE_H
