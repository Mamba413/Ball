/*!
*    Copyright 2017 by ChengFeng Liu, Jin Zhu<zhuj37mail2.sysu.edu.cn>
*     This program is free software: you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation, either version 3 of the License, or
* (at your option) any later version.

* This program is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
* GNU General Public License for more details.

* You should have received a copy of the GNU General Public License
* along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/
#ifndef BI_H_
#define BI_H_


void Ball_Information(double *bcov_stat, int *n, double **Dx, double **Dy, int **xidx, int **yidx, int *i_perm, int *i_perm_inv);
void Ball_Information_parallel(double *bcov_stat, int *n, double **Dx, double **Dy, int **xidx, int **yidx, int *i_perm, int *i_perm_inv, int *nthread);
void Ball_Information_wrapper(double *bcov_stat, int *n, double **Dx, double **Dy, int **xidx, int **yidx, int *i_perm, int *i_perm_inv, int *nthread);
void BI(double *bcov, double *permuted_bcov, double *x, double *y, int *n, int *R, int *thread);
void ranksort(int *n, int *zrank, double *z, int *zidx);
void U_Ball_Information(double *, int *n, int **Rank, int **lowxidx, int **higxidx, int **lowyidx, int **higyidx, int *i_perm);
void U_Ball_Information_parallel(double *, int *n, int **Rank, int **lowxidx, int **higxidx, int **lowyidx, int **higyidx, int *i_perm, int *nthread);
void U_Ball_Information_wrapper(double *, int *n, int **Rank, int **lowxidx, int **higxidx, int **lowyidx, int **higyidx, int *i_perm, int *nthread);
void UBI(double *bcor, double *permuted_bcor, double *x, double *y, int *n, int *R, int *thread);
double ubcov_value(double *x, double *y, int *n, int *weight, int *thread);
double bcov_value(double *x, double *y, int *n, int *weight, int *thread);
double bcor_value(double *x, double *y, int *n, int *weight, int *dst, int *thread);
// R API function:
void bcov_stat(double *bcor, double *x, double *y, int *n, int *weight, int *dst, int *type, int *thread);
void bcov_test(double *bcor, double *permuted_bcor, double *x, double *y, int *n, int *R, int *weight, int *dst, int *type, int *thread);


#endif /* BI_H_ */