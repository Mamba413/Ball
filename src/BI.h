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


#define MAX(a,b) (((a)>(b))?(a):(b));
void Merge(int *permutation, int *source, int *inversion_count, int dim, int n);
int Inversions(int *permutation, int *source, int *inversion_count,int dim, int n);
double Ball_Information(int *n, double **Dx, double **Dy, int **xidx, int **yidx, int *i_perm, int *i_perm_inv, int *weight);
void BI(double *bcor, double *permuted_bcor, double *x, double *y, int *n, int *R, int *weight);
void computeRank(int n, int **Rank);
void initRank(int n, int **Rank, int *xrank, int *yrank, int *i_perm);
void ranksort(int *n, int *zrank, double *z, int *zidx);
void sort(int *n, int *zidx, double *z, int **dzidx);
void createidx(int *n, int *zidx, double *z, int **lowzidx, int **higzidx);
double U_Ball_Information(int *n, int **Rank, int **lowxidx, int **higxidx, int **lowyidx, int **higyidx, int *i_perm, int *weight);
void UBI(double *bcor, double *permuted_bcor, double *x, double *y, int *n, int *R, int *weight);
double ubcov_value(double *x, double *y, int *n, int *weight);
double bcov_value(double *x, double *y, int *n, int *weight);
// R API function:
void bcov_stat(double *bcor, double *x, double *y, int *n, int *weight, int *dst, int *type);
void bcov_test(double *bcor, double *permuted_bcor, double *x, double *y, int *n, int *R, int *weight, int *dst, int *type);


#endif /* BI_H_ */