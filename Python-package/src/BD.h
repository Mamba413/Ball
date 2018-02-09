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
#ifndef BD_H_
#define BD_H_


void ranksort2(int n, int **Rxy, double **Dxy, int **Ixy);
void Findx2(int *Rxy, int *Ixy, int *i_perm, int *n1, int *n2, int *Rx);
void Findx(int **Rxy, int **Ixy, int *i_perm, int *n1, int *n2, int **Rx);
double Ball_Divergence(int **Rxy, int **Rx, int *i_perm_tmp, int *n1, int *n2, int *weight);
void resample3(int *i_perm, int *i_perm_tmp, int n, int *n1);
void ranksort3(int n, int *xyidx, double *xy, int **Rxy, int **Ixy);
void get_ij_dst(double *xy, double *ij_dst, int *cumulate_size, int *size, int *n, int *p, int *q);
void compute_cumulate_size(int *cumulate_size, int *size, int *k);
void permute_dst(double *xy, double *new_xy, int *index, int *N);
void get_ij_value(double *xy, double *ij_value, int *cumulate_size, int *size, int *p, int *q);
extern double bd_value(double *xy, int *n1, int *n2, int *weight);
extern double ubd_value(double *xy, int *n1, int *n2, int *weight);
double kbd_value(double *xy, int *size, int *n, int *k, int *weight);
double ukbd_value(double *xy, int *size, int *k, int *weight);
extern void BD(double *bd, double *permuted_bd, double *xy, int *n1, int *n2, int *p, int *dst, int *R, int *weight);
extern void UBD(double *bd, double *permuted_bd, double *xy, int *n1, int *n2, int *R, int *weight);
void KBD(double *kbd, double *permuted_kbd, double *xy, int *size, int *n, int *k, int *R, int *weight);
void UKBD(double *kbd, double *permuted_kbd, double *xy, int *size, int *n, int *k, int *R, int *weight);
// R API function:
//void bd_stat(double *bd, double *xy, int *n1, int *n2, int *weight, int *dst);
void kbd_stat(double *bd, double *xy, int *size, int *n, int *k, int *weight, int *dst);
//void bd_test(double *bd, double *permuted_bd, double *xy, int *n1, int *n2, int *p, int *dst, int *R, int *weight);
void kbd_test(double *kbd, double *permuted_kbd, double *xy, int *size, int *n, int *k, int *dst, int *R, int *weight);

#endif /* BD_H_ */