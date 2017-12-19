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
#ifndef UTILITIES_H_
#define UTILITIES_H_


void quicksort(double *a, int *idx, int l, int u);
void quicksort2(double *a, double *b, int *idx, int l, int u);
double **alloc_matrix(int r, int c);
int **alloc_int_matrix(int r, int c);
void free_matrix(double **matrix, int r, int c);
void free_int_matrix(int **matrix, int r, int c);
void vector2matrix(double *x, double **y, int N, int d, int isroworder);
void Euclidean_distance(double *x, double **Dx, int n, int d);
void distance(double *x, double *Dx, int *n, int *d);
void shuffle(int *array, int *N);
void shuffle_value(double *array, int *N);
void check_interrupt_fn(void *dummy);
int pending_interrupt();
int r_available_rand();


#endif /* UTILITIES_H_ */