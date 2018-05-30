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


#ifndef BCOR_H_
#define BCOR_H_

//double ubcov_value(double *x, double *y, int *n, int *weight, int *thread);
//double bcov_value(double *x, double *y, int *n, int *weight, int *thread);
//double bcor_value(double *x, double *y, int *n, int *weight, int *dst, int *thread);
//void bcov_stat(double *bcor, double *x, double *y, int *n, int *weight, int *dst, int *type, int *thread);

void _bcor_test(double *bcorsis_stat, double *y, double *x, int *x_number, int *f_number, int *n, int *p, int *nthread);
void _fast_bcor_test(double *bcorsis_stat, double *y, double *x, int *f_number, int *n, int *nthread);

// R API function:
void bcor_test(double *bcorsis_stat, double *y, double *x, int *x_number, int *f_number, int *size, int *n, int *p, int *k, int *dst, int *nthread);

#endif /* BCOR_H_ */
