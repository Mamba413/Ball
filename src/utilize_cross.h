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
#ifndef UTILIZE_CROSS_H_
#define UTILIZE_CROSS_H_

void print_stop_message_internal();
void R_check_interrupt_fn(void *dummy);
int pending_interrupt_status();
int r_available_rand();
int random_index(int n, int i);
int random_index2(int i);
int random_index_thread(int i);


#endif /* UTILIZE_CROSS_H_ */