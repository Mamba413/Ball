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

#include "R.h"
#include "Rinternals.h"


void print_stop_message_internal()
{
  Rprintf("Process stop due to user interruption! \n");
  return;
}


/* R_check_interrupt_fn and pending_interrupt_status
 * are used to check for interrupt without long jumping 
 * 
 * R_check_interrupt_fn is a R internal function;
 */
void R_check_interrupt_fn(void *dummy) {
  R_CheckUserInterrupt();
}


/*
 * API function for utilities.c
 */
int pending_interrupt_status()
{
  return !(R_ToplevelExec(R_check_interrupt_fn, NULL));
}


/*
 * A random integer generator available for R-package building.
 * Following the rule in "Writing R Extension"[https://cran.r-project.org/doc/manuals/r-release/R-exts.html#Random-numbers]!
 * output:
 * A random integer
 */
int r_available_rand()
{
  GetRNGstate();
  int random_value = round(RAND_MAX*unif_rand());
  PutRNGstate();
  return random_value;
}


/*
 * API function for utilities.c
 */
int random_index(int n, int i)
{
  int index = i + r_available_rand() / (RAND_MAX / (n - i) + 1);
  return index;
}


/*
 * API function for utilities.c
 */
int random_index2(int i)
{
  int index = r_available_rand() % (i + 1);
  return index;
}