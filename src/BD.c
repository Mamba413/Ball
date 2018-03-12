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
#ifdef _OPENMP
# include "omp.h"
#endif

#include "math.h"
#include "stdlib.h" 
#include "stdio.h"
#include "BD.h"
#include "utilities.h"


void ranksort2(int n, int **Rxy, double **Dxy, int **Ixy)
{
  int i, j, lastpos = n - 1;
  double lastval;
  for(i = 0; i < n; i++){
    lastval = -1;
    for(j = n - 1; j >= 0; j--){
      if(lastval != Dxy[i][j] )
        lastpos = j;
      lastval =  Dxy[i][j];
      Rxy[i][Ixy[i][j]] = lastpos;
    }
  }
}


void Findx2(int *Rxy, int *Ixy, int *i_perm, int *n1, int *n2, int *Rx)
{
  int j, lastpos, lastval, n, tmp;
  n = *n1 + *n2;
  
  lastpos = *n1 - 1;
  if(i_perm[Ixy[n-1]] == 1){
    tmp = 1;
    lastval = Rxy[Ixy[n-1]];
  }
  else{
    tmp = 0;
    lastval = -1;
  }
  Rx[Ixy[n-1]] = lastpos;
  
  for(j = n - 2; j>= 0; j--){
    if(i_perm[Ixy[j]] == 1){
      if(lastval != Rxy[Ixy[j]]){
        lastpos -= tmp;
        tmp = 0;
      }
      tmp++;
      lastval = Rxy[Ixy[j]];
      Rx[Ixy[j]] = lastpos;
    }
    else{
      if(Rxy[Ixy[j]] == Rxy[Ixy[j+1]])
        Rx[Ixy[j]] = Rx[Ixy[j+1]];
      else
        Rx[Ixy[j]] = lastpos - tmp;
    }
  }
}


void Findx(int **Rxy, int **Ixy, int *i_perm, int *n1, int *n2, int **Rx)
{
  int i, n;
  n = *n1 + *n2;
  for(i = 0; i < n; i++)
    Findx2(Rxy[i], Ixy[i], i_perm, n1, n2, Rx[i]);
}


void ranksort3(int n, int *xyidx, double *xy, int **Rxy, int **Ixy)
{
  int i, j, ileft, iright, lastpos;
  double lastval;
  for(i = 0; i < n; i++){
    lastval = -1;
    ileft = 0;
    iright = n - 1;
    j = n - 1;
    lastpos = n - 1;
    while(ileft<iright){
      if((lastval != xy[i] - xy[ileft]) && (lastval != xy[iright] - xy[i]))
        lastpos = j;
      if(ileft==i){
        lastval = xy[iright] - xy[i];
        Ixy[xyidx[i]][j] = xyidx[iright];
        Rxy[xyidx[i]][xyidx[iright]] = lastpos;
        iright--;
      }
      else if(iright==i){
        lastval = xy[i] - xy[ileft];
        Ixy[xyidx[i]][j] = xyidx[ileft];
        Rxy[xyidx[i]][xyidx[ileft]] = lastpos;
        ileft++;
      }
      else{
        if(xy[i] - xy[ileft] > xy[iright] - xy[i]){
          lastval = xy[i] - xy[ileft];
          Ixy[xyidx[i]][j] = xyidx[ileft];
          Rxy[xyidx[i]][xyidx[ileft]] = lastpos;
          ileft++;
        }
        else{
          lastval = xy[iright] - xy[i];
          Ixy[xyidx[i]][j] = xyidx[iright];
          Rxy[xyidx[i]][xyidx[iright]] = lastpos;
          iright--;
        }
      }
      j--;
    }
    Ixy[xyidx[i]][0] = xyidx[i];
    if(lastval==0)
      Rxy[xyidx[i]][xyidx[i]] = lastpos;
    else
      Rxy[xyidx[i]][xyidx[i]] = 0;
  }
}


double Ball_Divergence(int **Rxy, int **Rx, int *i_perm_tmp, int *n1, int *n2, int *weight)
{
  int i, j, n;
  double TS = 0, SS = 0, p1, p2, p3, ans;
  n = *n1 + *n2;
  double Dn1 = (*n1), Dn2 = (*n2);
  // Calculate: 
  for(i = 0; i < *n1; i++) {
    for(j = 0; j < *n1; j++) {
      p1 = Rx[i_perm_tmp[i]][i_perm_tmp[j]] + 1;
      p2 = Rxy[i_perm_tmp[i]][i_perm_tmp[j]] - p1 + 1;
      p3 = (p1 + p2) / n;
      if(p3*(1-p3)==0)
        continue;
      //TS++;
      //TS += pow(fabs(p1/(*n1)-p2/(*n2)),*gamma)/(pow(p3,*alpha)*pow(1-p3,*beta));
      ans = p1/(*n1) - p2/(*n2);
      if(*weight==0) {
        TS += (ans * ans);
      } else {
        TS += (ans * ans)*Dn1/(p1 + p2);
      }
    }
  }
  // Calculate
  for(i = *n1; i < n; i++) {
    for(j = *n1; j < n; j++) {
      p1 = Rx[i_perm_tmp[i]][i_perm_tmp[j]] + 1;
      p2 = Rxy[i_perm_tmp[i]][i_perm_tmp[j]] - p1 + 1;
      p3 = (p1 + p2) / n;
      if(p3*(1-p3)==0)
        continue;
      //SS++;
      //SS += pow(fabs(p1/(*n1)-p2/(*n2)),*gamma)/(pow(p3,*alpha)*pow(1-p3,*beta));
      ans = p1/(*n1)-p2/(*n2);
      if(*weight==0) {
        SS += (ans * ans);
      } else {
        SS += (ans * ans)*Dn2/(p1 + p2);
      }
    }
  }
  TS = TS/(1.0*(*n1)*(*n1));
  SS = SS/(1.0*(*n2)*(*n2));
  ans = TS + SS;
  return(ans);
}


double Ball_Divergence_parallel(int **Rxy, int **Rx, int *i_perm_tmp, int *n1, int *n2, int *weight, int *nthread)
{
  int n;
  double ball_divergence;
  n = *n1 + *n2;
  double Dn1 = (*n1), Dn2 = (*n2);
  double SS_value = 0.0, TS_value = 0.0;
#ifdef _OPENMP
  omp_set_num_threads(*nthread);
#endif
  #pragma omp parallel
  {
      int i, j;
      double ans;
      double TS = 0.0, SS = 0.0, p1, p2, p3;
      
      #pragma omp for
      for (i = 0; i < *n1; i++) {
// #ifdef _OPENMP
//         printf("i = %d, I am Thread %d\n", i, omp_get_thread_num());
// #endif
        for (j = 0; j < *n1; j++) {
          p1 = Rx[i_perm_tmp[i]][i_perm_tmp[j]] + 1;
          p2 = Rxy[i_perm_tmp[i]][i_perm_tmp[j]] - p1 + 1;
          p3 = (p1 + p2) / n;
          if (p3*(1 - p3) == 0)
            continue;
          ans = p1 / (*n1) - p2 / (*n2);
          if (*weight == 0) {
            TS += (ans * ans);
          }
          else {
            TS += (ans * ans)*Dn1 / (p1 + p2);
          }
        }
      }
    #pragma omp critical
    {
      TS_value += TS;
    }
    #pragma omp for
    for (i = *n1; i < n; i++) {
      for (j = *n1; j < n; j++) {
        p1 = Rx[i_perm_tmp[i]][i_perm_tmp[j]] + 1;
        p2 = Rxy[i_perm_tmp[i]][i_perm_tmp[j]] - p1 + 1;
        p3 = (p1 + p2) / n;
        if (p3*(1 - p3) == 0)
          continue;
        ans = p1 / (*n1) - p2 / (*n2);
        if (*weight == 0) {
          SS += (ans * ans);
        }
        else {
          SS += (ans * ans)*Dn2 / (p1 + p2);
        }
      }
    }
    #pragma omp critical
    {
      SS_value += SS;
    }
  }
  TS_value = TS_value/(1.0*(*n1)*(*n1));
  SS_value = SS_value/(1.0*(*n2)*(*n2));
  ball_divergence = TS_value + SS_value;
  return(ball_divergence);
}


double Ball_Divergence_wrapper(int **Rxy, int **Rx, int *i_perm_tmp, int *n1, int *n2, int *weight, int *nthread)
{
  if(*nthread > 1) {
    return Ball_Divergence_parallel(Rxy, Rx, i_perm_tmp, n1, n2, weight, nthread);
  } else {
    return Ball_Divergence(Rxy, Rx, i_perm_tmp, n1, n2, weight);
  }
}


/*
 * Two sample ball divergence test, based on permutation technique (for multivariate data)
 * input:
 * xy: vectorized distance matrix calculate used original data
 */
void BD(double *bd, double *permuted_bd, double *xy, int *n1, int *n2, int *p, int *dst, int *R, int *weight, int *nthread)
{
  //  computes TST(x,y)  	
  int i, j, n;
  /*
  Dxy: distance matrix
  Ixy: each row corresponding to sample index
  Rxy: each row corresponding the rank in each row (all data)
  Rx: each row corresponding the rank in each row (group 1)
  i_perm: group indicator (value 1 corresponding to group 1, and value 0 corresponding to group 0)
  i_perm_tmp: index of sample, value: 0, 1, 2, ..., n - 1;
  */
  int    **Ixy, **Rxy, **Rx, *i_perm, *i_perm_tmp;
  double **Dxy;
  double ans;
  
  n = *n1 + *n2;
  Dxy = alloc_matrix(n, n);
  Ixy = alloc_int_matrix(n, n);
  Rx  = alloc_int_matrix(n, n);
  Rxy = alloc_int_matrix(n, n);
  i_perm = (int *) malloc(n * sizeof(int));
  i_perm_tmp = (int *) malloc(n * sizeof(int));
  
  // get vectorize distance matrix Dxy:
  if(*dst == 1) {
    vector2matrix(xy, Dxy, n, n, 1);
  }
  else{
    Euclidean_distance(xy, Dxy, n, *p);
  }
  
  for(i = 0; i < n; i++)
    for(j = 0; j < n; j++)
      Ixy[i][j]=j;
  
  for(i = 0; i < n; i++){
    if(i<(*n1)) {
      i_perm[i] = 1;
    } else {
      i_perm[i] = 0;
    }
    i_perm_tmp[i] = i;
  }
  for(i = 0; i < n; i++){
    quicksort(Dxy[i], Ixy[i], 0, n-1);
  }
  ranksort2(n, Rxy, Dxy, Ixy);
  free_matrix(Dxy, n, n);
  Findx(Rxy, Ixy, i_perm, n1, n2, Rx);
  ans = Ball_Divergence_wrapper(Rxy, Rx, i_perm_tmp, n1, n2, weight, nthread);
  *bd = ans;
  if(*R > 0){
    for(i = 0; i < *R; i++){ 
      // stop permutation if user stop it manually:
      if(pending_interrupt()) {
        print_stop_message();
        break;
      }
      resample3(i_perm, i_perm_tmp, n, n1);
      Findx(Rxy, Ixy, i_perm, n1, n2, Rx);
      permuted_bd[i] = Ball_Divergence_wrapper(Rxy, Rx, i_perm_tmp, n1, n2, weight, nthread);
    }
  }
  // free matrix:
  free_int_matrix(Ixy, n, n);
  free_int_matrix(Rxy, n, n);
  free_int_matrix(Rx, n, n);
  free(i_perm);
  free(i_perm_tmp);
  return;
}	 


void BD_parallel(double *bd, double *permuted_bd, double *xy, int *n1, int *n2, int *p, int *dst, int *R, int *weight, int *nthread)
{
	int i, j, n;
	int    **Ixy, **Rxy, **Rx, *i_perm, *i_perm_tmp;
	double **Dxy;
	double ans;

	n = *n1 + *n2;
	Dxy = alloc_matrix(n, n);
	Ixy = alloc_int_matrix(n, n);
	Rx = alloc_int_matrix(n, n);
	Rxy = alloc_int_matrix(n, n);
	i_perm = (int *)malloc(n * sizeof(int));
	i_perm_tmp = (int *)malloc(n * sizeof(int));

	// get vectorize distance matrix Dxy:
	if (*dst == 1) {
		vector2matrix(xy, Dxy, n, n, 1);
	}

	for (i = 0; i < n; i++)
		for (j = 0; j < n; j++)
			Ixy[i][j] = j;

	for (i = 0; i < n; i++) {
		if (i<(*n1)) {
			i_perm[i] = 1;
		}
		else {
			i_perm[i] = 0;
		}
		i_perm_tmp[i] = i;
	}
	for (i = 0; i < n; i++) {
		quicksort(Dxy[i], Ixy[i], 0, n - 1);
	}
	ranksort2(n, Rxy, Dxy, Ixy);
	free_matrix(Dxy, n, n);
	Findx(Rxy, Ixy, i_perm, n1, n2, Rx);
	ans = Ball_Divergence(Rxy, Rx, i_perm_tmp, n1, n2, weight);
	*bd = ans;
	if (*R > 0) {
#ifdef _OPENMP
		omp_set_num_threads(*nthread);
#endif
		// Init parallel
#pragma omp parallel
		{
			int **Rx_thread, *i_perm_thread, *i_perm_tmp_thread;
			int k, i_thread;
			double ans = 0.0;
			Rx_thread = alloc_int_matrix(n, n);
			i_perm_thread = (int *)malloc(n * sizeof(int));
			i_perm_tmp_thread = (int *)malloc(n * sizeof(int));
#pragma omp critical
			{
				for (k = 0; k < n; k++) {
					//printf("In thread: %d\n", omp_get_thread_num());
					if (k < (*n1)) {
						i_perm_thread[k] = 1;
					}
					else {
						i_perm_thread[k] = 0;
					}
					i_perm_tmp_thread[k] = k;
					//printf("End assign\n");
				}
			}
#pragma omp for
			for (i_thread = 0; i_thread < (*R); i_thread++) {
				// stop permutation if user stop it manually:
				if (pending_interrupt()) {
					print_stop_message();
					break;
				}
				resample3(i_perm_thread, i_perm_tmp_thread, n, n1);
				Findx(Rxy, Ixy, i_perm_thread, n1, n2, Rx_thread);
				ans = Ball_Divergence(Rxy, Rx_thread, i_perm_tmp_thread, n1, n2, weight);
				permuted_bd[i_thread] = ans;
			}
		}

	}
	// free matrix:
	free_int_matrix(Ixy, n, n);
	free_int_matrix(Rxy, n, n);
	free_int_matrix(Rx, n, n);
	free(i_perm);
	free(i_perm_tmp);
	return;
}


/*
 * Two sample ball divergence test, based on permutation technique (for univariate data)
 * input:
 * xy: original data
 */
void UBD(double *bd, double *permuted_bd, double *xy, int *n1, int *n2, int *R, int *weight, int *nthread)
{
  //  computes TST(x,y)  	
  int i, j, n;
  int **Ixy, **Rxy, **Rx, *i_perm, *i_perm_tmp, *xyidx;
  double ans;
  
  n = *n1 + *n2;
  Ixy = alloc_int_matrix(n, n);
  Rx  = alloc_int_matrix(n, n);
  Rxy = alloc_int_matrix(n, n);
  i_perm = (int *) malloc(n * sizeof(int));
  i_perm_tmp = (int *) malloc(n * sizeof(int));
  xyidx = (int *) malloc(n * sizeof(int));
  for(i = 0; i < n; i++){
    xyidx[i] = i;
    for(j = 0; j < n; j++)
      Ixy[i][j]=j;
  }
  
  for(i = 0; i < n; i++){
    if(i<(*n1))
      i_perm[i] = 1;
    else
      i_perm[i] = 0;
    i_perm_tmp[i] = i;
  }
  
  quicksort(xy, xyidx, 0, n-1);
  ranksort3(n, xyidx, xy, Rxy, Ixy);
  Findx(Rxy, Ixy, i_perm, n1, n2, Rx);
  ans = Ball_Divergence_wrapper(Rxy, Rx, i_perm_tmp, n1, n2, weight, nthread);
  //printf("TST0=%f\n",ans);
  *bd = ans;
  if(*R > 0){
    for(i = 0; i < *R; i++){ 
      // stop permutation if user stop it manually:
      if(pending_interrupt()) {
        print_stop_message();
        break;
      }
      resample3(i_perm, i_perm_tmp, n, n1);
      Findx(Rxy, Ixy, i_perm, n1, n2, Rx);
      permuted_bd[i] = Ball_Divergence_wrapper(Rxy, Rx, i_perm_tmp, n1, n2, weight, nthread);
      /*
      printf("R=%d\n",i);
      printf("TST1=%f\n",ans1);
      printf("i_perm= %d",i_perm[0]);
      for(j=1;j<n;j++)
      printf(" %d",i_perm[j]);
      printf("\n");
      */
    }
  }
  
  free_int_matrix(Ixy, n, n);
  free_int_matrix(Rxy, n, n);
  free_int_matrix(Rx, n, n);
  free(i_perm);
  free(i_perm_tmp);
  free(xyidx);
  return;
}


void UBD_parallel(double *bd, double *permuted_bd, double *xy, int *n1, int *n2, int *R, int *weight, int *nthread)
{
	//  computes TST(x,y)  	
	int i, j, n;
	int **Ixy, **Rxy, **Rx, *i_perm, *i_perm_tmp, *xyidx;
	double ans;

	n = *n1 + *n2;
	Ixy = alloc_int_matrix(n, n);
	Rx = alloc_int_matrix(n, n);
	Rxy = alloc_int_matrix(n, n);
	i_perm = (int *) malloc(n * sizeof(int));
	i_perm_tmp = (int *) malloc(n * sizeof(int));
	xyidx = (int *) malloc(n * sizeof(int));
	for (i = 0; i < n; i++) {
		xyidx[i] = i;
		for (j = 0; j < n; j++)
			Ixy[i][j] = j;
	}

	for (i = 0; i < n; i++) {
		if (i<(*n1))
			i_perm[i] = 1;
		else
			i_perm[i] = 0;
		i_perm_tmp[i] = i;
	}

	quicksort(xy, xyidx, 0, n - 1);
	ranksort3(n, xyidx, xy, Rxy, Ixy);
	Findx(Rxy, Ixy, i_perm, n1, n2, Rx);
	ans = Ball_Divergence(Rxy, Rx, i_perm_tmp, n1, n2, weight);	
	*bd = ans;
	if (*R > 0) {
#ifdef _OPENMP
		omp_set_num_threads(*nthread);
#endif
		#pragma omp parallel
		{
			int **Rx_thread, *i_perm_thread, *i_perm_tmp_thread;
			int k, i_thread;
			double ans = 0.0;
			Rx_thread = alloc_int_matrix(n, n);
			i_perm_thread = (int *) malloc(n * sizeof(int));
			i_perm_tmp_thread = (int *) malloc(n * sizeof(int));
			#pragma omp critical
			{
				for (k = 0; k < n; k++) {
					//printf("In thread: %d\n", omp_get_thread_num());
					if (k < (*n1)) {
						i_perm_thread[k] = 1;
					}
					else {
						i_perm_thread[k] = 0;
					}
					i_perm_tmp_thread[k] = k;
					//printf("End assign\n");
				}
			}
			#pragma omp for
			for (i_thread = 0; i_thread < (*R); i_thread++) {
				// stop permutation if user stop it manually:
				if (pending_interrupt()) {
					print_stop_message();
					break;
				}
				resample3(i_perm_thread, i_perm_tmp_thread, n, n1);
				Findx(Rxy, Ixy, i_perm_thread, n1, n2, Rx_thread);
				ans = Ball_Divergence(Rxy, Rx_thread, i_perm_tmp_thread, n1, n2, weight);
				permuted_bd[i_thread] = ans;
			}
		}
	}
	free_int_matrix(Ixy, n, n);
	free_int_matrix(Rxy, n, n);
	free_int_matrix(Rx, n, n);
	free(i_perm);
	free(i_perm_tmp);
	free(xyidx);
	return;
}


/*
 * return ball divergence value of two sample with sample size n1, n2;
 * permutation procedure not consider at all
 * input:
 * xy: vectorized distance matrix
 */
double bd_value(double *xy, int *n1, int *n2, int *weight, int *nthread)
{
  int i, j, n;
  int    **Ixy, **Rxy, **Rx, *i_perm, *i_perm_tmp;
  double **Dxy;
  double ans;
  
  n = *n1 + *n2;
  Dxy = alloc_matrix(n, n);
  Ixy = alloc_int_matrix(n, n);
  Rx  = alloc_int_matrix(n, n);
  Rxy = alloc_int_matrix(n, n);
  i_perm = (int *) malloc(n * sizeof(int));
  i_perm_tmp = (int *) malloc(n * sizeof(int));
  vector2matrix(xy, Dxy, n, n, 1);
  for(i = 0; i < n; i++)
    for(j = 0; j < n; j++)
      Ixy[i][j]=j;
  
  for(i = 0; i < n; i++){
    if(i<(*n1))
      i_perm[i] = 1;
    else
      i_perm[i] = 0;
    i_perm_tmp[i] = i;
  }
  for(i = 0; i < n; i++){
    quicksort(Dxy[i], Ixy[i], 0, n-1);
  }
  ranksort2(n, Rxy, Dxy, Ixy);
  free_matrix(Dxy, n, n);
  Findx(Rxy, Ixy, i_perm, n1, n2, Rx);
  ans = Ball_Divergence_wrapper(Rxy, Rx, i_perm_tmp, n1, n2, weight, nthread);
  // free matrix
  free_int_matrix(Ixy, n, n);
  free_int_matrix(Rxy, n, n);
  free_int_matrix(Rx, n, n);
  free(i_perm);
  free(i_perm_tmp);
  return(ans);
}


// input:
// size: an array contain sample size in each group
// cumulate_size: the cumulative sums of size vector  
// k: group number
void compute_cumulate_size(int *cumulate_size, int *size, int *k)
{
  int i;
  for(i = 0; i < (*k); i++) {
    if(i == 0) {
      cumulate_size[i] = 0;
    } else {
      cumulate_size[i] = cumulate_size[i-1] + size[i-1];
    }
  }
  return;
}


// need to fetch submatrix and vectorize it
// xy: vectorize distance matrix
// ij_dst: vectorized distance matrix of group i and group j
// cumulate_size: the cumulative sums of size vector  
// size: an array contain sample size in each group
// *p < *q is needed and *p == i, *q == j
void get_ij_dst(double *xy, double *ij_dst, int *cumulate_size, int *size, int *n, int *p, int *q) 
{
//void get_ij_dst(double *xy, double *ij_dst, int *cumulate_size, int *size, int *n, int i, int j) {
  int k = 0;
  int k1 = 0;
  int k2 = 0;
  int i = *p;
  int j = *q;
  int n1 = size[i];
  int n2 = size[j];
  int num = n1 + n2;
  // for group1:
  int index_ii = (*n)*cumulate_size[i] + cumulate_size[i];
  int index_ij = (*n)*cumulate_size[i] + cumulate_size[j];
  for(k1 = 0; k1 < n1; k1++) {
    for(k2 = 0; k2 < num; k2++) {
      if(k2 < n1) {
        ij_dst[k] = xy[index_ii + k2];
      } else {
        ij_dst[k] = xy[index_ij + k2 - n1];
      }
      k = k + 1;
    }
    index_ii = index_ii + (*n);
    index_ij = index_ij + (*n);
  }
  // for group2:
  index_ii = (*n)*cumulate_size[j] + cumulate_size[i];
  index_ij = (*n)*cumulate_size[j] + cumulate_size[j];
  for(k1 = 0; k1 < n2; k1++) {
    for(k2 = 0; k2 < num; k2++) {
      if(k2 < n1) {
        ij_dst[k] = xy[index_ii + k2];
      } else {
        ij_dst[k] = xy[index_ij + k2 - n1];
      }
      k = k + 1;
    }
    index_ii = index_ii + (*n);
    index_ij = index_ij + (*n);
  }
  return;
}



// return k-sample ball divergence
// input: 
// xy: vectorized distance matrix
// size: an array contain sample size in each group
// n: sample size
// k: group number
// weight: if weight == TRUE, weight BD will be returned
// kbd: K-samples ball divergence statistic
// void kbd_value(double *xy, int *size, int *n, int *k, int *weight, double *kbd) {
double kbd_value(double *xy, int *size, int *n, int *k, int *weight, int *nthread)
{
  int K = *k;
  int i, j;
  int *cumulate_size;
  double *ij_dst;
  int two_group_size;
  int tmp_value1 = 0, tmp_value2 = 0;
  int *n1 = &tmp_value1;
  int *n2 = &tmp_value2;
  double kbd_stat_value = 0.0, bd_stat_value = 0.0;
  cumulate_size = (int *) malloc(K * sizeof(int));
  compute_cumulate_size(cumulate_size, size, k);
  for(i = 0; i < K; i++) {
    for(j = (i + 1); j < K; j++) {
      tmp_value1 = size[i]; 
      tmp_value2 = size[j]; 
      two_group_size = tmp_value1 + tmp_value2;
      two_group_size = two_group_size * two_group_size;
      ij_dst = (double *) malloc(two_group_size * sizeof(double)); 
      get_ij_dst(xy, ij_dst, cumulate_size, size, n, &i, &j);
      bd_stat_value =  bd_value(ij_dst, n1, n2, weight, nthread);
      kbd_stat_value = kbd_stat_value + bd_stat_value;
      free(ij_dst);
    }
  }
  free(cumulate_size);
  // *kbd = kbd_stat_value;
  // Rprintf("kbd value:%f;     ", kbd_stat_value);
  return(kbd_stat_value);
}


// return a distance matrix after permute index is given
// input: 
// xy: vectorized distance matrix
// new_xy: vectorized distance matrix store after sample was permuted
// index: permutation index
// N: sample size
void permute_dst(double *xy, double *new_xy, int *index, int *N)
{
  int n = (*N);
  int row_index = 0, col_index = 0, exchange_index1 = 0, exchange_index2 = 0;
  for(int j = 0; j < (*N); j++) {
    for(int k = 0; k < (*N); k++) {
      row_index = index[j];
      col_index = index[k];
      exchange_index2 = row_index*n + col_index;
      new_xy[exchange_index1] = xy[exchange_index2];
      exchange_index1 = exchange_index1 + 1;
    }
  }
  return;
}


// R function call this function to implement ball divergence based k-sample test
// if R = 0, k-sample ball divergence statistic will be returned, else, k-sample test p-value 
// base on ball divergence statistic will be return
// Input Arugement: 
// kbd: K-samples ball divergence statistic or p-value
// permuted_kbd: K-samples ball divergence statistic after permutation
// xy: vectorized distance matrix
// size: an array contain sample size in each group
// n: sample size
// k: group number
// weight: if weight == TRUE, weight BD will be returned
void KBD(double *kbd, double *permuted_kbd, double *xy, int *size, int *n, int *k, int *R, int *weight, int *nthread)
{
  double ans = 0.0;
  int N = (*n);
  int dst_size = N*N;
  int i_permute;
  ans = kbd_value(xy, size, n, k, weight, nthread);
  *kbd = ans;
  if((*R)>0) {
    int *index;
    double *new_xy;
    new_xy = (double *) malloc(dst_size * sizeof(double));
    index = (int *) malloc(N * sizeof(int));
    for(i_permute = 0; i_permute < N; i_permute++) {
      index[i_permute] = i_permute;
    }
    for(int j = 0; j < (*R); j++) {
      // stop permutation if user stop it manually:
      if(pending_interrupt()) {
        print_stop_message();
        break;
      }
      // permute data index:
      shuffle(index, n);
      // adjust vectorized distance matrix according to permuted index: 
      permute_dst(xy, new_xy, index, n);
      // K-sample BD after permutation:
      permuted_kbd[j] = kbd_value(new_xy, size, n, k, weight, nthread);
    }
    free(new_xy);
    free(index);
  }
  return;
}


/*
 * return ball divergence value of two sample(univariate) with sample size n1, n2;
 * permutation procedure not consider at all
 * input:
 * xy: value of sample
 * n1, n2: sample sizes of group1 and group2
 */
double ubd_value(double *xy, int *n1, int *n2, int *weight, int *nthread)
{
  //  computes TST(x,y)  	
  int i, j, n;
  int **Ixy, **Rxy, **Rx, *i_perm, *i_perm_tmp, *xyidx;
  double ans;
  
  n = *n1 + *n2;
  Ixy = alloc_int_matrix(n, n);
  Rx  = alloc_int_matrix(n, n);
  Rxy = alloc_int_matrix(n, n);
  i_perm = (int *) malloc(n * sizeof(int));
  i_perm_tmp = (int *) malloc(n * sizeof(int));
  xyidx = (int *) malloc(n * sizeof(int));
  for(i = 0; i < n; i++){
    xyidx[i] = i;
    for(j = 0; j < n; j++)
      Ixy[i][j]=j;
  }
  
  for(i = 0; i < n; i++){
    if(i<(*n1))
      i_perm[i] = 1;
    else
      i_perm[i] = 0;
    i_perm_tmp[i] = i;
  }
  
  quicksort(xy, xyidx, 0, n-1);
  ranksort3(n, xyidx, xy, Rxy, Ixy);
  Findx(Rxy, Ixy, i_perm, n1, n2, Rx);
  ans = Ball_Divergence_wrapper(Rxy, Rx, i_perm_tmp, n1, n2, weight, nthread);
  // free memory:
  free_int_matrix(Ixy, n, n);
  free_int_matrix(Rxy, n, n);
  free_int_matrix(Rx, n, n);
  free(i_perm);
  free(i_perm_tmp);
  free(xyidx);
  return ans;
}

/*
 * Given input data xy, fetch the i,j group data and insert into ij_value
 * input:
 * xy: value of all sample
 * ij_value: contain sample value of sample i and j
 * size: sample size of each group
 * cumulate_size: the cumulative sums of size vector
 * (for example, size: [1, 2, 3] ==> cumulate_size: [0, 1, 3])
 * p point the i group
 * q point the j group
 */
void get_ij_value(double *xy, double *ij_value, int *cumulate_size, int *size, int *p, int *q)
{
  int k1 = 0, k2 = 0;
  int i = *p;
  int j = *q;
  int n1 = size[i];
  int n2 = size[j];
  // for group1:
  int start_index_i = cumulate_size[i];
  int start_index_j = cumulate_size[j];
  for(k1 = 0; k1 < n1; k1++) {
    ij_value[k1] = xy[start_index_i + k1];
  }
  // for group2:
  for(k2 = 0; k2 < n2; k2++) {
    ij_value[k1] = xy[start_index_j + k2];
    k1 += 1;
  }
  return;
}


/*
 * Given input univariate data xy and sample size of each group, compute kbd;
 * input:
 * xy: value of all sample
 * size: sample size of each group
 * cumulate_size: the cumulative sums of size vector
 * (for example, size: [1, 2, 3] ==> cumulate_size: [0, 1, 3])
 * p point the i group
 * q point the j group
 */
double ukbd_value(double *xy, int *size, int *k, int *weight, int *nthread)
{
  int K = *k;
  int i, j;
  int *cumulate_size;
  int two_group_size;
  int tmp_value1 = 0, tmp_value2 = 0;
  int *n1 = &tmp_value1;
  int *n2 = &tmp_value2;
  double *ij_value;
  double kbd_stat_value = 0.0, bd_stat_value = 0.0;
  cumulate_size = (int *) malloc(K * sizeof(int));
  compute_cumulate_size(cumulate_size, size, k);
  for(i = 0; i < K; i++) {
    for(j = (i + 1); j < K; j++) {
      tmp_value1 = size[i]; 
      tmp_value2 = size[j]; 
      two_group_size = tmp_value1 + tmp_value2;
      ij_value = (double *) malloc(two_group_size * sizeof(double)); 
      get_ij_value(xy, ij_value, cumulate_size, size, &i, &j);
      bd_stat_value =  ubd_value(ij_value, n1, n2, weight, nthread);
      kbd_stat_value = kbd_stat_value + bd_stat_value;
      free(ij_value);
    }
  }
  free(cumulate_size);
  return(kbd_stat_value);
}


// R function call this function to implement ball divergence based k-sample test for 
// univariate random variable
// if R = 0, k-sample ball divergence statistic will be returned, else, k-sample test p-value 
// base on ball divergence statistic will be return
// Input Arugement: 
// kbd: K-samples ball divergence statistic
// permuted_kbd: K-samples ball divergence statistic after permutation
// xy: vectorized distance matrix
// size: an array contain sample size in each group
// n: sample size
// k: group number
// weight: if weight == TRUE, weight BD will be returned
void UKBD(double *kbd, double *permuted_kbd, double *xy, int *size, int *n, int *k, int *R, int *weight, int *nthread)
{
  double ans = 0.0;
  ans = ukbd_value(xy, size, k, weight, nthread);
  *kbd = ans;
  if((*R)>0) {
    // init permutation test:
    for(int j = 0; j < (*R); j++) {
      // stop permutation loop:
      if(pending_interrupt()) {
        print_stop_message();
        break;
      }
      // permute data value:
      shuffle_value(xy, n);
      // K-sample UBD after permutation:
      permuted_kbd[j] = ukbd_value(xy, size, k, weight, nthread);
    }
  }
  return;
}

/////////////////////////////////////////////////////////////////
//////////// four functions below wrap key function /////////////
// bd_stat, kbd_stat return ball divergence statistic 
// bd_stat, kbd_stat execute ball divergence statistic based test
/////////////////////////////////////////////////////////////////

/*
 * R function call this function to compute two sample ball divergence statistic
 * if k = 2, then n1 = size[1] and n2 = size[2];
 */
void bd_stat(double *bd, double *xy, int *size, int *n, int *k, int *weight, int *dst, int *nthread)
{
  double ans = 0.0;
  int n1 = 0, n2 = 0;
  if((*k) == 2) {
    n1 = size[0];
    n2 = size[1];
    if((*dst)) {
      ans = bd_value(xy, &n1, &n2, weight, nthread);
    } else {
      ans = ubd_value(xy, &n1, &n2, weight, nthread);
    }
  } else {
    if((*dst)) {
      ans = kbd_value(xy, size, n, k, weight, nthread);
    } else {
      ans = ukbd_value(xy, size, k, weight, nthread);
    }
  }
  *bd = ans;
  return;
}


/*
 * R function call this function to compute K sample ball divergence statistic
 */
/*
void kbd_stat(double *bd, double *xy, int *size, int *n, int *k, int *weight, int *dst, int *nthread)
{
  double ans = 0.0;
  if((*dst)) {
    ans = kbd_value(xy, size, n, k, weight, nthread);
  } else {
    ans = ukbd_value(xy, size, k, weight, nthread);
  }
  *bd = ans;
  return;
}
*/

/*
 * R function call this function to execute Two sample ball divergence test
 */
void bd_test(double *bd, double *permuted_bd, double *xy, int *size, int *n, int *k, int *dst, int *R, int *weight, int *nthread)
{
  int n1 = 0, n2 = 0;
  int p = 1;
  //parallel method
  // if parallel_type == 1, we parallel the computation through statistics.
  // if parallel_type == 2, we parallel the computation through permutation.
  int parallel_type = 2;
  if (((*n) > 500)) {
	  parallel_type = 1;
  }
  if ((*R) < 300) {
	  *nthread = 1;
  }
  //
  if((*k) == 2) {
    n1 = size[0];
    n2 = size[1];
    if(*dst) {
		if ((parallel_type) == 2) {
			BD_parallel(bd, permuted_bd, xy, &n1, &n2, &p, dst, R, weight, nthread);
		}
		else {
			BD(bd, permuted_bd, xy, &n1, &n2, &p, dst, R, weight, nthread); 
		}
    } else {
		if ((parallel_type) == 2) {
			UBD_parallel(bd, permuted_bd, xy, &n1, &n2, R, weight, nthread);
		}
		else {
			UBD(bd, permuted_bd, xy, &n1, &n2, R, weight, nthread);
		}
    } 
  } else {
    if(*dst) {
      KBD(bd, permuted_bd, xy, size, n, k, R, weight, nthread);
    } else {
      UKBD(bd, permuted_bd, xy, size, n, k, R, weight, nthread);
    }
  }
  return;
}


/*
 * R function call this function to execute Two sample ball divergence test
 */
/*
void kbd_test(double *kbd, double *permuted_kbd, double *xy, int *size, int *n, int *k, int *dst, int *R, int *weight, int *nthread)
{
  if(*dst) {
    KBD(kbd, permuted_kbd, xy, size, n, k, R, weight, nthread);
  } else {
    UKBD(kbd, permuted_kbd, xy, size, n, k, R, weight, nthread);
  }
  return;
}
*/