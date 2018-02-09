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
#include "math.h"
#include "stdlib.h"
#include "stdio.h"
#include "utilities.h"
#include "BD.h"


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
    if(i_perm[Ixy[j]]==1){
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



/*
 * permute group index: i_perm
 */
void resample3(int *i_perm, int *i_perm_tmp, int n, int *n1)
{
  int i, j, temp, tmp0, tmp1;
  
  // permute step:
  for (i = n - 1; i > 0; --i) {
    j = rand() % (i + 1);
    //j = r_available_rand() % (i + 1);
    temp = i_perm[j];
    i_perm[j] = i_perm[i];
    i_perm[i] = temp;
  }
  
  tmp0 = 0;
  tmp1 = 0;
  for(i = 0; i < n; i++)
    if(i_perm[i]==1){
      i_perm_tmp[tmp0++] = i;
    }
    else{
      i_perm_tmp[*n1 + tmp1] = i;
      tmp1++;
    }
}

void BD(double *bd, double *permuted_bd, double *xy, int *n1, int *n2, int *p, int *dst, int *R, int *weight)
{
  //  computes TST(x,y)  	
  int i, j, n;
  /*
   i_perm_tmp: index of sample, value: 0, 1, 2, ..., n - 1;
   i_perm: group index of sample, 0 or 1;
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
  ans = Ball_Divergence(Rxy, Rx, i_perm_tmp, n1, n2, weight);
  *bd = ans;
  if(*R > 0){
    for(i = 0; i < *R; i++){ 
      // stop permutation if user stop it manually:
      /*if(pending_interrupt()) {
        printf("Process stop due to user interruption! \n");
        break;
      }*/
      resample3(i_perm, i_perm_tmp, n, n1);
      Findx(Rxy, Ixy, i_perm, n1, n2, Rx);
      permuted_bd[i] = Ball_Divergence(Rxy, Rx, i_perm_tmp, n1, n2, weight);
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


void UBD(double *bd, double *permuted_bd, double *xy, int *n1, int *n2, int *R, int *weight)
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
  ans = Ball_Divergence(Rxy, Rx, i_perm_tmp, n1, n2, weight);
  //printf("TST0=%f\n",ans);
  *bd = ans;
  if(*R > 0){
    for(i = 0; i < *R; i++){ 
      // stop permutation if user stop it manually:
      /*if(pending_interrupt()) {
        printf("Process stop due to user interruption! \n");
        break;
      }*/
      resample3(i_perm, i_perm_tmp, n, n1);
      Findx(Rxy, Ixy, i_perm, n1, n2, Rx);
      permuted_bd[i] = Ball_Divergence(Rxy, Rx, i_perm_tmp, n1, n2, weight);
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


/*
 * return ball divergence value of two sample with sample size n1, n2;
 * permutation procedure not consider at all
 * input:
 * xy: vectorized distance matrix
 */
double bd_value(double *xy, int *n1, int *n2, int *weight)
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
  ans = Ball_Divergence(Rxy, Rx, i_perm_tmp, n1, n2, weight);
  // printf("%f", ans);
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
double kbd_value(double *xy, int *size, int *n, int *k, int *weight)
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
      bd_stat_value =  bd_value(ij_dst, n1, n2, weight);
      kbd_stat_value = kbd_stat_value + bd_stat_value;
      free(ij_dst);
    }
  }
  free(cumulate_size);
  // *kbd = kbd_stat_value;
  // printf("kbd value:%f;     ", kbd_stat_value);
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
      // printf("index1:%d;   ", exchange_index1);
      // printf("index2:%d \n", exchange_index2);
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
void KBD(double *kbd, double *permuted_kbd, double *xy, int *size, int *n, int *k, int *R, int *weight)
{
  double ans = 0.0;
  int N = (*n);
  int dst_size = N*N;
  int i_permute;
  ans = kbd_value(xy, size, n, k, weight);
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
     /* if(pending_interrupt()) {
        printf("Process stop due to user interruption! \n");
        break;
      }*/
      // permute data index:
      shuffle(index, n);
      // adjust vectorized distance matrix according to permuted index: 
      permute_dst(xy, new_xy, index, n);
      // K-sample BD after permutation:
      permuted_kbd[j] = kbd_value(new_xy, size, n, k, weight);
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
double ubd_value(double *xy, int *n1, int *n2, int *weight)
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
  ans = Ball_Divergence(Rxy, Rx, i_perm_tmp, n1, n2, weight);
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
double ukbd_value(double *xy, int *size, int *k, int *weight)
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
      bd_stat_value =  ubd_value(ij_value, n1, n2, weight);
      kbd_stat_value = kbd_stat_value + bd_stat_value;
      free(ij_value);
    }
  }
  free(cumulate_size);
  // printf("kbd value:%f \n", kbd_stat_value);
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
void UKBD(double *kbd, double *permuted_kbd, double *xy, int *size, int *n, int *k, int *R, int *weight)
{
  double ans = 0.0;
  ans = ukbd_value(xy, size, k, weight);
  *kbd = ans;
  if((*R)>0) {
    // init permutation test:
    for(int j = 0; j < (*R); j++) {
      // stop permutation loop:
     /* if(pending_interrupt()) {
        printf("Process stop due to user interruption! \n");
        break;
      }*/
      // permute data value:
      shuffle_value(xy, n);
      // K-sample UBD after permutation:
      permuted_kbd[j] = ukbd_value(xy, size, k, weight);
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
 */
/*void bd_stat(double *bd, double *xy, int *n1, int *n2, int *weight, int *dst)
{
  double ans = 0.0;
  if((*dst)) {
    ans = bd_value(xy, n1, n2, weight);
  } else {
    ans = ubd_value(xy, n1, n2, weight);
  }
  *bd = ans;
  return;
}
*/

/*
 * R function call this function to compute K sample ball divergence statistic
 */
void kbd_stat(double *bd, double *xy, int *size, int *n, int *k, int *weight, int *dst)
{
  double ans = 0.0;
  if((*dst)) {
    ans = kbd_value(xy, size, n, k, weight);
  } else {
    ans = ukbd_value(xy, size, k, weight);
  }
  *bd = ans;
  return;
}


/*
 * R function call this function to execute Two sample ball divergence test
 */
//void bd_test(double *bd, double *permuted_bd, double *xy, int *n1, int *n2, int *p, int *dst, int *R, int *weight)
//{
//  if(*dst) {
//    BD(bd, permuted_bd, xy, n1, n2, p, dst, R, weight);
//  } else {
//    UBD(bd, permuted_bd, xy, n1, n2, R, weight);
//  }
//  return;
//}


/*
 * R function call this function to execute Two sample ball divergence test
 */
void kbd_test(double *kbd, double *permuted_kbd, double *xy, int *size, int *n, int *k, int *dst, int *R, int *weight)
{
  if(*dst) {
    KBD(kbd, permuted_kbd, xy, size, n, k, R, weight);
  } else {
    UKBD(kbd, permuted_kbd, xy, size, n, k, R, weight);
  }
  return;
}