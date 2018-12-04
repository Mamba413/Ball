#include "Ball_omp.h"
#include "stdio.h"
#include "stdlib.h"
#include "math.h"
#include "utilities.h"
#include "R.h"

double Ball_Divergence_Value(int **Rxy, int **Rx, int *i_perm_tmp, int *n1, int *n2) {
  int i, j, n;
  double TS_weight0 = 0, SS_weight0 = 0;
  double p1, p2, p3, ans;
  n = *n1 + *n2;
  
  for (i = 0; i < *n1; i++) {
    for (j = 0; j < *n1; j++) {
      p1 = Rx[i_perm_tmp[i]][i_perm_tmp[j]] + 1;
      p2 = Rxy[i_perm_tmp[i]][i_perm_tmp[j]] - p1 + 1;
      p3 = (p1 + p2) / n;
      if (p3 * (1 - p3) == 0) { continue; }
      ans = p1 / (*n1) - p2 / (*n2);
      TS_weight0 += (ans * ans);
    }
  }
  for (i = *n1; i < n; i++) {
    for (j = *n1; j < n; j++) {
      p1 = Rx[i_perm_tmp[i]][i_perm_tmp[j]] + 1;
      p2 = Rxy[i_perm_tmp[i]][i_perm_tmp[j]] - p1 + 1;
      p3 = (p1 + p2) / n;
      if (p3 * (1 - p3) == 0) { continue; }
      ans = p1 / (*n1) - p2 / (*n2);
      SS_weight0 += (ans * ans);
    }
  }
  return TS_weight0 / (1.0 * (*n1) * (*n1)) + SS_weight0 / (1.0 * (*n2) * (*n2));
}

void vector2matrix_int(int *x, int **y, int N, int d, int isroworder)
{
  /* copy a d-variate sample into a matrix, N samples in rows */
  int i, k;
  if (isroworder == 1) {
    for (k = 0; k<d; k++)
      for (i = 0; i<N; i++)
        y[i][k] = (*(x + i*d + k));
  }
  else {
    for (k = 0; k<N; k++)
      for (i = 0; i<d; i++)
        y[i][k] = (*(x + k*N + i));
  }
  return;
}

int find_unique_group_num_index(int n1, int unique_group_num_length, int* unique_group_num_vec) {
  int i;
  for (i = 0; i < unique_group_num_length; i++)
  {
    if (n1 == unique_group_num_vec[i]) { break; }
  }
  return i;
}

void random_index_vec(int* randn_vec, int n)
{
  for (int i = n - 1; i > 0; --i) {
    randn_vec[i - 1] = ((int) (RAND_MAX * unif_rand())) % (i + 1);
  }
}

void findx2_gwas(int *Rxy, int *Ixy, int *i_perm, int n1, int n, int *Rx)
{
  int j, lastpos, lastval, tmp;
  
  lastpos = n1 - 1;
  Rx[Ixy[n - 1]] = lastpos;
  if (i_perm[Ixy[n - 1]] == 1) {
    tmp = 1;
    lastval = Rxy[Ixy[n - 1]];
  }
  else {
    tmp = 0;
    lastval = -1;
  }
  
  for (j = n - 2; j >= 0; j--) {
    if (i_perm[Ixy[j]] == 1) {
      if (lastval != Rxy[Ixy[j]]) {
        lastpos -= tmp;
        tmp = 0;
      }
      tmp++;
      lastval = Rxy[Ixy[j]];
      Rx[Ixy[j]] = lastpos;
    }
    else {
      if (Rxy[Ixy[j]] == Rxy[Ixy[j + 1]])
        Rx[Ixy[j]] = Rx[Ixy[j + 1]];
      else
        Rx[Ixy[j]] = lastpos - tmp;
    }
  }
}

void findx_gwas(int **Rxy, int **Ixy, int *i_perm, int n1, int n, int **Rx)
{
  for (int i = 0; i < n; i++)
    findx2_gwas(Rxy[i], Ixy[i], i_perm, n1, n, Rx[i]);
}

void find_i_perm_temp_gwas(int *snp, int *i_perm_temp, int n1, int n)
{
  int k = 0, l = 0;
  for (int i = 0; i < n; i++) {
    if (snp[i] == 1) {
      i_perm_temp[k++] = i;
    } else {
      i_perm_temp[n1 + l] = i;
      l++;
    }
  }
}

void permutation_gwas(int* randn_vec, int *i_perm, int *i_perm_tmp, int n, int *n1)
{
  int i, j, temp, tmp0, tmp1;
  
  for (i = n - 2; i >= 0; --i) {
    j = randn_vec[i];
    temp = i_perm[j];
    i_perm[j] = i_perm[i];
    i_perm[i] = temp;
  }
  
  tmp0 = 0;
  tmp1 = 0;
  for (i = 0; i < n; i++) {
    if (i_perm[i] == 1) {
      i_perm_tmp[tmp0++] = i;
    } else {
      i_perm_tmp[*n1 + tmp1] = i;
      tmp1++;
    }
  }
}

/*
* @param int** snp: snp-number X sample-size
* @param int* n1_num_vec: the sample size of group 1 with snp value 0. 
* Note: in high-level, must pre-process SNP data to obtain this vector.
*/
void bd_gwas(double* stats_vec, double** permuted_stat_mat, double* pvalue_vec,
             int** snp, int* unique_group_num_vec, int* unique_group_num_length,
             int* n1_num_vec, int* snp_number, double* xy, int* r, int n)
{
  int i, j, n2;
  int **Ixy, **Rxy, **Rx, *i_perm_tmp;
  double **Dxy;
  
  Dxy = alloc_matrix(n, n);
  vector2matrix(xy, Dxy, n, n, 1);
  Ixy = alloc_int_matrix(n, n);
  for (i = 0; i < n; i++) {
    for (j = 0; j < n; j++) {
      Ixy[i][j] = j;
    }
  }
  for (i = 0; i < n; i++) {
    quicksort(Dxy[i], Ixy[i], 0, n - 1);
  }
  Rxy = alloc_int_matrix(n, n);
  ranksort2(n, Rxy, Dxy, Ixy);
  free_matrix(Dxy, n, n);
  Rx = alloc_int_matrix(n, n);
  i_perm_tmp = (int *)malloc(n * sizeof(int));
  
  // compute ball divergence for multiple snp:
  for (i = 0; i < *snp_number; i++) {
    findx_gwas(Rxy, Ixy, snp[i], n1_num_vec[i], n, Rx);
    find_i_perm_temp_gwas(snp[i], i_perm_tmp, n1_num_vec[i], n);
    n2 = n - n1_num_vec[i];
    stats_vec[i] = Ball_Divergence_Value(Rxy, Rx, i_perm_tmp, &n1_num_vec[i], &n2);
  }
  free_int_matrix(Rx, n, n);
  free(i_perm_tmp);
  // printf("statistic computation done!\n");
  
  // permutation
  if (*r > 0) {
    int **i_perm_thread, **i_perm_tmp_thread;
    i_perm_thread = alloc_int_matrix(*unique_group_num_length, n);
    i_perm_tmp_thread = alloc_int_matrix(*unique_group_num_length, n);
    for (i = 0; i < *unique_group_num_length; i++) {
      for (j = 0; j < n; j++) {
        i_perm_thread[i][j] = j < unique_group_num_vec[i] ? 1 : 0;
        i_perm_tmp_thread[i][j] = j;
      }
    }
    
    int batch_num, fix_batch_num;
    fix_batch_num = 200;
    int **randn_mat;
    randn_mat = alloc_int_matrix(fix_batch_num, n - 1);
    
    // permutation replication:
    int k_thread = 0, r_res = *r;
    while (r_res > 0) {
      // random exchange method
      batch_num = r_res >= fix_batch_num ? fix_batch_num : r_res;
      for (i = 0; i < batch_num; i++) {
        random_index_vec(randn_mat[i], n);
      }
#pragma omp parallel
      {
          int i_thread, **Rx_thread;
          Rx_thread = alloc_int_matrix(n, n);
#pragma omp for
          for (i_thread = 0; i_thread < *unique_group_num_length; i_thread++) {
            int l_thread = 0;
            int r_thread = k_thread * fix_batch_num;
            int n1_thread = unique_group_num_vec[i_thread];
            int n2_thread = n - n1_thread;
            while (l_thread < batch_num)
            {
              permutation_gwas(randn_mat[l_thread], i_perm_thread[i_thread],
                               i_perm_tmp_thread[i_thread], n, &n1_thread);
              findx_gwas(Rxy, Ixy, i_perm_thread[i_thread], n1_thread, n, Rx_thread);
              permuted_stat_mat[i_thread][r_thread] = Ball_Divergence_Value(Rxy, Rx_thread,
                                                                            i_perm_tmp_thread[i_thread], &n1_thread, &n2_thread);
              l_thread++;
              r_thread++;
            }
          }
          free_int_matrix(Rx_thread, n, n);
      }
      // printf("epoch: %d, batch: %d, done!\n", k_thread, batch_num);
      k_thread++;
      r_res -= batch_num;
    }
    // printf("Permutation done!\n");
    free_int_matrix(i_perm_thread, *unique_group_num_length, n);
    free_int_matrix(i_perm_tmp_thread, *unique_group_num_length, n);
    free_int_matrix(randn_mat, batch_num, n - 1);
    
    // compute p_value
    int permuated_mat_ind;
    for (i = 0; i < *snp_number; i++) {
      permuated_mat_ind = find_unique_group_num_index(n1_num_vec[i], *unique_group_num_length, unique_group_num_vec);
      pvalue_vec[i] = compute_pvalue(stats_vec[i], permuted_stat_mat[permuated_mat_ind], *r);
    }
  }
  free_int_matrix(Rxy, n, n);
  free_int_matrix(Ixy, n, n);
}


void ubd_gwas() {
  
}

void bd_gwas_test(double* stats_vec, double* permuted_stat_vec, double* pvalue_vec,
                  int* snp, int* unique_group_num_vec, int* unique_group_num_length, 
                  int* n1_num_vec, int* snp_number, double* xy, int *r, int *n,
                  int* dst, int* nthread) 
{
#ifdef Ball_OMP_H_
  omp_set_num_threads(*nthread);
#endif
  int **snp_data;
  snp_data = alloc_int_matrix(*snp_number, *n);
  vector2matrix_int(snp, snp_data, *snp_number, *n, 1);
  double **permuted_stat_mat;
  permuted_stat_mat = alloc_matrix(*unique_group_num_length, *r);
  if (*dst) {
    bd_gwas(stats_vec, permuted_stat_mat, pvalue_vec, 
            snp_data, unique_group_num_vec, unique_group_num_length,
            n1_num_vec, snp_number, xy, r, *n);
  } else {
    ubd_gwas();
  }
  int k = 0;
  for (int i = 0; i < *unique_group_num_length; i++) {
    for (int j = 0; j < *r; j++) {
      permuted_stat_vec[k++] = permuted_stat_mat[i][j];
    }
  }
  free_int_matrix(snp_data, *snp_number, *n);
  free_matrix(permuted_stat_mat, *unique_group_num_length, *r);
}