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
#include "R.h"
#include "Rmath.h"
#include "stdio.h"
#include "utilities.h"
#include "BI.h"


void Merge(int *permutation, int *source, int *inversion_count, int dim, int n)
{
  int *left=(int *) malloc(n * sizeof(int));
  int *right=(int *) malloc(n * sizeof(int));
  int *left_source=(int *) malloc(n * sizeof(int));
  int *right_source=(int *) malloc(n * sizeof(int));
  int left_index = 0, right_index = 0;
  int i, half_dim = dim / 2, nleft = half_dim, nright = dim - half_dim;
  for (i = 0; i < half_dim; i++) {
    left[i] = permutation[i];
    left_source[i] = source[i];
    right[i] = permutation[i + half_dim];
    right_source[i] = source[i + half_dim];
  }
  
  if (nleft < nright) {
    right[i] = permutation[i + half_dim];
    right_source[i] = source[i + half_dim];
  }
  
  for (i = 0; i < dim; i++) {
    if ((left_index < half_dim) && (right_index < dim - half_dim)) {
      if (left[left_index] <= right[right_index]) { // I added "=" in order to support ties
        permutation[i] = left[left_index];
        source[i] = left_source[left_index];
        left_index++;
      } else {
        permutation[i] = right[right_index];
        source[i] = right_source[right_index];
        inversion_count[source[i]] += (half_dim - left_index);
        right_index++;
      }
    } else {
      if (left_index < half_dim) {
        permutation[i] = left[left_index];
        source[i] = left_source[left_index];
        left_index++;
      }
      
      if (right_index < dim - half_dim) {
        permutation[i] = right[right_index];
        source[i] = right_source[right_index];
        right_index++;
      }
    }
  }
  free(left);
  free(right);
  free(left_source);
  free(right_source);
}

int Inversions(int *permutation, int *source, int *inversion_count,int dim, int n)
{
  if (dim==1)
    return 0;
  else{
    Inversions(permutation, source, inversion_count, dim/2, n);
    Inversions(&permutation[dim/2], &source[dim/2], inversion_count,dim-dim/2, n);
    Merge(permutation, source, inversion_count, dim, n);
  }
  return 0;
}

void resample(int *i_perm, int *i_perm_inv, int *n)
{
  int i, j, temp;
  for (i = *n - 1; i > 0; --i) {
    // j = rand() % (i + 1);
    j = r_available_rand() % (i + 1);
    temp = i_perm[j];
    i_perm[j] = i_perm[i];
    i_perm[i] = temp;
  }
  for (i = 0; i < *n; ++i) {
    i_perm_inv[i_perm[i]] = i;
  }
}

double Ball_Information(int *n, double **Dx, double **Dy, int **xidx, int **yidx, int *i_perm, int *i_perm_inv, int *weight)
{
  int i, j, k, pi, src, lastpos, *yrank, *isource, *icount, *xy_index, *xy_temp, **xyidx;
  double pxy, px, py, lastval, rct0=0, *xx_cpy, *yy_cpy;
  yrank = (int *) malloc(*n * sizeof(int));
  isource = (int *) malloc(*n * sizeof(int));
  icount = (int *) malloc(*n * sizeof(int));
  xy_index = (int *) malloc(*n * sizeof(int));
  xy_temp = (int *) malloc(*n * sizeof(int));
  xyidx = alloc_int_matrix(*n, *n);
  xx_cpy = (double *) malloc(*n * sizeof(double));
  yy_cpy = (double *) malloc(*n * sizeof(double));
  
  for(i = 0; i<*n; i++ )
    for(j = 0; j<*n; j++)
      xyidx[i][j] = j;
  
  for(i = 0; i<(*n); i++)
  {
    memcpy(xx_cpy, Dx[i], *n * sizeof(double));
    for(j = 0; j<*n; j++)
      yy_cpy[j] = Dy[i_perm[i]][i_perm[j]];
    quicksort2(xx_cpy, yy_cpy, xyidx[i], 0, *n-1);
  }
  free(xx_cpy);
  free(yy_cpy);
  
  for(i=0;i<(*n);i++)
  {
    pi = i_perm[i];
    lastval = 0;
    lastpos = -1;
    for (j = *n - 1, k = *n - 1; j >= 1; --j, --k)
    {
      k -= (yidx[pi][k] == pi);
      if (lastpos == -1 || Dy[pi][yidx[pi][k]] != lastval)
      {
        lastval = Dy[pi][yidx[pi][k]];
        lastpos = j;
      }
      src = i_perm_inv[yidx[pi][k]];
      src -= (src > i);
      yrank[src] = lastpos;
    }
    
    for (j = 0, k = 0; j < *n - 1; ++j, ++k)
    {
      k += (xyidx[i][k] == i);
      src = xyidx[i][k]; // NOTE: k may be different than from the line above
      src -= (src > i);
      xy_index[j] = yrank[src];
      isource[j] = j;
      icount[j] = 0;
      xy_temp[j] = xy_index[j];
    }
    Inversions(xy_temp, isource, icount, *n - 1, *n);
    lastval = 0;
    lastpos = -1;
    
    for(j = *n - 2, k = *n - 1; j >= 0; --j, --k)
    {
      k -= (xidx[i][k] == i);
      if (lastpos == -1 || Dx[i][xidx[i][k]] != lastval) {
        lastval = Dx[i][xidx[i][k]];
        lastpos = j;
      }
      
      pxy = lastpos - icount[j] + 2;
      px = lastpos + 2;
      py = xy_index[j] + 1;
      px /= (*n);
      py /= (*n);
      pxy /= (*n);
      if(*weight==0)
        rct0 += pow(pxy-px*py,2);
      else
        rct0 += pow(pxy-px*py,2)/(px*py);
    }
    pxy = 0;
    px = 0;
    py = 0;
    for(j=0;j<*n;j++)
    {
      if(Dx[i][xidx[i][j]]==0)
      {
        px+=1;
        if(Dy[pi][i_perm[xidx[i][j]]]==0)
        {
          pxy+=1;
          py+=1;
        }
      }
      else if(Dy[pi][i_perm[xidx[i][j]]]==0)
        py+=1;
    }
    px/=(*n);
    py/=(*n);
    pxy/=(*n);
    if(*weight==0)
      rct0 += pow(pxy-px*py,2);
    else
      rct0 += pow(pxy-px*py,2)/(px*py);
  }
  rct0 = rct0/(1.0*(*n)*(*n));
  free(isource);
  free(icount);
  free(xy_index);
  free(yrank);
  free(xy_temp);
  free_int_matrix(xyidx, *n, *n);
  return(rct0);
}

void BI(double *bcov, double *permuted_bcov, double *x, double *y, int *n, int *R, int *weight)
{
  /*  computes RCT(x,y)  */	
  int    i, j, **xidx, **yidx, *i_perm, *i_perm_inv;
  double **Dx, **Dy, *x_cpy, *y_cpy, RCTV0;
  
  Dx = alloc_matrix(*n, *n);
  Dy = alloc_matrix(*n, *n);
  xidx = alloc_int_matrix(*n, *n);
  yidx = alloc_int_matrix(*n, *n);
  i_perm = (int *) malloc(*n * sizeof(int));
  i_perm_inv = (int *) malloc(*n * sizeof(int));
  x_cpy = (double *) malloc(*n * sizeof(double));
  y_cpy = (double *) malloc(*n * sizeof(double));
  
  vector2matrix(x, Dx, *n, *n, 1);
  vector2matrix(y, Dy, *n, *n, 1);

  for(i = 0; i<*n; i++ )
  {
    for(j = 0; j<*n; j++)
    {
      xidx[i][j] = j;
      yidx[i][j] = j;
    }
    i_perm[i] = i;
    i_perm_inv[i] = i;
  }
  
  for(i=0;i<(*n);i++)
  {
    // copy site to x_cpy and y_cpy
    memcpy(x_cpy, Dx[i], *n * sizeof(double));
    memcpy(y_cpy, Dy[i], *n * sizeof(double));  
    quicksort(x_cpy, xidx[i], 0, *n-1);
    quicksort(y_cpy, yidx[i], 0, *n-1);
  }
  free(x_cpy);
  free(y_cpy);
  
  RCTV0 = Ball_Information( n, Dx, Dy, xidx, yidx, i_perm, i_perm_inv, weight);
  *bcov = RCTV0;
  for(i=0; i<*R; i++)
  {
    if(pending_interrupt()) {
      Rprintf("Process stop due to user interruption! \n");
      break;
    }
    resample(i_perm, i_perm_inv, n);
    permuted_bcov[i] = Ball_Information(n, Dx, Dy, xidx, yidx, i_perm, i_perm_inv, weight);
  }
  
  free_matrix(Dx, *n, *n);
  free_matrix(Dy, *n, *n);
  free_int_matrix(xidx, *n, *n);
  free_int_matrix(yidx, *n, *n);
  free(i_perm);
  free(i_perm_inv);
  return;
}

void computeRank(int n, int **Rank)
{
  int i, j;
  for(i = 1; i < n; i++)
    for(j = 1; j < n; j++)
      Rank[i][j] += (Rank[i][j-1] + Rank[i-1][j] - Rank[i-1][j-1]);
}

void initRank(int n, int **Rank, int *xrank, int *yrank, int *i_perm)
{
  int i, j;
  for(i = 0; i < n + 1; i++)
    for(j = 0; j < n + 1; j++)
      Rank[i][j] = 0;
  for(i = 0; i < n; i++)
    Rank[xrank[i]+1][yrank[i_perm[i]]+1] += 1;
  computeRank(n + 1, Rank);
}

void ranksort(int *n, int *zrank, double *z, int *zidx)
{
  int i, lastpos = 0;
  double lastval = -1.0;
  
  for(i = *n -1; i >= 0; i--){
    if(lastval != z[i])
      lastpos = i;
    lastval = z[i];
    zrank[zidx[i]] = lastpos;
  }
}

void sort(int *n, int *zidx, double *z, int **dzidx)
{
  // the z[i] is the center
  int i, j, zi, ileft, iright, lastpos;
  double lastval, tmp;
  for(i = 0; i < *n; i++){
    zi = zidx[i];
    j = *n - 1;
    lastpos = *n -1;
    ileft = 0;
    iright = *n - 1;
    lastval = -1.0;
    
    while((ileft!=i)||(iright!=i)){
      if(i==ileft){
        tmp = z[iright] - z[i];
        if(lastval != tmp)
          lastpos = j;
        dzidx[zi][zidx[iright]] = lastpos;
        lastval = tmp;
        iright--;
      }
      else if(i==iright){
        tmp = z[i] - z[ileft];
        if(lastval != tmp)
          lastpos = j;
        dzidx[zi][zidx[ileft]] = lastpos;
        lastval = tmp;
        ileft++;
      }
      else{
        if(z[i]-z[ileft]>z[iright]-z[i]){
          tmp = z[i] - z[ileft];
          if(lastval != tmp)
            lastpos = j;
          dzidx[zi][zidx[ileft]] = lastpos;
          lastval = tmp;
          ileft++;
        }
        else{
          tmp = z[iright] - z[i];
          if(lastval != tmp)
            lastpos = j;
          dzidx[zi][zidx[iright]] = lastpos;
          lastval = tmp;
          iright--;
        }
      }
      j--;
    }
    
    if(lastval==0)
      dzidx[zi][zi] = lastpos;
    else
      dzidx[zi][zi] = 0;
  }
}

void createidx(int *n, int *zidx, double *z, int **lowzidx, int **higzidx)
{
  int i, zi, ileft, iright, jleft, jright, lowpos, higpos;
  double lastval, tmp1, tmp2, tmp;
  for(i = 0; i < *n; i++){
    zi = zidx[i];
    lowpos = 1;
    higpos = *n;
    ileft = 0;
    iright = *n - 1;
    jleft = 0;
    jright = 0;
    
    tmp1 = z[iright] - z[i];
    tmp2 = z[i] - z[ileft];
    if(tmp1 > tmp2){
      lastval = tmp1;
      lowzidx[zi][zidx[iright]]=lowpos;
      higzidx[zi][zidx[iright]]=higpos;
      iright--;
      jright++;
    }
    else{
      lastval = tmp2;
      if(ileft == i){
        lowzidx[zi][zidx[iright]]=lowpos;
        higzidx[zi][zidx[iright]]=higpos;
        iright--;
        jright++;
      }
      else{
        lowzidx[zi][zidx[ileft]]=lowpos;
        higzidx[zi][zidx[ileft]]=higpos;
        ileft++;
        jleft++;
      }
    }
    
    while(ileft<=iright){
      tmp1 = z[iright] - z[i];
      tmp2 = z[i] - z[ileft];
      tmp = MAX(tmp1, tmp2);
      while(lastval==tmp){
        if(tmp1 > tmp2){
          lowzidx[zi][zidx[iright]]=lowpos;
          higzidx[zi][zidx[iright]]=higpos;
          iright--;
					jright++;
				}
				else{
					if(ileft == i){
						lowzidx[zi][zidx[iright]]=lowpos;
						higzidx[zi][zidx[iright]]=higpos;
						iright--;
						jright++;
					}
					else{
						lowzidx[zi][zidx[ileft]]=lowpos;
						higzidx[zi][zidx[ileft]]=higpos;
						ileft++;
						jleft++;
					}
				}
				if(iright < ileft)
					break;
				tmp1 = z[iright] - z[i];
				tmp2 = z[i] - z[ileft];
				tmp = MAX(tmp1, tmp2);
			}
			if(iright < ileft)
				break;
			lowpos += jleft;
			higpos -= jright;
			jleft = 0;
			jright = 0;
			if(tmp1 > tmp2){
				lastval = tmp1;
				lowzidx[zi][zidx[iright]]=lowpos;
				higzidx[zi][zidx[iright]]=higpos;
				iright--;
				jright++;
			}
			else{
				lastval = tmp2;
				if(ileft == i){
					lowzidx[zi][zidx[iright]]=lowpos;
					higzidx[zi][zidx[iright]]=higpos;
					iright--;
					jright++;
				}
				else{
					lowzidx[zi][zidx[ileft]]=lowpos;
					higzidx[zi][zidx[ileft]]=higpos;
					ileft++;
					jleft++;
				}
			}
		}
	}
}

void resample2(int *i_perm, int *n)
{
	int i, j, temp;
	for (i = *n - 1; i > 0; --i) {
		// j = rand() % (i + 1);
		j = r_available_rand() % (i + 1);
		temp = i_perm[j];
		i_perm[j] = i_perm[i];
		i_perm[i] = temp;
	}
}

double U_Ball_Information(int *n, int **Rank, int **lowxidx, int **higxidx, int **lowyidx, int **higyidx, int *i_perm, int *weight)
{
	int i, j, pi, pj;
	double px, py, pxy, ans = 0;
	for(i = 0; i < *n; i++){
		for(j = 0; j < *n; j++){
			pi = i_perm[i];
			pj = i_perm[j];
			px = higxidx[i][j] - lowxidx[i][j] + 1;
			py = higyidx[pi][pj] - lowyidx[pi][pj] + 1;
			pxy = Rank[higxidx[i][j]][higyidx[pi][pj]] + Rank[lowxidx[i][j]-1][lowyidx[pi][pj]-1];
			pxy -= (Rank[higxidx[i][j]][lowyidx[pi][pj]-1] + Rank[lowxidx[i][j]-1][higyidx[pi][pj]]);

			pxy /= (*n);
			px /= (*n);
			py /= (*n);
			if(*weight==0)
				ans += pow(pxy-px*py,2);
			else
				ans += pow(pxy-px*py,2)/(px*py);
		}
	}
	ans /= (1.0*(*n)*(*n));
	return(ans);
}


void UBI(double *bcov, double *permuted_bcov, double *x, double *y, int *n, int *R, int *weight) 
{
	int i, j, *xidx, *yidx, *xrank, *yrank, *i_perm, **Rank, **lowxidx, **higxidx, **lowyidx, **higyidx;
	
	xidx = (int *) malloc(*n * sizeof(int));
	yidx = (int *) malloc(*n * sizeof(int));
	xrank = (int *) malloc(*n * sizeof(int));
	yrank = (int *) malloc(*n * sizeof(int));
	i_perm = (int *) malloc(*n * sizeof(int));
	Rank = alloc_int_matrix((*n)+1, (*n)+1);
	lowxidx = alloc_int_matrix(*n, *n);
	higxidx = alloc_int_matrix(*n, *n);
	lowyidx = alloc_int_matrix(*n, *n);
	higyidx = alloc_int_matrix(*n, *n);
	
	for(i = 0; i<*n; i++ )
	{
		xidx[i] = i;
		yidx[i] = i;
		i_perm[i] = i;
	}
	
	// first step: sort the x
	quicksort(x, xidx, 0, *n-1);
	quicksort(y, yidx, 0, *n-1);
	ranksort(n, xrank, x, xidx);
	ranksort(n, yrank, y, yidx);
	createidx(n, xidx, x, lowxidx, higxidx);
	createidx(n, yidx, y, lowyidx, higyidx);
	
	
	initRank(*n, Rank, xrank, yrank, i_perm);
	*bcov = U_Ball_Information(n, Rank, lowxidx, higxidx, lowyidx, higyidx, i_perm, weight);
	for(j=0; j<*R; j++)
	{
	  if(pending_interrupt()) {
	    Rprintf("Process stop due to user interruption! \n");
	    break;
	  }
		resample2(i_perm, n);
		initRank(*n, Rank, xrank, yrank, i_perm);
		permuted_bcov[j] = U_Ball_Information(n, Rank, lowxidx, higxidx, lowyidx, higyidx, i_perm, weight);
	}

	free_int_matrix(Rank, (*n) + 1, (*n) + 1);
	free_int_matrix(lowxidx, *n, *n);
	free_int_matrix(higxidx, *n, *n);
	free_int_matrix(lowyidx, *n, *n);
	free_int_matrix(higyidx, *n, *n);
	free(xidx);
	free(yidx);
	free(xrank);
	free(yrank);
	free(i_perm);
	return;
}


double ubcov_value(double *x, double *y, int *n, int *weight)
{
  int    i, *xidx, *yidx, *xrank, *yrank, *i_perm, **Rank, **lowxidx, **higxidx, **lowyidx, **higyidx;
  double RCTV0;
  
  xidx = (int *) malloc(*n * sizeof(int));
  yidx = (int *) malloc(*n * sizeof(int));
  xrank = (int *) malloc(*n * sizeof(int));
  yrank = (int *) malloc(*n * sizeof(int));
  i_perm = (int *) malloc(*n * sizeof(int));
  Rank = alloc_int_matrix((*n)+1, (*n)+1);
  lowxidx = alloc_int_matrix(*n, *n);
  higxidx = alloc_int_matrix(*n, *n);
  lowyidx = alloc_int_matrix(*n, *n);
  higyidx = alloc_int_matrix(*n, *n);
  
  for(i = 0; i<*n; i++ )
  {
    xidx[i] = i;
    yidx[i] = i;
    i_perm[i] = i;
  }
  
  // first step: sort the x
  quicksort(x, xidx, 0, *n-1);
  quicksort(y, yidx, 0, *n-1);
  ranksort(n, xrank, x, xidx);
  ranksort(n, yrank, y, yidx);
  createidx(n, xidx, x, lowxidx, higxidx);
  createidx(n, yidx, y, lowyidx, higyidx);
  
  
  initRank(*n, Rank, xrank, yrank, i_perm);
  RCTV0 = U_Ball_Information(n, Rank, lowxidx, higxidx, lowyidx, higyidx, i_perm, weight);
  
  free_int_matrix(Rank, (*n) + 1, (*n) + 1);
  free_int_matrix(lowxidx, *n, *n);
  free_int_matrix(higxidx, *n, *n);
  free_int_matrix(lowyidx, *n, *n);
  free_int_matrix(higyidx, *n, *n);
  free(xidx);
  free(yidx);
  free(xrank);
  free(yrank);
  free(i_perm);
  return(RCTV0);
}


double bcov_value(double *x, double *y, int *n, int *weight) 
{
  /*  computes RCT(x,y)  */	
  int    i, j, **xidx, **yidx, *i_perm, *i_perm_inv;
  double **Dx, **Dy, *x_cpy, *y_cpy, RCTV0;
  
  Dx = alloc_matrix(*n, *n);
  Dy = alloc_matrix(*n, *n);
  xidx = alloc_int_matrix(*n, *n);
  yidx = alloc_int_matrix(*n, *n);
  i_perm = (int *) malloc(*n * sizeof(int));
  i_perm_inv = (int *) malloc(*n * sizeof(int));
  x_cpy = (double *) malloc(*n * sizeof(double));
  y_cpy = (double *) malloc(*n * sizeof(double));
  
  vector2matrix(x, Dx, *n, *n, 1);
  vector2matrix(y, Dy, *n, *n, 1);
  
  for(i = 0; i<*n; i++ )
  {
    for(j = 0; j<*n; j++)
    {
      xidx[i][j] = j;
      yidx[i][j] = j;
    }
    i_perm[i] = i;
    i_perm_inv[i] = i;
  }
  
  for(i=0;i<(*n);i++)
  {
    // copy site to x_cpy and y_cpy
    memcpy(x_cpy, Dx[i], *n * sizeof(double));
    memcpy(y_cpy, Dy[i], *n * sizeof(double));  
    quicksort(x_cpy, xidx[i], 0, *n-1);
    quicksort(y_cpy, yidx[i], 0, *n-1);
  }
  free(x_cpy);
  free(y_cpy);
  
  RCTV0 = Ball_Information(n, Dx, Dy, xidx, yidx, i_perm, i_perm_inv, weight);
  
  free_matrix(Dx, *n, *n);
  free_matrix(Dy, *n, *n);
  free_int_matrix(xidx, *n, *n);
  free_int_matrix(yidx, *n, *n);
  free(i_perm);
  free(i_perm_inv);
  return(RCTV0);
}


double bcor_value(double *x, double *y, int *n, int *weight, int *dst)
{
  double bcov_stat_xy, bcov_stat_xx, bcov_stat_yy;
  double bcor_stat;
  if((*dst)) {
    bcov_stat_xy = bcov_value(x, y, n, weight);
    bcov_stat_xx = bcov_value(x, x, n, weight);
    bcov_stat_yy = bcov_value(x, x, n, weight);
  } else {
    bcov_stat_xy = ubcov_value(x, y, n, weight);
    bcov_stat_xx = ubcov_value(x, x, n, weight);
    bcov_stat_yy = ubcov_value(x, y, n, weight);
  }
  bcor_stat = bcov_stat_xy/sqrt(bcov_stat_xx*bcov_stat_yy);
  return(bcor_stat);
}


void UBCOR(double *bcor, double *permuted_bcor, double *x, double *y, int *n, int *R, int *weight)
{
  int i;
  double bcov_stat_xx, bcov_stat_yy;
  UBI(bcor, permuted_bcor, x, y, n, R, weight);
  bcov_stat_xx = ubcov_value(x, x, n, weight);
  bcov_stat_yy = ubcov_value(y, y, n, weight);
  (*bcor) = (*bcor) / bcov_stat_xx / bcov_stat_yy;
  for(i = 0; i < (*R); i++) {
    permuted_bcor[i] = permuted_bcor[i] / bcov_stat_xx / bcov_stat_yy;
  }
}


void BCOR(double *bcor, double *permuted_bcor, double *x, double *y, int *n, int *R, int *weight)
{
  int i;
  double bcov_stat_xx, bcov_stat_yy;
  BI(bcor, permuted_bcor, x, y, n, R, weight);
  bcov_stat_xx = bcov_value(x, x, n, weight);
  bcov_stat_yy = bcov_value(y, y, n, weight);
  (*bcor) = (*bcor) / bcov_stat_xx / bcov_stat_yy;
  for(i = 0; i < (*R); i++) {
    permuted_bcor[i] = permuted_bcor[i] / bcov_stat_xx / bcov_stat_yy;
  }
}


/////////////////////////////////////////////////////////////////
//////////// two functions below wrap key function //////////////
// bcov_stat return ball covariance statistic 
// bcov_test execute ball covariance statistic based test
/////////////////////////////////////////////////////////////////

void bcov_stat(double *bcov, double *x, double *y, int *n, int *weight, int *dst, int *type)
{
  double bcov_stat;
  if((*type) == 1) {
    if((*dst)) {
      bcov_stat = bcov_value(x, y, n, weight);
    } else {
      bcov_stat = ubcov_value(x, y, n, weight);
    }
    *bcov = bcov_stat;
  } else {
    *bcov = bcor_value(x, y, n, weight, dst);
  }
  return;
}


void bcov_test(double *bcov, double *permuted_bcov, double *x, double *y, int *n, int *R, int *weight, int *dst, int *type)
{
  if((*type) == 1) {
    if((*dst)) {
      BI(bcov, permuted_bcov, x, y, n, R, weight);
    } else {
      UBI(bcov, permuted_bcov, x, y, n, R, weight);
    }
  } else {
    if((*dst)) {
      BCOR(bcov, permuted_bcov, x, y, n, R, weight);
    } else {
      UBCOR(bcov, permuted_bcov, x, y, n, R, weight);
    }
  }
  return;
}

