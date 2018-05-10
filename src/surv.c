#include "math.h"
#include "stdlib.h"
#include "stdio.h"
#include "surv.h"
#include "utilities.h"

// modify
/*
Input
x (covariate) : { 0.40542, 2.14187, -0.71706, 0.73307, -1.09226, -1.86642, 0.55529, -0.61887, -0.49071, -1.20182, -0.70259, 0.61007, -0.63128, 0.41266, 0.8385, -0.43179, -0.74084, -1.07779, 1.41215, 1.3739, 0.3186, -0.592, 0.61007, 1.52067, -0.70776, -0.50311, 1.869, -0.70569, 0.68139, 1.95582 };
t (sorted time) : { 0, 0.2, 0.2, 0.3, 0.3, 0.4, 0.5, 0.5, 0.5, 0.6, 0.6, 0.6, 0.7, 0.7, 0.7, 0.8, 0.8, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 1, 1, 1, 1, 1.1, 1.1, 1.1 };
delta (censored index) : { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1 };
Sc (KM estimator for censor time) : {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0.92308, 0.92308, 0.92308, 0.92308, 0.92308, 0.92308, 0.92308, 0.92308, 0.92308, 0.92308, 0.92308, 0.92308, 0.92308};
n (sample size) : 30
RCTV (BCor survival extension) : 0.001264425
Computation complexity: O(n^2)
*/
void SRCT(double *x, double *t, int *delta, double *Sc, int *n, double *RCTV)
{
	/*  computes RCT(x,y)  */
  int    i, j;
  double jp = 0.0, p1 = 0.0, p2 = 0.0;
  *RCTV = 0;

  //printf("------- Oringinal Algorithm --------\n");
  for(i = 0; i < (*n); i++){
	  if (delta[i] == 1) {
		for(j=0;j<(*n);j++){
            if((x[j]>x[i]) && (t[j]>t[i]))
			    jp += 1.0;
            if(x[j]>x[i])
                p1 += 1.0;
            if(t[j]>t[i])
                p2 += 1.0;
		}
		//printf("sample-%d --- p1: %f, p2: %f, jp: %f\n", i, p1, p2, jp);
		*RCTV += pow(jp / (*n) - p1*p2 / ((*n)*(*n)), 2) / pow(Sc[i], 3);
	  }
	  jp = 0.0;
	  p1 = 0.0;
	  p2 = 0.0;
	}
  *RCTV = *RCTV/(1.0*(*n));
  return;
}


/*
Input:
x (covariate rank) : { 16, 29, 5, 22, 2, 0, 18, 10, 13, 1, 8, 20, 9, 17, 23, 14, 4, 3, 25, 24, 15, 11, 20, 26, 6, 12, 27, 7, 21, 28 };
t (time rank) : { 0, 2, 2, 4, 4, 5, 8, 8, 8, 11, 11, 11, 14, 14, 14, 16, 16, 22, 22, 22, 22, 22, 22, 26, 26, 26, 26, 29, 29, 29 };
delta (censored index) : { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1 };
Sc (KM estimator for censor time) : {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0.92308, 0.92308, 0.92308, 0.92308, 0.92308, 0.92308, 0.92308, 0.92308, 0.92308, 0.92308, 0.92308, 0.92308, 0.92308};
n (sample size) : 30
RCTV (BCor survival extension) : 0.001264425

Test code:
int surv_x_rank[30] = { 16, 29, 5, 22, 2, 0, 18, 10, 13, 1, 8, 20, 9, 17, 23, 14, 4, 3, 25, 24, 15, 11, 20, 26, 6, 12, 27, 7, 21, 28 };
int surv_t_rank[30] = { 0, 2, 2, 4, 4, 5, 8, 8, 8, 11, 11, 11, 14, 14, 14, 16, 16, 22, 22, 22, 22, 22, 22, 26, 26, 26, 26, 29, 29, 29 };
SRCT_fast(surv_x_rank, surv_t_rank, surv_delta, surv_sc, &surv_num, surv_bcor);
printf("Computation result: %f; ", *surv_bcor);
printf("True value: %f\nTest result: ", 0.001264425);
test_value(surv_bcor, 0.001264425);

*/
void SRCT_new(int *x, int *t, int *delta, double *Sc, int *n, double *RCTV)
{
	int i, j;
	double jp = 0.0, p1 = 0.0, p2 = 0.0;
	*RCTV = 0.0;

	//printf("------- Speedup Algorithm --------\n");
	for (i = 0; i < (*n); i++) {
		if (delta[i] == 1) {
			p1 = *n - x[i] - 1;
			p2 = *n - t[i] - 1;
			for (j = i + 1; j<(*n); j++) {
				if ( (x[j] > x[i]) & (t[j] > t[i]) )
					jp += 1.0;
			}
			*RCTV += pow(jp / (*n) - p1*p2 / ((*n)*(*n)), 2) / pow(Sc[i], 3);
		}
		jp = 0.0;
		p1 = 0.0;
		p2 = 0.0;
	}
	*RCTV = *RCTV / (1.0*(*n));
	return;
}



/*
Input:
x (covariate rank) : { 16, 29, 5, 22, 2, 0, 18, 10, 13, 1, 8, 20, 9, 17, 23, 14, 4, 3, 25, 24, 15, 11, 20, 26, 6, 12, 27, 7, 21, 28 };
x_sorted_rank: { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 20, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29 };
t (time rank) : { 0, 2, 2, 4, 4, 5, 8, 8, 8, 11, 11, 11, 14, 14, 14, 16, 16, 22, 22, 22, 22, 22, 22, 26, 26, 26, 26, 29, 29, 29 };
t_rank2 : { 0, 1, 1, 3, 3, 5, 6, 6, 6, 9, 9, 9, 12, 12, 12, 15, 15, 17, 17, 17, 17, 17, 17, 23, 23, 23, 23, 27, 27, 27 };
delta (censored index) : { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1 };
Sc (KM estimator for censor time) : {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0.92308, 0.92308, 0.92308, 0.92308, 0.92308, 0.92308, 0.92308, 0.92308, 0.92308, 0.92308, 0.92308, 0.92308, 0.92308};
n (sample size) : 30
RCTV (BCor survival extension) : 0.001264425

Test code:
int surv_x_rank[30] = { 16, 29, 5, 22, 2, 0, 18, 10, 13, 1, 8, 20, 9, 17, 23, 14, 4, 3, 25, 24, 15, 11, 20, 26, 6, 12, 27, 7, 21, 28 };
int surv_x_sorted_rank[30] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 20, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29 };
int surv_t_rank[30] = { 0, 2, 2, 4, 4, 5, 8, 8, 8, 11, 11, 11, 14, 14, 14, 16, 16, 22, 22, 22, 22, 22, 22, 26, 26, 26, 26, 29, 29, 29 };
int surv_t_rank2[30] = { 0, 1, 1, 3, 3, 5, 6, 6, 6, 9, 9, 9, 12, 12, 12, 15, 15, 17, 17, 17, 17, 17, 17, 23, 23, 23, 23, 27, 27, 27 };
SRCT_fast(surv_x_rank, surv_x_sorted_rank, surv_t_rank, surv_t_rank2, surv_delta, surv_sc, &surv_num, surv_bcor);
printf("Computation result: %f; ", *surv_bcor);
printf("True value: %f\nTest result: ", 0.001264425);
test_value(surv_bcor, 0.001264425);

*/
//void SRCT_fast(int *x, int * x_sorted_rank, int *t, int *t_rank_min, int *delta, double *Sc, int *n, double *RCTV)
//{
//	int i, *i_perm, **Rank;
//	double jp = 0.0, p1 = 0.0, p2 = 0.0;
//	*RCTV = 0;
//
//	i_perm = (int *) malloc(*n * sizeof(int));
//	Rank = alloc_int_matrix((*n) + 1, (*n) + 1);
//
//	for (i = 0; i<*n; i++)
//	{
//		i_perm[i] = i;
//	}
//	initRank_surv(*n, Rank, x_sorted_rank, t_rank_min, i_perm);
//	for (i = *n-1; i >= 0; i--)
//	{
//		for (int j = 0; j < (*n); j++)
//			printf("%d ", Rank[i][j]);
//		printf("\n");
//	}
//
//	printf("------- Speedup Algorithm --------\n");
//	for (i = 0; i < (*n); i++) {
//		if (delta[i] == 1) {
//			p1 = *n - x[i] - 1;
//			p2 = *n - t[i] - 1;
//			jp = *n + Rank[t[i]][x[i]] - 1;
//			jp -= (Rank[t[i]][*n] + Rank[*n][x[i]]);
//			printf("sample-%d --- p1: %f, p2: %f, jp: %f\n", t[i], p1, p2, jp);
//			// 
//			*RCTV += pow(jp / (*n) - p1*p2 / ((*n)*(*n)), 2) / pow(Sc[i], 3);
//		}
//		jp = 0.0;
//		p1 = 0.0;
//		p2 = 0.0;
//	}
//	*RCTV = *RCTV / (1.0*(*n));
//	return;
//}
