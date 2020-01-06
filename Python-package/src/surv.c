#include "math.h"
#include "stdlib.h"
#include "stdio.h"
#include "surv.h"
#include "utilities.h"

/**
 * @param x : covariate rank
 * @param t : time rank
 * @param delta : censored index
 * @param Sc : KM estimator for censor time
 * @param n : sample size
 * @param RCTV : BCor survival extension
 */
void SRCT_new(double *x, int *t, int *delta, double *Sc, int *n, double *RCTV) {
    int i, j, *xidx, *xrank;
    double *xcopy;
    double inv_n = 1.0 / *n, jp = 0.0, p1, p2;
    *RCTV = 0.0;

    xidx = (int *) malloc(*n * sizeof(int));
    xrank = (int *) malloc(*n * sizeof(int));
    xcopy = (double *) malloc(*n * sizeof(double));

    for (i = 0; i < (*n); i++) {
        xidx[i] = i;
        xcopy[i] = x[i];
    }

    quicksort(xcopy, xidx, 0, *n - 1);
    ranksort(n, xrank, xcopy, xidx);
    free(xidx);
    free(xcopy);

    for (i = 0; i < (*n); i++) {
        if (delta[i] == 1) {
            p1 = *n - xrank[i] - 1;
            p2 = *n - t[i] - 1;
            for (j = i + 1; j < (*n); j++) {
                if ((x[j] > x[i]) & (t[j] > t[i])) {
                    jp += 1.0;
                }
            }
            *RCTV += pow(jp * inv_n - p1 * p2 * inv_n * inv_n, 2) / pow(Sc[i], 3);
        }
        jp = 0.0;
    }
    *RCTV = *RCTV * inv_n;

    free(xrank);
}