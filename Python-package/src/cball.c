#include <stdlib.h>
#include <stdio.h>
#include "BD.h"
#include "cball.h"

double bd_stat(double *xy, int n1, int n2, int weight, int dst)
{
	double ans = 0.0;
	if ((dst)) {
		ans = bd_value(xy, &n1, &n2, &weight);
	}
	else {
		ans = ubd_value(xy, &n1, &n2, &weight);
	}

	return ans;
}


void bd_test(double *bd, double *permuted_bd, double *xy, int n1, int n2, int p, int dst, int R, int weight)
{
	if (dst) {
		BD(bd, permuted_bd, xy, &n1, &n2, &p, &dst, &R, &weight);
	}
	else {
		UBD(bd, permuted_bd, xy, &n1, &n2, &R, &weight);
	}

	/*printf("%f \n",*bd);
	printf("%f \n",permuted_bd[0]);
	printf("%f \n", permuted_bd[1]);
	printf("%f \n", permuted_bd[2]);
	printf("%f \n", permuted_bd[3]);
	printf("%f \n", permuted_bd[4]);*/

	return;
}