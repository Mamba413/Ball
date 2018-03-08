#include <stdlib.h>
#include <stdio.h>
#include "BI.h"
#include "BD.h"

void test_value(double *computed_value, double true_value)
{
	const double EPSINON = 0.0000001;
	double difference = (*computed_value) - true_value;
	if ((difference < EPSINON) & (difference > -EPSINON)) {
		printf("Success!\n\n");
	}
	else {
		printf("Fail!\n\n");
	}
	return;
}

void main()
{
	double tmp = 0.0;
	double *bcov_stat_value = &tmp, *bd_stat_value = &tmp;
	double x[10] = { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10 };
	double y[10] = { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10 };
	double z[15] = { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15 };
	int n = 10, weight = 0, dst = 0, type = 1, K = 2;
	int size[2] = { 5, 5 };
	int nth = 1;
	printf("----------- Test Computation for ball divergence (Two sample) -----------\n");
	bd_stat(bd_stat_value, x, size, &n, &K, &weight, &dst, &nth);
	printf("Computation result: %f; ", *bd_stat_value);
	printf("True value: %f\nTest result: ", 0.7232);
	test_value(bd_stat_value, 0.7232);
	printf("----------- Test Computation for ball divergence (K sample) -----------\n");
	int size1[3] = { 5, 5, 5 }; 
	int n1 = 15;
	K = 3;
	bd_stat(bd_stat_value, z, size1, &n1, &K, &weight, &dst, &nth);
	printf("Computation result: %f; ", *bd_stat_value);
	printf("True value: %f\nTest result: ", 2.4032);
	test_value(bd_stat_value, 2.4032);
	printf("----------- Test Computation for ball covariance -----------\n");
	bcov_stat(bcov_stat_value, x, y, &n, &weight, &dst, &type);
	printf("Computation result: %f; ", *bcov_stat_value);
	printf("True value: %f\nTest result: ", 0.034214);
	test_value(bcov_stat_value, 0.034214);
	printf("----------- Test Computation for ball correlation -----------\n");
	type = 2;
	bcov_stat(bcov_stat_value, x, y, &n, &weight, &dst, &type);
	printf("Computation result: %f; ", *bcov_stat_value);
	printf("True value: %f\nTest result: ", 1.0000);
	test_value(bcov_stat_value, 1.0000);
	system("pause");
	return;
}