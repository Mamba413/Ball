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


double pvalue(double ball_stat_value, double *permuted_ball, int R)
{
	double larger_num = 0.0;
	for (int i = 0; i < R; i++)
	{
		if (permuted_ball[i] < ball_stat_value)
		{
			larger_num += 1.0;
		}
	}
	double R_double = R;
	double p_value = (1.0 + larger_num) / (1.0 + R_double);
	return(p_value);
}

void test_multithread_permutation()
{
	double ball_stat_value = 0.0;
	double permuted_ball[50000];
	//double x[500] = { -0.28286, 0.38531, 2.22329, -0.76211, 1.40028, -0.55494, 0.35822, -1.37954, -0.21038, 0.27986, -0.60309, -0.70182, 0.93769, -0.01463, -1.55417, 2.39335, 0.12592, -0.87153, -1.16284, -1.28369, -0.4243, 0.31727, 0.55048, -0.03054, 0.96479, 1.09353, 0.41857, -0.24966, 0.14913, 0.86621, 0.70147, -0.15248, -1.01218, 0.46842, 1.17143, -0.05809, 1.03764, 0.57203, -1.25455, 0.69609, -0.43046, 0.12532, 1.09387, 1.4611, -0.56635, -0.94728, -2.26491, -0.3478, 0.86459, 0.52679, -0.6007, -2.35246, 0.50167, 0.36458, -0.30171, -0.48328, -0.16348, -2.07051, 2.77379, -2.13392, 1.74078, -1.337, 0.69181, 1.63559, -0.04504, 0.01383, -0.91121, -0.25594, -1.19633, -2.73983, 0.12689, -0.7549, 1.39338, 1.39374, 0.16725, -1.53045, 0.4178, 0.40497, -2.577, -0.68729, 0.2925, 1.01131, -0.25546, -0.82506, 0.51708, 0.94347, -1.21907, -0.28145, -0.20382, 1.04423, -0.80029, 0.16019, 1.35386, -0.24158, -0.46558, 1.29068, 1.09765, 0.26908, -0.21747, -0.15166, -1.5853, 0.25974, 1.3792, 0.74593, -0.2694, 0.2522, 0.99806, 1.87947, -0.16804, -1.07538, 0.29829, -0.51545, 0.67964, 0.03924, -0.53052, 2.17787, 1.38056, 1.35234, -0.06366, -0.03202, -1.46319, -0.91758, 0.48526, 0.50822, 0.90617, 1.50162, -0.62176, 0.94162, 0.89386, 0.86655, -1.67444, -0.19613, 0.07439, 1.31515, 0.8651, -1.13264, 0.48882, -0.86943, -0.15036, 0.74736, 0.35769, -1.25198, -0.58194, 0.53914, 0.60833, 0.16008, -0.27818, 0.13994, -0.87743, -0.56001, 0.34178, 0.79299, 1.38646, 0.37407, 0.89687, -0.10525, 0.6059, -1.32785, 0.64818, 1.27565, 0.1485, 0.9363, -0.69653, 0.47559, -1.27983, -0.11413, 0.19235, -1.59595, 0.45205, 0.63413, -0.20855, -1.76373, 0.05714, -0.21221, 1.33361, 1.07148, 0.78614, -0.03844, -0.65694, -1.5841, -0.39861, 0.84106, 1.27957, -0.71609, 0.69837, 0.81938, 0.71218, -0.30424, 0.41076, -0.68297, -0.67915, -0.34759, 0.20584, -0.76885, -0.95456, -1.04695, 0.10225, -0.35122, 0.74137, -0.74975, -2.59687, -0.35466, -0.53445, 0.42892, -2.01333, 0.75063, 1.0034, 0.23518, -1.2277, -0.59236, -0.18962, 0.03932, 2.12523, -0.10181, 0.25381, 0.15975, 0.42216, -1.81404, 0.95111, 0.55219, 0.56889, 0.2287, -0.19493, 0.2728, 0.94426, 0.21065, -0.7898, 0.78541, -0.24545, -0.47833, 0.5362, -0.43754, -2.28452, -0.08657, 1.0411, 1.15262, 1.53016, -0.39607, 0.3459, -0.94028, -0.12632, -0.42809, 0.24935, -0.3393, -1.85463, -0.30735, 0.25021, 1.3685, 0.20503, 1.27274, 0.77683, 0.57041, 0.25935, -2.70602, -1.96288, 0.33073, 2.63683, 0.69188, -2.42177, 1.13245, -0.05395, 0.87197, -0.82445, -0.13079, 0.67253, 0.42959, 0.03598, -0.80142, 1.7489, 0.40219, 0.95098, -0.6845, 0.14139, 0.30405, -1.32127, -0.43797, 0.23635, 0.98814, -1.96802, -0.80206, 0.70918, -0.12896, -0.08991, 0.35476, -1.34868, -0.29685, -1.11314, 0.25082, 0.9235, 0.15557, 0.3918, 1.11155, 0.32854, -0.26187, -0.77111, 1.3424, 0.77516, 1.07405, -0.78632, 0.19546, -0.62053, 0.38465, -1.15366, 0.16446, -0.40097, 1.56468, -0.08712, 1.90791, -2.07696, -0.18683, -0.83539, -0.28957, 1.87706, 0.23155, -1.58951, -0.01193, -0.13118, -1.05663, -1.13767, 0.54567, 0.50842, -0.77564, -1.69976, 1.55247, 0.52171, -0.06777, -0.15177, -0.12102, 0.67503, 0.91593, 0.77352, -1.28832, 0.29937, 0.66758, 0.14692, 0.47003, 1.15367, -0.91624, 1.73044, 0.79869, 0.27583, -0.62377, -0.53435, -1.02462, -1.11374, 0.95712, -0.04835, -1.61926, 1.06452, 0.58052, -0.43699, -0.18544, 1.15503, -1.21649, -0.49999, -0.48779, 1.2659, -0.91662, 0.87498, -0.32051, -0.48184, -0.51414, -0.60025, -1.56989, -0.26628, -0.52593, -1.37707, 0.42175, -0.06456, -0.11621, 0.82846, -0.22175, 0.38681, -0.92702, 0.59844, 0.35485, 0.82865, 1.70262, 0.85593, -1.23114, -0.46816, 1.15169, 1.17061, 2.64063, -1.43942, -1.81606, -0.63889, -0.26714, 0.31959, -1.09909, -0.43546, -1.6659, -0.20418, -0.90284, -0.55799, 0.65804, -2.21253, 0.03993, 1.56496, 1.11105, 0.06034, 0.18298, -1.02883, 1.445, 0.44101, -2.06828, -2.52204, 1.2123, -0.41204, 0.10966, -1.94188, 1.55886, -0.81227, 2.77714, 0.64138, 1.81328, 0.9062, -0.02516, -1.23612, 0.4609, 0.27104, 0.29702, -0.32621, -0.43322, -0.31227, -1.14037, 0.30055, -0.92976, 1.3539, -0.87522, 1.47265, 0.88548, -1.46541, 0.63836, -0.54284, -1.07616, -0.74759, -2.04795, 1.25068, -0.1889, -0.45503, -0.0764, 0.86397, -0.26969, 0.95822, 0.35275, -0.53595, -1.06847, 1.05095, -0.41824, -0.05256, 0.12689, 0.5121, 1.97475, 0.35614, -1.34633, 2.03104, -0.38114, 0.90931, 0.71736, -1.85646, 1.12352, 0.28742, -1.37027, 0.75223, 0.34743, -0.05956, -0.96693, 0.08652, 1.15613, -0.6594, 2.18718, -0.48262, 1.16866, -0.35341, -0.70089, 0.05677, -0.43704, 0.07423, -0.79479, -0.62094, 0.206, 0.43173, -0.31525, -1.82816, -0.37565, -1.13833, 0.15032, -0.16239, 1.03816, -0.34908, 0.68789, 0.312, 2.05022, -0.03686, -0.39797, 0.08494, -0.75486, -1.8527, 0.08957 };
	double x[15] = { 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0 };
	//int size[2] = { 250, 250 };
	//int k = 2;
	int size[5] = { 5, 5, 5 };
	int k = 3;
	int n = 15, dst = 0, R = 200, weight = 0, nth = 4;
	int parallel_type = 2;
	bd_test(&ball_stat_value, permuted_ball, x, size, &n, &k, &dst, &R, &weight, &nth, &parallel_type);
	double p_value = 0.0;
	p_value = pvalue(ball_stat_value, permuted_ball, R);
	printf("Ball statistics: %f; ", ball_stat_value);
	printf("p-value: %f \n", p_value);
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
	printf("----------- Test Computation for ball divergence (Two sample, Univariate) -----------\n");
	bd_stat(bd_stat_value, x, size, &n, &K, &weight, &dst, &nth);
	printf("Computation result: %f; ", *bd_stat_value);
	printf("True value: %f\nTest result: ", 0.7232);
	test_value(bd_stat_value, 0.7232);
	printf("----------- Test Computation for ball divergence (Two sample, Multivariate) -----------\n");
	double x_dst[100] = { 0, 1.41421, 2.82843, 4.24264, 5.65685, 7.07107, 8.48528, 9.89949, 11.31371, 12.72792, 1.41421, 0, 1.41421, 2.82843, 4.24264, 5.65685, 7.07107, 8.48528, 9.89949, 11.31371, 2.82843, 1.41421, 0, 1.41421, 2.82843, 4.24264, 5.65685, 7.07107, 8.48528, 9.89949, 4.24264, 2.82843, 1.41421, 0, 1.41421, 2.82843, 4.24264, 5.65685, 7.07107, 8.48528, 5.65685, 4.24264, 2.82843, 1.41421, 0, 1.41421, 2.82843, 4.24264, 5.65685, 7.07107, 7.07107, 5.65685, 4.24264, 2.82843, 1.41421, 0, 1.41421, 2.82843, 4.24264, 5.65685, 8.48528, 7.07107, 5.65685, 4.24264, 2.82843, 1.41421, 0, 1.41421, 2.82843, 4.24264, 9.89949, 8.48528, 7.07107, 5.65685, 4.24264, 2.82843, 1.41421, 0, 1.41421, 2.82843, 11.31371, 9.89949, 8.48528, 7.07107, 5.65685, 4.24264, 2.82843, 1.41421, 0, 1.41421, 12.72792, 11.31371, 9.89949, 8.48528, 7.07107, 5.65685, 4.24264, 2.82843, 1.41421, 0 };
	dst = 1;
	bd_stat(bd_stat_value, x_dst, size, &n, &K, &weight, &dst, &nth);
	printf("Computation result: %f; ", *bd_stat_value);
	printf("True value: %f\nTest result: ", 0.7232);
	test_value(bd_stat_value, 0.7232);
	printf("----------- Test Computation for ball divergence (K sample) -----------\n");
	int size1[3] = { 5, 5, 5 }; 
	int n1 = 15;
	K = 3;
	dst = 0;
	bd_stat(bd_stat_value, z, size1, &n1, &K, &weight, &dst, &nth);
	printf("Computation result: %f; ", *bd_stat_value);
	printf("True value: %f\nTest result: ", 2.4032);
	test_value(bd_stat_value, 2.4032);
	printf("----------- Test Computation for ball covariance (Univariate) -----------\n");
	bcov_stat(bcov_stat_value, x, y, &n, &weight, &dst, &type, &nth);
	printf("Computation result: %f; ", *bcov_stat_value);
	printf("True value: %f\nTest result: ", 0.034214);
	test_value(bcov_stat_value, 0.034214);
	printf("----------- Test Computation for ball covariance (Multivariate, Case 1) -----------\n");
	double y_dst[100] = { 0, 1.41421, 2.82843, 4.24264, 5.65685, 7.07107, 8.48528, 9.89949, 11.31371, 12.72792, 1.41421, 0, 1.41421, 2.82843, 4.24264, 5.65685, 7.07107, 8.48528, 9.89949, 11.31371, 2.82843, 1.41421, 0, 1.41421, 2.82843, 4.24264, 5.65685, 7.07107, 8.48528, 9.89949, 4.24264, 2.82843, 1.41421, 0, 1.41421, 2.82843, 4.24264, 5.65685, 7.07107, 8.48528, 5.65685, 4.24264, 2.82843, 1.41421, 0, 1.41421, 2.82843, 4.24264, 5.65685, 7.07107, 7.07107, 5.65685, 4.24264, 2.82843, 1.41421, 0, 1.41421, 2.82843, 4.24264, 5.65685, 8.48528, 7.07107, 5.65685, 4.24264, 2.82843, 1.41421, 0, 1.41421, 2.82843, 4.24264, 9.89949, 8.48528, 7.07107, 5.65685, 4.24264, 2.82843, 1.41421, 0, 1.41421, 2.82843, 11.31371, 9.89949, 8.48528, 7.07107, 5.65685, 4.24264, 2.82843, 1.41421, 0, 1.41421, 12.72792, 11.31371, 9.89949, 8.48528, 7.07107, 5.65685, 4.24264, 2.82843, 1.41421, 0 };
	dst = 1;
	bcov_stat(bcov_stat_value, x_dst, y_dst, &n, &weight, &dst, &type, &nth);
	printf("Computation result: %f; ", *bcov_stat_value);
	printf("True value: %f\nTest result: ", 0.034214);
	test_value(bcov_stat_value, 0.034214);
	printf("----------- Test Computation for ball covariance (Multivariate, Case 2) -----------\n");
	double x_dst1[100] = { 0, 2.94102, 2.17977, 1.75621, 4.17487, 2.2276, 2.75585, 0.63125, 4.43155, 5.0511, 2.94102, 0, 4.84657, 2.26971, 1.53009, 1.11437, 0.24297, 2.79024, 4.16784, 2.11475, 2.17977, 4.84657, 0, 2.83416, 6.23703, 4.34783, 4.62333, 2.07098, 4.4642, 6.93359, 1.75621, 2.26971, 2.83416, 0, 3.79152, 2.30828, 2.02694, 1.18126, 2.77045, 4.21881, 4.17487, 1.53009, 6.23703, 3.79152, 0, 1.9614, 1.77209, 4.1671, 5.54001, 1.34456, 2.2276, 1.11437, 4.34783, 2.30828, 1.9614, 0, 1.09868, 2.30866, 4.76509, 2.99578, 2.75585, 0.24297, 4.62333, 2.02694, 1.77209, 1.09868, 0, 2.57511, 3.97415, 2.32058, 0.63125, 2.79024, 2.07098, 1.18126, 4.1671, 2.30866, 2.57511, 0, 3.80297, 4.89524, 4.43155, 4.16784, 4.4642, 2.77045, 5.54001, 4.76509, 3.97415, 3.80297, 0, 5.28156, 5.0511, 2.11475, 6.93359, 4.21881, 1.34456, 2.99578, 2.32058, 4.89524, 5.28156, 0 };
	double y_dst1[100] = { 0, 0.43717, 0.38432, 0.58089, 0.8515, 0.9343, 0.82414, 0.87082, 0.82319, 0.61563, 0.43717, 0, 0.82149, 1.01806, 1.28867, 0.49713, 0.38697, 1.30799, 1.26036, 1.0528, 0.38432, 0.82149, 0, 0.19657, 0.46718, 1.31862, 1.20846, 0.48651, 0.43887, 0.23132, 0.58089, 1.01806, 0.19657, 0, 0.27061, 1.51519, 1.40503, 0.28994, 0.2423, 0.03475, 0.8515, 1.28867, 0.46718, 0.27061, 0, 1.7858, 1.67564, 0.01933, 0.02831, 0.23586, 0.9343, 0.49713, 1.31862, 1.51519, 1.7858, 0, 0.11016, 1.80512, 1.75749, 1.54993, 0.82414, 0.38697, 1.20846, 1.40503, 1.67564, 0.11016, 0, 1.69497, 1.64733, 1.43977, 0.87082, 1.30799, 0.48651, 0.28994, 0.01933, 1.80512, 1.69497, 0, 0.04763, 0.25519, 0.82319, 1.26036, 0.43887, 0.2423, 0.02831, 1.75749, 1.64733, 0.04763, 0, 0.20756, 0.61563, 1.0528, 0.23132, 0.03475, 0.23586, 1.54993, 1.43977, 0.25519, 0.20756, 0 };
	dst = 1;
	nth = 2;
	bcov_stat(bcov_stat_value, x_dst1, y_dst1, &n, &weight, &dst, &type, &nth);
	printf("Computation result: %f; ", *bcov_stat_value);
	printf("\n\n");
	//printf("True value: %f\nTest result: ", 0.034214);
	//test_value(bcov_stat_value, 0.034214);
	printf("----------- Test Computation for ball correlation -----------\n");
	dst = 0;
	type = 2;
	bcov_stat(bcov_stat_value, x, y, &n, &weight, &dst, &type, &nth);
	printf("Computation result: %f; ", *bcov_stat_value);
	printf("True value: %f\nTest result: ", 1.0000);
	test_value(bcov_stat_value, 1.0000);

	//printf("----------- Test Multi-thread -----------\n");
	//double permuted_bcov[10000];
	//double m1[500] = { 1.08837, -0.97366, -0.14463, -1.74737, -0.33314, -0.59724, -0.43556, 0.90097, -0.17926, -0.50293, -0.46313, 1.34861, 0.69907, -0.14161, -0.01974, -0.28245, -0.58671, -0.19274, -1.19092, -1.0305, 0.21805, 2.11794, -0.52646, 1.8246, 0.81556, -0.36491, 0.71661, -0.83224, 0.21435, 2.12217, 0.26824, 0.54749, -0.90835, 0.63278, -1.89402, -0.78611, 0.65132, 2.34997, -0.20414, -0.7907, 0.00118, 0.72869, -0.62324, -0.58815, 1.39272, -0.01284, -1.16032, -0.59133, -0.50471, -0.66761, 0.74484, -0.63312, 0.70599, 0.37722, -1.96205, 2.68086, -0.56146, -0.09524, -0.82673, 0.90261, 0.59155, -0.76322, -0.34845, -0.40582, 1.63254, 0.10744, 0.97182, -1.60693, 0.34145, -0.44532, 0.06644, -1.14121, 1.38032, 0.18508, 1.37783, 0.34288, 0.37413, 0.02954, 1.65979, -0.63205, 0.39136, 1.16937, 1.4187, 0.33232, 0.09469, 0.18154, -0.57579, 0.05015, -0.76944, -0.06505, 0.65398, 1.52831, -0.78746, 0.77654, 0.37783, 1.38389, -0.03591, -0.48486, 0.64004, -1.31675, 0.36724, -0.055, -0.98871, -1.05039, 1.02977, -0.12941, 0.90065, -0.31512, -0.39973, -0.5927, 2.21613, 3.27482, -1.06057, 1.1904, -0.20048, -1.78073, 0.10895, -0.34348, -1.06109, -0.03309, -0.64515, 0.40684, -0.11016, -0.98614, 1.16335, -1.11748, 0.45911, -0.48384, -1.78631, 2.05555, 0.80558, 0.22366, -0.12919, -0.40823, -1.70383, 0.10302, -0.07558, -0.49967, 0.15184, -0.09365, -0.6384, -0.10343, -1.31835, -0.43143, -1.0503, -0.53429, 1.24063, 1.68421, 1.19785, -1.03705, 0.97933, -0.08066, -1.40248, 0.47796, -0.14243, 0.57394, 0.41081, 0.33121, -1.19774, -1.51809, -0.11089, -1.6925, 0.8778, -0.26082, 1.74338, -1.12757, 1.06057, -0.29711, 1.06931, 1.60571, -1.08619, -0.51822, 0.75872, -0.57019, -0.21889, -0.51164, 2.18044, 1.12248, 0.43728, 0.15939, -1.22468, 0.97876, -0.93903, -0.58198, -0.06572, -1.39376, -2.57313, -0.31769, -0.92293, -0.3034, -0.26359, -0.05068, 0.59819, -0.5007, -0.79025, -0.31031, 0.01873, -1.03714, -1.43518, 0.44748, 1.31825, -2.54427, 0.29688, 0.69988, 2.45557, 1.46687, -0.17652, -0.1661, 1.12561, -0.99097, 1.74534, -0.34647, 0.31116, 0.51422, -0.01489, 0.60155, 2.30481, -0.12264, -0.97583, -0.50407, -1.98115, 0.30902, 0.73245, 1.99584, 2.04332, -0.61705, 0.04842, 1.36524, 1.18251, 2.94207, -1.9307, 0.79862, 1.63472, 0.82461, -1.99281, -0.37531, -0.0894, -0.20632, -0.65629, -0.55689, 0.26395, 0.96604, -0.2776, -0.03345, 0.69473, 0.10187, 0.5109, -0.87169, -1.93284, -2.01018, -1.91614, 1.24273, 1.58289, -0.43579, -0.3023, -0.26309, -1.37543, -0.71144, -1.10916, -1.28979, -1.39923, 1.58925, -1.37161, -0.53388, 0.34251, -1.50549, 1.55516, -0.25351, 1.76456, 1.49681, 0.27167, -0.06571, -0.73207, -1.465, 1.68213, -0.7596, 0.65023, -0.46029, 0.07222, -0.59627, 1.00329, 0.8737, 1.03012, -0.67591, 0.24637, 0.01737, 0.28878, -0.41112, 0.73276, -0.26209, -1.80494, -0.47782, 0.19314, 0.38121, -0.52558, -1.13714, -0.31697, 1.46324, 0.8065, -0.57543, -0.33364, -1.3644, 1.18588, 0.83886, 1.00207, 0.33923, -2.64925, -0.53554, 1.09432, 0.1108, -0.22685, -0.19691, -0.0756, 0.11181, -0.67366, -0.12541, 0.67958, 1.51639, 1.08561, -0.47248, -0.25048, 0.13187, 1.92332, -0.24112, 0.16005, -0.31166, -0.27018, 0.83693, 1.16999, -1.07348, 0.65459, 1.06001, -0.78706, -0.09592, -1.48394, -0.29017, 0.37091, -1.67488, 0.62966, -0.18485, -0.19833, -0.80221, -0.05065, -1.4575, 0.74765, -0.71136, 1.0386, 0.4415, -1.39913, 2.78556, 0.49137, 0.40963, 0.35413, 0.19823, -1.34214, -0.03025, 1.78693, -0.74889, 0.78155, 0.7801, -1.06375, -0.1348, -0.5655, -1.57548, 0.94146, -0.72029, -1.57808, 0.13325, 1.49418, -0.41514, -1.39872, 1.6656, -0.36637, 1.68952, -0.5361, 1.42374, -0.12491, -0.07955, -1.24456, -0.41542, -0.45017, -0.26945, 2.18117, -1.77597, 0.43219, 0.37185, 0.23178, -0.36794, -1.08646, 0.72553, -0.06628, 0.1867, -0.35604, 1.70114, 0.18321, -0.2381, -0.41968, 0.61621, 0.07051, 0.36651, -2.19726, -0.00682, 1.60931, 1.40876, 0.8708, 1.38364, 0.07981, 0.65166, -0.87882, 0.00109, -0.57891, 0.75607, 0.03181, 0.70086, -1.07109, 1.50089, 1.09098, 0.77801, 0.74075, 1.04036, 0.24449, 0.59291, 0.53829, -2.19931, 1.53816, 0.48325, -0.18343, -0.0382, -0.82952, 0.36279, 2.36051, 0.10282, -0.23883, 0.70969, 1.68131, -0.74803, 0.91505, 1.37091, 0.04017, -0.47114, -1.12378, -2.13866, -0.114, -0.34239, -0.19135, 0.24358, -0.69458, 1.33927, 1.47471, -0.0444, -0.12923, -0.5318, 1.99953, 0.47496, 0.64181, -1.09011, 1.33332, -0.0255, 1.73117, -0.02157, 0.7817, 1.68187, -0.62497, 0.40832, 0.24626, 1.02874, -0.05426, -0.78812, -1.36266, 2.07591, 1.17601, -0.05496, 1.14451, 0.0735, 0.27128, 0.48998, 1.30076, 2.10199, 0.26988, -0.09657, -0.89352, 1.34594, 0.29536, -0.35202, -0.02665, 1.40445, 0.10872, -1.46292, 1.18304, 1.12494, -0.2134, 0.5972, -0.42949, 2.31906, -0.05024, -1.23675, -0.11681, -0.22015, 0.70427, 0.30114 };
	//double m2[500] = { 0.59043, 0.73902, 1.35969, 0.4717, -0.55702, 0.28087, -0.79588, 0.39794, 1.57425, 0.61559, 0.77355, -1.85854, -0.17649, 1.06796, -0.87918, -1.25702, -0.32604, -0.23308, -0.35538, 0.48762, 0.21494, 0.69872, 2.54709, -1.02198, -0.86901, 0.92635, 0.11074, 0.03724, -0.40435, 0.40193, 0.29703, 1.18189, -1.41418, 0.10286, -0.3811, 0.47089, 0.30077, -0.03538, 0.03005, -1.10409, 1.13098, 0.24908, -0.9922, -2.40896, -0.04575, 0.06339, 1.00315, -0.74341, 0.26145, 0.32623, 0.62598, 0.29423, 0.68507, 1.10518, 0.70028, -0.5827, 0.02686, -0.96166, 0.22425, 0.54776, -1.57828, 1.90769, -0.03892, 0.10218, -1.34784, 0.76384, -1.03568, -0.01839, 0.78073, -0.18943, 0.02429, 0.66767, -0.55091, -0.48733, 0.57697, -0.52754, 1.20016, 1.03704, -1.2767, -2.78332, -0.09774, 1.37233, -0.28072, 1.24304, -0.44657, -0.60707, 0.80708, -0.95652, -0.73488, -0.1663, 0.29311, 0.42848, 0.88589, 0.18805, -2.01079, -1.25916, 1.03503, -1.08605, 0.04879, -0.18096, -0.11604, -0.19174, 1.34531, -1.46607, 0.55778, -0.24264, 1.11762, -0.84022, -0.18617, -1.66689, -2.44636, 0.21207, -1.45472, 1.06371, 0.05923, 0.49381, 1.35934, -0.19161, -2.02131, 0.03412, 1.19291, -0.73306, 0.76295, -1.07641, 0.85412, 1.16764, 0.99569, 0.14673, 0.91969, 1.32688, -0.08978, -0.95189, -0.346, 0.71143, 0.00157, 1.86995, 1.1947, -0.16909, 0.52508, 0.31753, -0.78823, 0.64893, -0.66965, -0.29468, -0.38501, 1.27336, 0.72214, 0.38427, 1.60488, -0.79293, -0.91361, -0.25728, 0.3311, -1.24378, 0.82229, 0.62758, -2.27309, 0.77812, -0.46307, 0.48409, -0.19249, -0.49241, 0.40484, 0.52682, 0.11376, -1.96371, 0.65966, 1.31904, 0.46204, 1.51428, 0.22505, 0.55186, 1.54479, -0.32964, 0.00956, -0.37837, -0.12939, 1.28034, -0.93739, 0.3324, -0.65837, -0.89723, 0.46224, 0.2937, -0.04323, -0.12029, 0.80673, -1.84476, 0.15184, -1.31891, 0.82292, 0.96672, 0.35302, -0.05274, -0.48815, -0.86286, -1.86579, 0.75751, 0.32389, 0.36092, 0.46473, 0.25634, 0.45608, -2.00679, 0.6134, -1.06691, 0.20428, 0.33744, -1.10731, 0.83024, -0.72591, -0.2822, -0.6684, 0.2243, -1.92791, -0.77893, -1.60536, -0.92753, -1.10094, -0.85466, 0.48333, -0.57182, 0.75615, -0.06397, -1.48506, 1.6857, 1.20692, -1.43756, 0.0759, -0.12987, -0.80226, -1.68567, -0.23085, 0.12065, 1.62835, 0.52099, 0.92947, -1.69064, 0.7264, -1.55475, 0.27229, 0.02162, -0.38901, -0.40211, 0.04171, 0.30869, 1.81743, -0.84171, -0.34589, 1.67376, 0.14069, 0.65859, -0.27017, -0.22909, -0.03824, 0.36274, -0.83248, -1.93866, -0.06766, -0.96975, 1.49781, 0.85687, -0.51741, -1.2366, -0.50588, -0.23323, 0.06055, -0.28679, 0.42638, 0.41716, -0.85856, 1.24501, -0.54944, -2.14373, -0.70306, -1.29784, 0.30303, -1.12085, -0.9563, -0.44079, -1.62482, 1.1295, 1.0375, -0.76894, -1.02874, -0.04127, -0.20329, 2.16398, 1.26837, -1.67849, -0.77004, -0.43606, -0.04388, 0.15593, -0.22608, 0.69059, -0.36136, 0.38308, -0.1214, -0.33405, 0.97822, 0.37216, 0.7743, -0.2597, 0.17904, 0.86038, 0.08493, -1.84777, -0.62562, -0.89226, -0.17838, 0.66459, -0.31715, 0.19586, 1.37068, -0.16929, -0.49169, -0.26572, -0.39948, -1.50093, -1.50519, -1.6248, -0.59873, -0.74251, -0.9206, -0.67699, 0.93957, -1.00901, 0.27976, -0.30027, 0.91222, 0.13222, 0.09407, -0.00448, 0.87424, 0.52401, -0.20452, 1.38131, 0.75715, 0.62505, -0.31664, -1.08908, 0.67025, 0.32057, 0.67536, 0.20215, 2.48281, -0.14602, -1.46856, -0.11344, -0.73773, -1.46799, 0.35057, 0.72066, 0.80357, 0.02462, -0.61924, -0.50905, 0.58542, -0.97892, -0.42476, -0.14205, -0.87877, -0.3683, 0.67675, 1.83454, 0.52829, 0.27148, 1.36463, -0.78017, -0.48732, 0.60565, -0.97053, 0.28047, -0.44951, -0.72792, 0.69981, -0.28562, -0.74792, 0.15289, 0.9952, 0.41952, 0.48554, -2.19451, 0.01642, 0.57162, 2.07191, 2.44355, -0.14499, -1.40803, -0.69705, -0.02143, -0.06746, 1.63641, -0.69513, 0.95718, 0.05925, 0.05848, 0.02399, 2.07293, -0.04626, -0.34551, 1.1117, -2.62044, -0.26545, 1.62817, -0.03607, 0.20276, -0.53158, 1.44436, -0.86126, -1.38459, 0.56972, 0.22826, -1.19282, -1.88144, -0.06037, -0.86134, 0.12943, -0.33455, -1.27462, 1.00022, -1.01642, -2.14409, -0.20117, 1.82033, 0.14003, -0.11077, 0.69637, -0.64455, -0.60274, -0.00191, 0.78019, 0.45039, -0.11985, -0.24543, 1.24644, -1.09524, 0.28411, -0.14997, -0.46216, 0.6948, -0.88775, 0.39214, 0.42642, -1.65611, -1.42844, 0.43892, -1.34412, 0.79907, 1.15016, 1.39097, -0.23473, -1.02691, 0.39496, -0.85243, 0.8836, 0.34019, 0.63432, 0.48672, -0.04512, 0.57848, 0.70554, 0.45539, -0.24069, 1.02068, 0.78303, -0.83741, 0.42064, 0.06586, -0.95693, -0.14122, -0.50722, -2.47783, -0.48123, 0.5985, 0.71655, -0.40602, -0.26169, -0.99237, 1.21474, -0.19264, 0.46304, 0.85428, 0.10044, -0.8947, -0.95167, -0.99662, 0.90768, -0.00923, 0.55399, 0.3596, 0.19553, 0.02627, 0.75873, -0.16044, 1.13287, -0.55721, -0.6246, 1.30546 };
	//n = 500;
	//int R = 10000;
	//type = 1;
	//dst = 0;
	//bcov_test(bd_stat_value, permuted_bcov, m1, m2, &n, &R, &weight, &dst, &type);

	test_multithread_permutation();
	system("pause");
	return;
}

