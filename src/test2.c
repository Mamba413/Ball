#include <stdlib.h>
#include <stdio.h>
#include "BD.h"
#include "BI.h"
#include "bcor.h"


void test_value(double computed_value, double true_value)
{
	const double EPSINON = 0.00001;
	double difference = computed_value - true_value;
	if ((difference < EPSINON) & (difference > -EPSINON)) {
		printf("Success!\n\n");
	}
	else {
		printf("Fail!\n\n");
	}
	return;
}


void test_bd_value()
{
	printf("\n");
	double ball_stat_value[2], p_value[2];
	double xy[20] = { -0.26162, -1.17788, 0.07105, 0.31527, 0.26728, 1.02388, -1.33132, 1.5708, 2.05537, -0.60046, 0.32773, 0.74089, 0.40954, 1.95479, -1.24204, 0.80392, -0.87819, -0.87891, -0.25194, -0.28557 };
	double xy_dst[400] = { 0, 0.91626, 0.33267, 0.57689, 0.5289, 1.2855, 1.0697, 1.83242, 2.31699, 0.33884, 0.58935, 1.00251, 0.67116, 2.21641, 0.98042, 1.06554, 0.61657, 0.61729, 0.00968000000000002, 0.02395, 0.91626, 0, 1.24893, 1.49315, 1.44516, 2.20176, 0.15344, 2.74868, 3.23325, 0.57742, 1.50561, 1.91877, 1.58742, 3.13267, 0.06416, 1.9818, 0.29969, 0.29897, 0.92594, 0.89231, 0.33267, 1.24893, 0, 0.24422, 0.19623, 0.95283, 1.40237, 1.49975, 1.98432, 0.67151, 0.25668, 0.66984, 0.33849, 1.88374, 1.31309, 0.73287, 0.94924, 0.94996, 0.32299, 0.35662, 0.57689, 1.49315, 0.24422, 0, 0.04799, 0.70861, 1.64659, 1.25553, 1.7401, 0.91573, 0.01246, 0.42562, 0.09427, 1.63952, 1.55731, 0.48865, 1.19346, 1.19418, 0.56721, 0.60084, 0.5289, 1.44516, 0.19623, 0.04799, 0, 0.7566, 1.5986, 1.30352, 1.78809, 0.86774, 0.06045, 0.47361, 0.14226, 1.68751, 1.50932, 0.53664, 1.14547, 1.14619, 0.51922, 0.55285, 1.2855, 2.20176, 0.95283, 0.70861, 0.7566, 0, 2.3552, 0.54692, 1.03149, 1.62434, 0.69615, 0.28299, 0.61434, 0.93091, 2.26592, 0.21996, 1.90207, 1.90279, 1.27582, 1.30945, 1.0697, 0.15344, 1.40237, 1.64659, 1.5986, 2.3552, 0, 2.90212, 3.38669, 0.73086, 1.65905, 2.07221, 1.74086, 3.28611, 0.08928, 2.13524, 0.45313, 0.45241, 1.07938, 1.04575, 1.83242, 2.74868, 1.49975, 1.25553, 1.30352, 0.54692, 2.90212, 0, 0.48457, 2.17126, 1.24307, 0.82991, 1.16126, 0.38399, 2.81284, 0.76688, 2.44899, 2.44971, 1.82274, 1.85637, 2.31699, 3.23325, 1.98432, 1.7401, 1.78809, 1.03149, 3.38669, 0.48457, 0, 2.65583, 1.72764, 1.31448, 1.64583, 0.10058, 3.29741, 1.25145, 2.93356, 2.93428, 2.30731, 2.34094, 0.33884, 0.57742, 0.67151, 0.91573, 0.86774, 1.62434, 0.73086, 2.17126, 2.65583, 0, 0.92819, 1.34135, 1.01, 2.55525, 0.64158, 1.40438, 0.27773, 0.27845, 0.34852, 0.31489, 0.58935, 1.50561, 0.25668, 0.01246, 0.06045, 0.69615, 1.65905, 1.24307, 1.72764, 0.92819, 0, 0.41316, 0.08181, 1.62706, 1.56977, 0.47619, 1.20592, 1.20664, 0.57967, 0.6133, 1.00251, 1.91877, 0.66984, 0.42562, 0.47361, 0.28299, 2.07221, 0.82991, 1.31448, 1.34135, 0.41316, 0, 0.33135, 1.2139, 1.98293, 0.0630299999999999, 1.61908, 1.6198, 0.99283, 1.02646, 0.67116, 1.58742, 0.33849, 0.09427, 0.14226, 0.61434, 1.74086, 1.16126, 1.64583, 1.01, 0.08181, 0.33135, 0, 1.54525, 1.65158, 0.39438, 1.28773, 1.28845, 0.66148, 0.69511, 2.21641, 3.13267, 1.88374, 1.63952, 1.68751, 0.93091, 3.28611, 0.38399, 0.10058, 2.55525, 1.62706, 1.2139, 1.54525, 0, 3.19683, 1.15087, 2.83298, 2.8337, 2.20673, 2.24036, 0.98042, 0.06416, 1.31309, 1.55731, 1.50932, 2.26592, 0.08928, 2.81284, 3.29741, 0.64158, 1.56977, 1.98293, 1.65158, 3.19683, 0, 2.04596, 0.36385, 0.36313, 0.9901, 0.95647, 1.06554, 1.9818, 0.73287, 0.48865, 0.53664, 0.21996, 2.13524, 0.76688, 1.25145, 1.40438, 0.47619, 0.0630299999999999, 0.39438, 1.15087, 2.04596, 0, 1.68211, 1.68283, 1.05586, 1.08949, 0.61657, 0.29969, 0.94924, 1.19346, 1.14547, 1.90207, 0.45313, 2.44899, 2.93356, 0.27773, 1.20592, 1.61908, 1.28773, 2.83298, 0.36385, 1.68211, 0, 0.000719999999999943, 0.62625, 0.59262, 0.61729, 0.29897, 0.94996, 1.19418, 1.14619, 1.90279, 0.45241, 2.44971, 2.93428, 0.27845, 1.20664, 1.6198, 1.28845, 2.8337, 0.36313, 1.68283, 0.000719999999999943, 0, 0.62697, 0.59334, 0.00968000000000002, 0.92594, 0.32299, 0.56721, 0.51922, 1.27582, 1.07938, 1.82274, 2.30731, 0.34852, 0.57967, 0.99283, 0.66148, 2.20673, 0.9901, 1.05586, 0.62625, 0.62697, 0, 0.03363, 0.02395, 0.89231, 0.35662, 0.60084, 0.55285, 1.30945, 1.04575, 1.85637, 2.34094, 0.31489, 0.6133, 1.02646, 0.69511, 2.24036, 0.95647, 1.08949, 0.59262, 0.59334, 0.03363, 0 };
	int size[2] = { 10, 10 };
	int n = 20, R = 0, K = 2, dst = 0, nth = 1;

	// Software of ChengFeng Liu output result (golden standard): 0.0274 
	// univariate case:
	bd_test(ball_stat_value, p_value, xy, size, &n, &K, &dst, &R, &nth);
	printf("Univariate Ball Divergence: %f; \n", ball_stat_value[0]);
	test_value(ball_stat_value[0], 0.0274);

	// multivariate case:
	dst = 1;
	bd_test(ball_stat_value, p_value, xy_dst, size, &n, &K, &dst, &R, &nth);
	printf("Multivariate Ball Divergence: %f; \n", ball_stat_value[0]);
	test_value(ball_stat_value[0], 0.0274);
}


void test_bcov_value()
{
	printf("\n");
	double ball_stat_value[3], p_value[3];
	double x[10] = { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10 };
	double y[10] = { 11, 12, 13, 14, 15, 16, 17, 18, 19, 20 };
	double x_dst[100] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 2, 1, 0, 1, 2, 3, 4, 5, 6, 7, 3, 2, 1, 0, 1, 2, 3, 4, 5, 6, 4, 3, 2, 1, 0, 1, 2, 3, 4, 5, 5, 4, 3, 2, 1, 0, 1, 2, 3, 4, 6, 5, 4, 3, 2, 1, 0, 1, 2, 3, 7, 6, 5, 4, 3, 2, 1, 0, 1, 2, 8, 7, 6, 5, 4, 3, 2, 1, 0, 1, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0 };
	double y_dst[100] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 2, 1, 0, 1, 2, 3, 4, 5, 6, 7, 3, 2, 1, 0, 1, 2, 3, 4, 5, 6, 4, 3, 2, 1, 0, 1, 2, 3, 4, 5, 5, 4, 3, 2, 1, 0, 1, 2, 3, 4, 6, 5, 4, 3, 2, 1, 0, 1, 2, 3, 7, 6, 5, 4, 3, 2, 1, 0, 1, 2, 8, 7, 6, 5, 4, 3, 2, 1, 0, 1, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0 };
	int n = 10;
	int R = 0, weight = 0, dst = 0, type = 1, nth = 1;

	// Software of ChengFeng Liu output result (golden standard): 0.034214
	// univariate case:
	bcov_test(ball_stat_value, p_value, x, y, &n, &R, &weight, &dst, &type, &nth);
	printf("Univariate Ball Covariance: %f; \n", ball_stat_value[0]);
	test_value(ball_stat_value[0], 0.034214);

	// multivariate case:
	dst = 1;
	bcov_test(ball_stat_value, p_value, x_dst, y_dst, &n, &R, &weight, &dst, &type, &nth);
	printf("Multivariate Ball Covariance: %f; \n", ball_stat_value[0]);
	test_value(ball_stat_value[0], 0.034214);
}


void test_bcor_value() {
	printf("\n");
	double ball_stat_value[6];
	double x[20] = { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1 };
	double y[10] = { 11, 12, 13, 14, 15, 16, 17, 18, 19, 20 };
	int x_number[2] = { 1, 1 };
	int f_number = 2, n = 10, weight = 0, dst_y = 0, dst_x = 0, k = 0, p = 1, nth = 2, size = 2;

	// univariate case:
	printf("----- First Covariance ----- \n");
	bcor_test(ball_stat_value, y, x, &x_number, &f_number, &size, &n, &p, &k, &dst_y, &dst_x, &nth);
	printf("Univariate Ball Correlation: %f; \n", ball_stat_value[0]);
	test_value(ball_stat_value[0], 1);
	printf("Univariate Weight Ball Correlation: %f; \n", ball_stat_value[1]);
	test_value(ball_stat_value[1], 1);
	printf("Univariate HHG Ball Correlation: %f; \n", ball_stat_value[2]);
	test_value(ball_stat_value[2], 1);

	printf("----- Second Covariance ----- \n");
	printf("Univariate Ball Correlation: %f; \n", ball_stat_value[3]);
	test_value(ball_stat_value[3], 1);
	printf("Univariate Weight Ball Correlation: %f; \n", ball_stat_value[4]);
	test_value(ball_stat_value[4], 1);
	printf("Univariate HHG Ball Correlation: %f; \n", ball_stat_value[5]);
	test_value(ball_stat_value[5], 1);

	// multivariate case:
	double y_dst[100] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 2, 1, 0, 1, 2, 3, 4, 5, 6, 7, 3, 2, 1, 0, 1, 2, 3, 4, 5, 6, 4, 3, 2, 1, 0, 1, 2, 3, 4, 5, 5, 4, 3, 2, 1, 0, 1, 2, 3, 4, 6, 5, 4, 3, 2, 1, 0, 1, 2, 3, 7, 6, 5, 4, 3, 2, 1, 0, 1, 2, 8, 7, 6, 5, 4, 3, 2, 1, 0, 1, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0 };
	dst_y = 1;
	p = 0;
	printf("----- First Covariance ----- \n");
	bcor_test(ball_stat_value, y_dst, x, &x_number, &f_number, &size, &n, &p, &k, &dst_y, &dst_x, &nth);
	printf("Multivariate Ball Correlation: %f; \n", ball_stat_value[0]);
	test_value(ball_stat_value[0], 1);
	printf("Multivariate Weight Ball Correlation: %f; \n", ball_stat_value[1]);
	test_value(ball_stat_value[1], 1);
	printf("Multivariate HHG Ball Correlation: %f; \n", ball_stat_value[2]);
	test_value(ball_stat_value[2], 1);

	printf("----- Second Covariance ----- \n");
	printf("Multivariate Ball Correlation: %f; \n", ball_stat_value[3]);
	test_value(ball_stat_value[3], 1);
	printf("Multivariate Weight Ball Correlation: %f; \n", ball_stat_value[4]);
	test_value(ball_stat_value[4], 1);
	printf("Multivariate HHG Ball Correlation: %f; \n", ball_stat_value[5]);
	test_value(ball_stat_value[5], 1);

	// distance matrix case:
	double x_dst[100] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 2, 1, 0, 1, 2, 3, 4, 5, 6, 7, 3, 2, 1, 0, 1, 2, 3, 4, 5, 6, 4, 3, 2, 1, 0, 1, 2, 3, 4, 5, 5, 4, 3, 2, 1, 0, 1, 2, 3, 4, 6, 5, 4, 3, 2, 1, 0, 1, 2, 3, 7, 6, 5, 4, 3, 2, 1, 0, 1, 2, 8, 7, 6, 5, 4, 3, 2, 1, 0, 1, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0 };
	dst_x = 1;
	bcor_test(ball_stat_value, y_dst, x_dst, &x_number, &f_number, &size, &n, &p, &k, &dst_y, &dst_x, &nth);
	printf("Distance Based Ball Correlation: %f; \n", ball_stat_value[0]);
	test_value(ball_stat_value[0], 1);
	printf("Distance Based Weight Ball Correlation: %f; \n", ball_stat_value[1]);
	test_value(ball_stat_value[1], 1);
	printf("Distance Based HHG Ball Correlation: %f; \n", ball_stat_value[2]);
	test_value(ball_stat_value[2], 1);

	return;
}


void bd_test_2sample_multithread_permutation()
{
	printf("----- Test Ball Divergence based Two Sample Test -----\n");
	double ball_stat_value[2], p_value[2];
	double x[500] = { 2.25175, 1.55615, 0.54524, 2.85418, 1.10683, 0.72145, 1.97587, 0.84344, -1.70603, -0.78955, 2.32068, 0.97631, 0.77115, 1.82555, -0.40834, 0.70672, 1.35385, 1.65415, 0.54051, 1.76525, 0.64433, 1.7401, 0.16161, 1.713, 1.08914, 0.39928, 1.47068, 1.04822, 0.74917, 0.2367, 1.22338, 1.93922, -0.84661, 0.0847, 2.92419, 0.57285, 0.84396, -1.94022, 0.14856, 1.23254, 1.07918, 1.79337, 2.72291, 0.63817, 1.01389, 0.95096, 1.5309, 1.28475, 0.56053, 0.50629, 0.18289, 1.03813, 1.24354, 1.53103, 0.14486, 0.97716, 0.91053, 0.52578, 0.72429, 0.43408, 1.67109, 1.27183, 1.69888, 0.76495, 0.67365, 0.35637, 0.01117, 1.18014, -0.3262, 0.64187, -0.02459, 0.31699, 2.20765, -0.61115, 0.16935, 0.96341, -0.0671, 1.43816, 0.43653, 2.5803, 1.36074, 0.10524, -0.4532, 0.908, 1.71979, 2.05462, 0.28815, 0.63903, 0.75906, 2.10066, 1.64896, 2.14244, 1.01404, 2.49235, 2.23282, 1.34103, 0.66155, 2.235, -0.35812, 1.51349, 0.90529, 1.14303, 0.36791, 1.37934, -1.43776, 1.07899, 0.80347, 2.65376, 1.75877, 1.57104, 0.31, 2.09237, 0.38631, 1.6475, 1.99678, 0.81316, 0.57989, -0.02215, 0.99435, 0.40954, 2.20811, 1.479, 1.49441, 0.0701, 0.87551, 0.39109, 1.29121, 0.54331, 1.66132, 1.1505, 3.23168, 1.58024, 2.08576, 2.66405, 1.79259, -0.953, 1.70187, -0.11666, 0.59768, -0.06712, 1.8215, 0.48039, -0.75134, 0.13491, 1.53026, 0.38776, 0.54211, 2.80549, 2.10628, 1.96142, 1.22193, 3.45196, 0.52233, 1.14369, 2.39459, 0.12202, -0.04252, 1.64233, 1.86858, 1.62409, 1.74682, -1.05572, 1.42309, 0.70916, -0.75358, 0.99209, -0.88119, -0.0581, -0.76142, 1.67469, 1.88033, 1.38972, 1.00137, 1.3897, 1.98954, 2.12393, 0.49024, 1.65855, 0.49964, 1.92137, 2.68253, 1.55927, -0.07893, 0.5343, -0.91064, 1.17404, 2.10121, 0.43845, 0.73953, 0.66175, 1.7428, 0.53864, 1.02051, 2.66725, 0.92197, 0.22764, -0.25545, -0.17605, -0.4281, -0.24533, 0.62213, 0.91474, -0.58603, 0.36518, 1.19489, 0.18535, 0.08253, 0.77608, 2.37213, 0.58334, 0.3253, 1.99076, 0.49968, 1.79846, -0.45763, 1.63836, 1.30264, 1.92182, 2.48508, 1.13147, -0.19058, 0.88212, -1.01598, 0.77181, -0.31026, 0.62621, 0.16373, 0.2404, 1.06689, -0.71961, 0.45225, 2.05092, 1.7834, 0.76702, -0.41365, 0.98748, -1.33853, 0.24561, 0.06651, -0.27744, -0.76812, -0.41024, -0.38618, 0.32252, 1.63112, -0.06034, 0.25776, -1.54519, 0.24131, 1.54197, 0.51826, 1.97117, 1.91757, -0.67011, 0.41305, -0.05486, -0.35572, 2.02861, 0.23772, 1.64604, -0.29029, -0.18933, 1.12945, -0.29203, -0.12381, 0.38669, -0.48716, 0.18902, 1.0314, -0.64156, -0.99785, -0.95451, -1.00394, -1.07203, 0.79471, 3.1304, 0.87868, 0.62735, -0.00596, 0.56395, -0.16823, 0.581, 0.6356, 0.94813, -0.92096, 1.66221, 0.16091, 0.48368, -0.16963, 0.30523, -1.48444, -1.43197, 0.00041, -0.96845, 0.22431, 0.42396, -1.2707, 0.30458, 2.58439, 1.01013, 0.02986, -0.16052, 0.73194, -0.83472, 0.03374, -0.4647, -0.49565, 0.73078, 1.3263, -1.43352, 1.14522, 0.87447, 0.59477, -2.91983, 0.9821, 0.04228, -0.7579, 0.72066, 0.63426, -1.98498, -1.84243, 0.23127, 1.36689, 0.65751, -1.32793, 0.76658, 0.91408, -0.36984, 0.24676, -1.32251, 1.11424, 0.28732, -1.71797, 0.68414, 0.35187, -0.38769, 0.23127, -1.53572, 1.3105, 0.84261, 0.84039, -0.39682, -0.59393, -0.4193, 0.78651, -1.87045, -0.05073, 1.12713, -1.13818, 0.41322, -0.85635, 0.85415, 0.09091, -1.02064, -0.86351, 0.97995, -0.44872, 1.56491, 1.01164, -0.31884, 0.50394, 1.38312, -0.29272, 1.17322, 1.4049, -0.33293, 2.18284, -0.1389, 1.5242, 1.25374, -0.55504, -1.91363, 0.45471, -1.07961, -0.69629, -1.45287, -0.48276, -0.29733, -0.1064, 1.7312, 0.64748, -1.12057, 0.6634, -0.48684, -0.43941, 0.40728, -0.07195, 0.61082, 0.68196, -2.13277, 0.30082, -0.58186, 0.29594, -0.93074, 0.95474, -0.80438, 2.12899, 0.48902, 0.81091, 0.4143, -0.12889, 1.50433, 1.88676, 0.89654, 1.01675, 0.00424, 1.42727, -0.72508, 1.05808, -0.02322, 0.87603, 1.03567, 0.21079, -0.09184, -0.90076, -0.4944, -1.22989, 0.78523, -0.14292, 0.1735, -0.87049, 1.01657, -0.07105, -0.77789, -1.67587, -0.4794, -0.50866, 1.81928, -0.78227, 0.93179, 1.11011, 0.53856, 1.00276, -1.01718, -0.72113, -1.22834, 0.16478, -1.49958, -0.01627, 0.61458, 0.30104, -0.15816, -1.32077, -2.45525, 0.43569, 0.37217, -0.9515, 0.30982, -1.27753, 0.48857, -0.33971, 0.01417, -1.51227, -1.08933, 1.11125, 1.00337, -0.55228, -2.49924, -0.56945, 1.857, -0.76726, -0.5512, 0.33967, 0.45384, -0.42878, 0.70092, 0.84725, -0.08575, 0.96277, -0.88731, -0.05063, 0.56627, -0.18277, 0.00336, 0.60146, -1.4978, 0.99308, -0.59198, 0.00807, 0.33094, 0.67669, -1.58132, -0.05644, 0.9292, -0.94089, -0.62233, -0.14153, -1.2477, -0.13975, -1.28383, 0.59095, 0.42241, -1.60643, 0.10822, -0.50155, -0.3301, 1.86959, 0.30328, -0.15598, -2.63078 };
	int size[2] = { 250, 250 };
	int k = 2;
	int n = 500;
	int dst = 0, R = 199, nth = 2;
	int parallel_type = 2;
	bd_test(ball_stat_value, p_value, x, size, &n, &k, &dst, &R, &nth);
	printf("Ball statistics: %f; ", ball_stat_value[0]);
	printf("p-value: %f \n", p_value[0]);
	printf("Ball statistics: %f; ", ball_stat_value[1]);
	printf("p-value: %f \n", p_value[1]);

	// Software of ChengFeng Liu output result (golden standard): 0.059821 
	test_value(ball_stat_value[0], 0.059821);
	return;
}


void bd_test_ksample_multithread_permutation()
{
	printf("Test Ball Divergence based K Sample Test\n");
	double ball_stat_value[4], p_value[4];
	double x[15] = { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15 };
	int size[3] = { 5, 5, 5};
	int k = 3;
	int n = 15;
	int dst = 0, R = 199, nth = 2;
	int parallel_type = 2;
	bd_test(ball_stat_value, p_value, x, size, &n, &k, &dst, &R, &nth);
	printf("Ball statistics: %f; ", ball_stat_value[0]);
	printf("p-value: %f \n", p_value[0]);
	printf("Ball statistics: %f; ", ball_stat_value[1]);
	printf("p-value: %f \n", p_value[1]);
	printf("Ball statistics: %f; ", ball_stat_value[2]);
	printf("p-value: %f \n", p_value[2]);
	printf("Ball statistics: %f; ", ball_stat_value[3]);
	printf("p-value: %f \n", p_value[3]);
	printf("No need to check the K-sample test statistic! \n\n");
	return;
}


void bcov_test_multithread_permutaion()
{
	printf("----- Test Ball Covariance based Independence Test -----\n");
	double ball_stat_value[3], p_value[3];
	double x_dst[900] = { 0, 1.4437, 1.59269, 1.96331, 1.66324, 0.24915, 0.7639, 1.02009, 0.99903, 2.8877, 1.1999, 1.2857, 1.2278, 1.79689, 0.89516, 1.68374, 0.19083, 0.0557, 2.55586, 1.36259, 0.70611, 0.31423, 0.14226, 0.54568, 0.1566, 0.92097, 0.58167, 0.5148, 1.75494, 2.16156, 1.4437, 0, 0.14899, 0.51961, 0.21954, 1.19455, 0.6798, 0.42361, 2.44273, 1.444, 2.6436, 0.158, 0.2159, 0.35319, 0.54854, 0.24004, 1.25287, 1.4994, 1.11216, 0.08111, 2.14981, 1.75793, 1.30144, 1.98938, 1.6003, 0.52273, 0.86203, 0.9289, 0.31124, 0.71786, 1.59269, 0.14899, 0, 0.37062, 0.07055, 1.34354, 0.82879, 0.5726, 2.59172, 1.29501, 2.79259, 0.30699, 0.36489, 0.2042, 0.69753, 0.09105, 1.40186, 1.64839, 0.96317, 0.2301, 2.2988, 1.90692, 1.45043, 2.13837, 1.74929, 0.67172, 1.01102, 1.07789, 0.16225, 0.56887, 1.96331, 0.51961, 0.37062, 0, 0.30007, 1.71416, 1.19941, 0.94322, 2.96234, 0.92439, 3.16321, 0.67761, 0.73551, 0.16642, 1.06815, 0.27957, 1.77248, 2.01901, 0.59255, 0.60072, 2.66942, 2.27754, 1.82105, 2.50899, 2.11991, 1.04234, 1.38164, 1.44851, 0.20837, 0.19825, 1.66324, 0.21954, 0.07055, 0.30007, 0, 1.41409, 0.89934, 0.64315, 2.66227, 1.22446, 2.86314, 0.37754, 0.43544, 0.13365, 0.76808, 0.0205, 1.47241, 1.71894, 0.89262, 0.30065, 2.36935, 1.97747, 1.52098, 2.20892, 1.81984, 0.74227, 1.08157, 1.14844, 0.0917, 0.49832, 0.24915, 1.19455, 1.34354, 1.71416, 1.41409, 0, 0.51475, 0.77094, 1.24818, 2.63855, 1.44905, 1.03655, 0.97865, 1.54774, 0.64601, 1.43459, 0.05832, 0.30485, 2.30671, 1.11344, 0.95526, 0.56338, 0.10689, 0.79483, 0.40575, 0.67182, 0.33252, 0.26565, 1.50579, 1.91241, 0.7639, 0.6798, 0.82879, 1.19941, 0.89934, 0.51475, 0, 0.25619, 1.76293, 2.1238, 1.9638, 0.5218, 0.4639, 1.03299, 0.13126, 0.91984, 0.57307, 0.8196, 1.79196, 0.59869, 1.47001, 1.07813, 0.62164, 1.30958, 0.9205, 0.15707, 0.18223, 0.2491, 0.99104, 1.39766, 1.02009, 0.42361, 0.5726, 0.94322, 0.64315, 0.77094, 0.25619, 0, 2.01912, 1.86761, 2.21999, 0.26561, 0.20771, 0.7768, 0.12493, 0.66365, 0.82926, 1.07579, 1.53577, 0.3425, 1.7262, 1.33432, 0.87783, 1.56577, 1.17669, 0.09912, 0.43842, 0.50529, 0.73485, 1.14147, 0.99903, 2.44273, 2.59172, 2.96234, 2.66227, 1.24818, 1.76293, 2.01912, 0, 3.88673, 0.20087, 2.28473, 2.22683, 2.79592, 1.89419, 2.68277, 1.18986, 0.94333, 3.55489, 2.36162, 0.29292, 0.6848, 1.14129, 0.45335, 0.84243, 1.92, 1.5807, 1.51383, 2.75397, 3.16059, 2.8877, 1.444, 1.29501, 0.92439, 1.22446, 2.63855, 2.1238, 1.86761, 3.88673, 0, 4.0876, 1.602, 1.6599, 1.09081, 1.99254, 1.20396, 2.69687, 2.9434, 0.33184, 1.52511, 3.59381, 3.20193, 2.74544, 3.43338, 3.0443, 1.96673, 2.30603, 2.3729, 1.13276, 0.72614, 1.1999, 2.6436, 2.79259, 3.16321, 2.86314, 1.44905, 1.9638, 2.21999, 0.20087, 4.0876, 0, 2.4856, 2.4277, 2.99679, 2.09506, 2.88364, 1.39073, 1.1442, 3.75576, 2.56249, 0.49379, 0.88567, 1.34216, 0.65422, 1.0433, 2.12087, 1.78157, 1.7147, 2.95484, 3.36146, 1.2857, 0.158, 0.30699, 0.67761, 0.37754, 1.03655, 0.5218, 0.26561, 2.28473, 1.602, 2.4856, 0, 0.0579, 0.51119, 0.39054, 0.39804, 1.09487, 1.3414, 1.27016, 0.07689, 1.99181, 1.59993, 1.14344, 1.83138, 1.4423, 0.36473, 0.70403, 0.7709, 0.46924, 0.87586, 1.2278, 0.2159, 0.36489, 0.73551, 0.43544, 0.97865, 0.4639, 0.20771, 2.22683, 1.6599, 2.4277, 0.0579, 0, 0.56909, 0.33264, 0.45594, 1.03697, 1.2835, 1.32806, 0.13479, 1.93391, 1.54203, 1.08554, 1.77348, 1.3844, 0.30683, 0.64613, 0.713, 0.52714, 0.93376, 1.79689, 0.35319, 0.2042, 0.16642, 0.13365, 1.54774, 1.03299, 0.7768, 2.79592, 1.09081, 2.99679, 0.51119, 0.56909, 0, 0.90173, 0.11315, 1.60606, 1.85259, 0.75897, 0.4343, 2.503, 2.11112, 1.65463, 2.34257, 1.95349, 0.87592, 1.21522, 1.28209, 0.04195, 0.36467, 0.89516, 0.54854, 0.69753, 1.06815, 0.76808, 0.64601, 0.13126, 0.12493, 1.89419, 1.99254, 2.09506, 0.39054, 0.33264, 0.90173, 0, 0.78858, 0.70433, 0.95086, 1.6607, 0.46743, 1.60127, 1.20939, 0.7529, 1.44084, 1.05176, 0.02581, 0.31349, 0.38036, 0.85978, 1.2664, 1.68374, 0.24004, 0.09105, 0.27957, 0.0205, 1.43459, 0.91984, 0.66365, 2.68277, 1.20396, 2.88364, 0.39804, 0.45594, 0.11315, 0.78858, 0, 1.49291, 1.73944, 0.87212, 0.32115, 2.38985, 1.99797, 1.54148, 2.22942, 1.84034, 0.76277, 1.10207, 1.16894, 0.0712, 0.47782, 0.19083, 1.25287, 1.40186, 1.77248, 1.47241, 0.05832, 0.57307, 0.82926, 1.18986, 2.69687, 1.39073, 1.09487, 1.03697, 1.60606, 0.70433, 1.49291, 0, 0.24653, 2.36503, 1.17176, 0.89694, 0.50506, 0.04857, 0.73651, 0.34743, 0.73014, 0.39084, 0.32397, 1.56411, 1.97073, 0.0557, 1.4994, 1.64839, 2.01901, 1.71894, 0.30485, 0.8196, 1.07579, 0.94333, 2.9434, 1.1442, 1.3414, 1.2835, 1.85259, 0.95086, 1.73944, 0.24653, 0, 2.61156, 1.41829, 0.65041, 0.25853, 0.19796, 0.48998, 0.1009, 0.97667, 0.63737, 0.5705, 1.81064, 2.21726, 2.55586, 1.11216, 0.96317, 0.59255, 0.89262, 2.30671, 1.79196, 1.53577, 3.55489, 0.33184, 3.75576, 1.27016, 1.32806, 0.75897, 1.6607, 0.87212, 2.36503, 2.61156, 0, 1.19327, 3.26197, 2.87009, 2.4136, 3.10154, 2.71246, 1.63489, 1.97419, 2.04106, 0.80092, 0.3943, 1.36259, 0.08111, 0.2301, 0.60072, 0.30065, 1.11344, 0.59869, 0.3425, 2.36162, 1.52511, 2.56249, 0.07689, 0.13479, 0.4343, 0.46743, 0.32115, 1.17176, 1.41829, 1.19327, 0, 2.0687, 1.67682, 1.22033, 1.90827, 1.51919, 0.44162, 0.78092, 0.84779, 0.39235, 0.79897, 0.70611, 2.14981, 2.2988, 2.66942, 2.36935, 0.95526, 1.47001, 1.7262, 0.29292, 3.59381, 0.49379, 1.99181, 1.93391, 2.503, 1.60127, 2.38985, 0.89694, 0.65041, 3.26197, 2.0687, 0, 0.39188, 0.84837, 0.16043, 0.54951, 1.62708, 1.28778, 1.22091, 2.46105, 2.86767, 0.31423, 1.75793, 1.90692, 2.27754, 1.97747, 0.56338, 1.07813, 1.33432, 0.6848, 3.20193, 0.88567, 1.59993, 1.54203, 2.11112, 1.20939, 1.99797, 0.50506, 0.25853, 2.87009, 1.67682, 0.39188, 0, 0.45649, 0.23145, 0.15763, 1.2352, 0.8959, 0.82903, 2.06917, 2.47579, 0.14226, 1.30144, 1.45043, 1.82105, 1.52098, 0.10689, 0.62164, 0.87783, 1.14129, 2.74544, 1.34216, 1.14344, 1.08554, 1.65463, 0.7529, 1.54148, 0.04857, 0.19796, 2.4136, 1.22033, 0.84837, 0.45649, 0, 0.68794, 0.29886, 0.77871, 0.43941, 0.37254, 1.61268, 2.0193, 0.54568, 1.98938, 2.13837, 2.50899, 2.20892, 0.79483, 1.30958, 1.56577, 0.45335, 3.43338, 0.65422, 1.83138, 1.77348, 2.34257, 1.44084, 2.22942, 0.73651, 0.48998, 3.10154, 1.90827, 0.16043, 0.23145, 0.68794, 0, 0.38908, 1.46665, 1.12735, 1.06048, 2.30062, 2.70724, 0.1566, 1.6003, 1.74929, 2.11991, 1.81984, 0.40575, 0.9205, 1.17669, 0.84243, 3.0443, 1.0433, 1.4423, 1.3844, 1.95349, 1.05176, 1.84034, 0.34743, 0.1009, 2.71246, 1.51919, 0.54951, 0.15763, 0.29886, 0.38908, 0, 1.07757, 0.73827, 0.6714, 1.91154, 2.31816, 0.92097, 0.52273, 0.67172, 1.04234, 0.74227, 0.67182, 0.15707, 0.09912, 1.92, 1.96673, 2.12087, 0.36473, 0.30683, 0.87592, 0.02581, 0.76277, 0.73014, 0.97667, 1.63489, 0.44162, 1.62708, 1.2352, 0.77871, 1.46665, 1.07757, 0, 0.3393, 0.40617, 0.83397, 1.24059, 0.58167, 0.86203, 1.01102, 1.38164, 1.08157, 0.33252, 0.18223, 0.43842, 1.5807, 2.30603, 1.78157, 0.70403, 0.64613, 1.21522, 0.31349, 1.10207, 0.39084, 0.63737, 1.97419, 0.78092, 1.28778, 0.8959, 0.43941, 1.12735, 0.73827, 0.3393, 0, 0.06687, 1.17327, 1.57989, 0.5148, 0.9289, 1.07789, 1.44851, 1.14844, 0.26565, 0.2491, 0.50529, 1.51383, 2.3729, 1.7147, 0.7709, 0.713, 1.28209, 0.38036, 1.16894, 0.32397, 0.5705, 2.04106, 0.84779, 1.22091, 0.82903, 0.37254, 1.06048, 0.6714, 0.40617, 0.06687, 0, 1.24014, 1.64676, 1.75494, 0.31124, 0.16225, 0.20837, 0.0917, 1.50579, 0.99104, 0.73485, 2.75397, 1.13276, 2.95484, 0.46924, 0.52714, 0.04195, 0.85978, 0.0712, 1.56411, 1.81064, 0.80092, 0.39235, 2.46105, 2.06917, 1.61268, 2.30062, 1.91154, 0.83397, 1.17327, 1.24014, 0, 0.40662, 2.16156, 0.71786, 0.56887, 0.19825, 0.49832, 1.91241, 1.39766, 1.14147, 3.16059, 0.72614, 3.36146, 0.87586, 0.93376, 0.36467, 1.2664, 0.47782, 1.97073, 2.21726, 0.3943, 0.79897, 2.86767, 2.47579, 2.0193, 2.70724, 2.31816, 1.24059, 1.57989, 1.64676, 0.40662, 0 };
	double y_dst[900] = { 0, 1.82482, 0.88996, 2.01347, 0.77075, 1.97642, 1.9998, 0.96979, 0.82214, 1.70271, 0.82916, 0.46743, 0.12047, 0.92949, 0.57493, 1.24174, 2.48968, 1.26487, 2.59411, 3.70077, 3.28537, 0.73034, 1.20357, 0.86433, 0.47562, 1.7931, 2.20678, 1.75609, 1.03666, 2.79814, 1.82482, 0, 0.93486, 0.18865, 1.05407, 0.1516, 0.17498, 0.85503, 1.00268, 0.12211, 0.99566, 1.35739, 1.70435, 0.89533, 1.24989, 3.06656, 0.66486, 0.55995, 0.76929, 1.87595, 1.46055, 1.09448, 0.62125, 0.96049, 1.3492, 0.03172, 0.38196, 0.06873, 0.78816, 0.97332, 0.88996, 0.93486, 0, 1.12351, 0.11921, 1.08646, 1.10984, 0.07983, 0.06782, 0.81275, 0.0608, 0.42253, 0.76949, 0.03953, 0.31503, 2.1317, 1.59972, 0.37491, 1.70415, 2.81081, 2.39541, 0.15962, 0.31361, 0.02563, 0.41434, 0.90314, 1.31682, 0.86613, 0.1467, 1.90818, 2.01347, 0.18865, 1.12351, 0, 1.24272, 0.03705, 0.01367, 1.04368, 1.19133, 0.31076, 1.18431, 1.54604, 1.893, 1.08398, 1.43854, 3.25521, 0.47621, 0.7486, 0.58064, 1.6873, 1.2719, 1.28313, 0.8099, 1.14914, 1.53785, 0.22037, 0.19331, 0.25738, 0.97681, 0.78467, 0.77075, 1.05407, 0.11921, 1.24272, 0, 1.20567, 1.22905, 0.19904, 0.05139, 0.93196, 0.05841, 0.30332, 0.65028, 0.15874, 0.19582, 2.01249, 1.71893, 0.49412, 1.82336, 2.93002, 2.51462, 0.04041, 0.43282, 0.09358, 0.29513, 1.02235, 1.43603, 0.98534, 0.26591, 2.02739, 1.97642, 0.1516, 1.08646, 0.03705, 1.20567, 0, 0.02338, 1.00663, 1.15428, 0.27371, 1.14726, 1.50899, 1.85595, 1.04693, 1.40149, 3.21816, 0.51326, 0.71155, 0.61769, 1.72435, 1.30895, 1.24608, 0.77285, 1.11209, 1.5008, 0.18332, 0.23036, 0.22033, 0.93976, 0.82172, 1.9998, 0.17498, 1.10984, 0.01367, 1.22905, 0.02338, 0, 1.03001, 1.17766, 0.29709, 1.17064, 1.53237, 1.87933, 1.07031, 1.42487, 3.24154, 0.48988, 0.73493, 0.59431, 1.70097, 1.28557, 1.26946, 0.79623, 1.13547, 1.52418, 0.2067, 0.20698, 0.24371, 0.96314, 0.79834, 0.96979, 0.85503, 0.07983, 1.04368, 0.19904, 1.00663, 1.03001, 0, 0.14765, 0.73292, 0.14063, 0.50236, 0.84932, 0.0403, 0.39486, 2.21153, 1.51989, 0.29508, 1.62432, 2.73098, 2.31558, 0.23945, 0.23378, 0.10546, 0.49417, 0.82331, 1.23699, 0.7863, 0.06687, 1.82835, 0.82214, 1.00268, 0.06782, 1.19133, 0.05139, 1.15428, 1.17766, 0.14765, 0, 0.88057, 0.00702, 0.35471, 0.70167, 0.10735, 0.24721, 2.06388, 1.66754, 0.44273, 1.77197, 2.87863, 2.46323, 0.0918, 0.38143, 0.04219, 0.34652, 0.97096, 1.38464, 0.93395, 0.21452, 1.976, 1.70271, 0.12211, 0.81275, 0.31076, 0.93196, 0.27371, 0.29709, 0.73292, 0.88057, 0, 0.87355, 1.23528, 1.58224, 0.77322, 1.12778, 2.94445, 0.78697, 0.43784, 0.8914, 1.99806, 1.58266, 0.97237, 0.49914, 0.83838, 1.22709, 0.09039, 0.50407, 0.05338, 0.66605, 1.09543, 0.82916, 0.99566, 0.0608, 1.18431, 0.05841, 1.14726, 1.17064, 0.14063, 0.00702, 0.87355, 0, 0.36173, 0.70869, 0.10033, 0.25423, 2.0709, 1.66052, 0.43571, 1.76495, 2.87161, 2.45621, 0.09882, 0.37441, 0.03517, 0.35354, 0.96394, 1.37762, 0.92693, 0.2075, 1.96898, 0.46743, 1.35739, 0.42253, 1.54604, 0.30332, 1.50899, 1.53237, 0.50236, 0.35471, 1.23528, 0.36173, 0, 0.34696, 0.46206, 0.1075, 1.70917, 2.02225, 0.79744, 2.12668, 3.23334, 2.81794, 0.26291, 0.73614, 0.3969, 0.00819, 1.32567, 1.73935, 1.28866, 0.56923, 2.33071, 0.12047, 1.70435, 0.76949, 1.893, 0.65028, 1.85595, 1.87933, 0.84932, 0.70167, 1.58224, 0.70869, 0.34696, 0, 0.80902, 0.45446, 1.36221, 2.36921, 1.1444, 2.47364, 3.5803, 3.1649, 0.60987, 1.0831, 0.74386, 0.35515, 1.67263, 2.08631, 1.63562, 0.91619, 2.67767, 0.92949, 0.89533, 0.03953, 1.08398, 0.15874, 1.04693, 1.07031, 0.0403, 0.10735, 0.77322, 0.10033, 0.46206, 0.80902, 0, 0.35456, 2.17123, 1.56019, 0.33538, 1.66462, 2.77128, 2.35588, 0.19915, 0.27408, 0.06516, 0.45387, 0.86361, 1.27729, 0.8266, 0.10717, 1.86865, 0.57493, 1.24989, 0.31503, 1.43854, 0.19582, 1.40149, 1.42487, 0.39486, 0.24721, 1.12778, 0.25423, 0.1075, 0.45446, 0.35456, 0, 1.81667, 1.91475, 0.68994, 2.01918, 3.12584, 2.71044, 0.15541, 0.62864, 0.2894, 0.09931, 1.21817, 1.63185, 1.18116, 0.46173, 2.22321, 1.24174, 3.06656, 2.1317, 3.25521, 2.01249, 3.21816, 3.24154, 2.21153, 2.06388, 2.94445, 2.0709, 1.70917, 1.36221, 2.17123, 1.81667, 0, 3.73142, 2.50661, 3.83585, 4.94251, 4.52711, 1.97208, 2.44531, 2.10607, 1.71736, 3.03484, 3.44852, 2.99783, 2.2784, 4.03988, 2.48968, 0.66486, 1.59972, 0.47621, 1.71893, 0.51326, 0.48988, 1.51989, 1.66754, 0.78697, 1.66052, 2.02225, 2.36921, 1.56019, 1.91475, 3.73142, 0, 1.22481, 0.10443, 1.21109, 0.79569, 1.75934, 1.28611, 1.62535, 2.01406, 0.69658, 0.2829, 0.73359, 1.45302, 0.30846, 1.26487, 0.55995, 0.37491, 0.7486, 0.49412, 0.71155, 0.73493, 0.29508, 0.44273, 0.43784, 0.43571, 0.79744, 1.1444, 0.33538, 0.68994, 2.50661, 1.22481, 0, 1.32924, 2.4359, 2.0205, 0.53453, 0.0613, 0.40054, 0.78925, 0.52823, 0.94191, 0.49122, 0.22821, 1.53327, 2.59411, 0.76929, 1.70415, 0.58064, 1.82336, 0.61769, 0.59431, 1.62432, 1.77197, 0.8914, 1.76495, 2.12668, 2.47364, 1.66462, 2.01918, 3.83585, 0.10443, 1.32924, 0, 1.10666, 0.69126, 1.86377, 1.39054, 1.72978, 2.11849, 0.80101, 0.38733, 0.83802, 1.55745, 0.20403, 3.70077, 1.87595, 2.81081, 1.6873, 2.93002, 1.72435, 1.70097, 2.73098, 2.87863, 1.99806, 2.87161, 3.23334, 3.5803, 2.77128, 3.12584, 4.94251, 1.21109, 2.4359, 1.10666, 0, 0.4154, 2.97043, 2.4972, 2.83644, 3.22515, 1.90767, 1.49399, 1.94468, 2.66411, 0.90263, 3.28537, 1.46055, 2.39541, 1.2719, 2.51462, 1.30895, 1.28557, 2.31558, 2.46323, 1.58266, 2.45621, 2.81794, 3.1649, 2.35588, 2.71044, 4.52711, 0.79569, 2.0205, 0.69126, 0.4154, 0, 2.55503, 2.0818, 2.42104, 2.80975, 1.49227, 1.07859, 1.52928, 2.24871, 0.48723, 0.73034, 1.09448, 0.15962, 1.28313, 0.04041, 1.24608, 1.26946, 0.23945, 0.0918, 0.97237, 0.09882, 0.26291, 0.60987, 0.19915, 0.15541, 1.97208, 1.75934, 0.53453, 1.86377, 2.97043, 2.55503, 0, 0.47323, 0.13399, 0.25472, 1.06276, 1.47644, 1.02575, 0.30632, 2.0678, 1.20357, 0.62125, 0.31361, 0.8099, 0.43282, 0.77285, 0.79623, 0.23378, 0.38143, 0.49914, 0.37441, 0.73614, 1.0831, 0.27408, 0.62864, 2.44531, 1.28611, 0.0613, 1.39054, 2.4972, 2.0818, 0.47323, 0, 0.33924, 0.72795, 0.58953, 1.00321, 0.55252, 0.16691, 1.59457, 0.86433, 0.96049, 0.02563, 1.14914, 0.09358, 1.11209, 1.13547, 0.10546, 0.04219, 0.83838, 0.03517, 0.3969, 0.74386, 0.06516, 0.2894, 2.10607, 1.62535, 0.40054, 1.72978, 2.83644, 2.42104, 0.13399, 0.33924, 0, 0.38871, 0.92877, 1.34245, 0.89176, 0.17233, 1.93381, 0.47562, 1.3492, 0.41434, 1.53785, 0.29513, 1.5008, 1.52418, 0.49417, 0.34652, 1.22709, 0.35354, 0.00819, 0.35515, 0.45387, 0.09931, 1.71736, 2.01406, 0.78925, 2.11849, 3.22515, 2.80975, 0.25472, 0.72795, 0.38871, 0, 1.31748, 1.73116, 1.28047, 0.56104, 2.32252, 1.7931, 0.03172, 0.90314, 0.22037, 1.02235, 0.18332, 0.2067, 0.82331, 0.97096, 0.09039, 0.96394, 1.32567, 1.67263, 0.86361, 1.21817, 3.03484, 0.69658, 0.52823, 0.80101, 1.90767, 1.49227, 1.06276, 0.58953, 0.92877, 1.31748, 0, 0.41368, 0.03701, 0.75644, 1.00504, 2.20678, 0.38196, 1.31682, 0.19331, 1.43603, 0.23036, 0.20698, 1.23699, 1.38464, 0.50407, 1.37762, 1.73935, 2.08631, 1.27729, 1.63185, 3.44852, 0.2829, 0.94191, 0.38733, 1.49399, 1.07859, 1.47644, 1.00321, 1.34245, 1.73116, 0.41368, 0, 0.45069, 1.17012, 0.59136, 1.75609, 0.06873, 0.86613, 0.25738, 0.98534, 0.22033, 0.24371, 0.7863, 0.93395, 0.05338, 0.92693, 1.28866, 1.63562, 0.8266, 1.18116, 2.99783, 0.73359, 0.49122, 0.83802, 1.94468, 1.52928, 1.02575, 0.55252, 0.89176, 1.28047, 0.03701, 0.45069, 0, 0.71943, 1.04205, 1.03666, 0.78816, 0.1467, 0.97681, 0.26591, 0.93976, 0.96314, 0.06687, 0.21452, 0.66605, 0.2075, 0.56923, 0.91619, 0.10717, 0.46173, 2.2784, 1.45302, 0.22821, 1.55745, 2.66411, 2.24871, 0.30632, 0.16691, 0.17233, 0.56104, 0.75644, 1.17012, 0.71943, 0, 1.76148, 2.79814, 0.97332, 1.90818, 0.78467, 2.02739, 0.82172, 0.79834, 1.82835, 1.976, 1.09543, 1.96898, 2.33071, 2.67767, 1.86865, 2.22321, 4.03988, 0.30846, 1.53327, 0.20403, 0.90263, 0.48723, 2.0678, 1.59457, 1.93381, 2.32252, 1.00504, 0.59136, 1.04205, 1.76148, 0 };
	double x[30] = { 0.55857, -0.88513, -1.03412, -1.40474, -1.10467, 0.30942, -0.20533, -0.46152, 1.5576, -2.32913, 1.75847, -0.72713, -0.66923, -1.23832, -0.33659, -1.12517, 0.36774, 0.61427, -1.99729, -0.80402, 1.26468, 0.8728, 0.41631, 1.10425, 0.71517, -0.3624, -0.0231, 0.04377, -1.19637, -1.60299 };
	double y[30] = { 1.26563, -0.55919, 0.37567, -0.74784, 0.49488, -0.71079, -0.73417, 0.29584, 0.44349, -0.43708, 0.43647, 0.7982, 1.14516, 0.33614, 0.6907, 2.50737, -1.22405, 0.00076, -1.32848, -2.43514, -2.01974, 0.53529, 0.06206, 0.4013, 0.79001, -0.52747, -0.94115, -0.49046, 0.22897, -1.53251 };
	int n = 30;
	int R = 1000, weight = 0, dst = 0, type = 1, nth = 2;
	
	// univariate case:
	bcov_test(ball_stat_value, p_value, x, y, &n, &R, &weight, &dst, &type, &nth);
	printf("Ball statistics: %f; ", ball_stat_value[0]);
	printf("p-value: %f \n", p_value[0]);
	printf("Ball statistics: %f; ", ball_stat_value[1]);
	printf("p-value: %f \n", p_value[1]);
	printf("Ball statistics: %f; ", ball_stat_value[2]);
	printf("p-value: %f \n", p_value[2]);
	test_value(ball_stat_value[0], 0.001256868);

	// multivariate case:
	dst = 1;
	bcov_test(ball_stat_value, p_value, x_dst, y_dst, &n, &R, &weight, &dst, &type, &nth);
	printf("Ball statistics: %f; ", ball_stat_value[0]);
	printf("p-value: %f \n", p_value[0]);
	printf("Ball statistics: %f; ", ball_stat_value[1]);
	printf("p-value: %f \n", p_value[1]);
	printf("Ball statistics: %f; ", ball_stat_value[2]);
	printf("p-value: %f \n", p_value[2]);
	test_value(ball_stat_value[0], 0.001256868);
	return;
}


void main()
{
	printf("-----------------------------\n");
	printf("Test ball statistic value: \n");
	printf("-----------------------------\n");
	test_bd_value();
	test_bcov_value();
	test_bcor_value();
	printf("-------------------------------- \n");
	printf("Test ball permutation procedure: \n");
	printf("-------------------------------- \n");
	bd_test_2sample_multithread_permutation();
	bd_test_ksample_multithread_permutation();
	bcov_test_multithread_permutaion();
	system("pause");
	return;
}


//void test_Find2()
//{
//	double ball_stat_value[2], p_value[2];
//	double x[100] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 2, 1, 0, 1, 2, 3, 4, 5, 6, 7, 3, 2, 1, 0, 1, 2, 3, 4, 5, 6, 4, 3, 2, 1, 0, 1, 2, 3, 4, 5, 5, 4, 3, 2, 1, 0, 1, 2, 3, 4, 6, 5, 4, 3, 2, 1, 0, 1, 2, 3, 7, 6, 5, 4, 3, 2, 1, 0, 1, 2, 8, 7, 6, 5, 4, 3, 2, 1, 0, 1, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0 };
//	int size[2] = { 5, 5 };
//	int k = 2;
//	int n = 10;
//	int dst = 1, R = 10, nth = 1;
//	int parallel_type = 1;
//	bd_test(ball_stat_value, p_value, x, size, &n, &k, &dst, &R, &nth);
//}


//int main()
//{
//	test_Find2();
//	system("pause");
//	return;
//}