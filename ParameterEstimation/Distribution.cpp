#include "Distribution.h"



Distribution::Distribution(vector<double> series):
	time_series(series) {}

vector<double> Distribution::calcGradients(vector<double> x) {

	double dx = 2 * (200 * pow(x(0), 3) - 200 * x(0) * x(1) + x(0) - 1);
	double dy = 200 * (x(1) - pow(x(0), 2));

	vector<double> gradients(2);
	gradients(0) = dx;
	gradients(1) = dy;
	return gradients;
}


double Distribution::function_value(vector<double> x) {

	double function_value = pow(1 - x(0), 2) + 100 * pow(x(1) - x(0) * x(0), 2);
	return function_value;
}


Distribution::~Distribution(void) {}