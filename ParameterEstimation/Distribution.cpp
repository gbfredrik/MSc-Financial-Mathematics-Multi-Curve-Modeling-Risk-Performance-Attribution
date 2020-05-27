#include "Distribution.h"



Distribution::Distribution(vector<double> series):
	time_series(series) {}

vector<double> Distribution::calcNumGradients(vector<double> const& x) {



	double increment{ 0.00001 };
	vector<double> num_gradients(4);
	vector<double> x_0diff(x);
	x_0diff(0) += increment;
	vector<double> x_1diff(x);
	x_1diff(1) += increment;

	num_gradients(0) = (function_value(x_0diff) - function_value(x)) / increment;
	num_gradients(1) = (function_value(x_1diff) - function_value(x)) / increment;

	return num_gradients;
}

matrix<double> Distribution::calcNumHessian(vector<double> const& x) {
	matrix<double> v(0,0);
	return v;
}



vector<double> Distribution::calcGradients(vector<double> const& x) {

	std::cout << "calcGradients in Distribution \n";

	double dx = 2 * (200 * pow(x(0), 3) - 200 * x(0) * x(1) + x(0) - 1);
	double dy = 200 * (x(1) - pow(x(0), 2));

	vector<double> gradients(2);
	gradients(0) = dx;
	gradients(1) = dy;
	return gradients;
}

void Distribution::getSeries() {
	
	std::cout << "In distribution: " << time_series;
}


double Distribution::function_value(vector<double> const& x) {

	double function_value = pow(1 - x(0), 2) + 100 * pow(x(1) - x(0) * x(0), 2);
	return function_value;
}

double Distribution::calcStepSize(vector<double> const& x, vector<double> const& d) {

	double a = 1;
	double c1 = pow(10, -4);
	double c2 = 0.9;


	//while (dist->function_value(x + a * d) > dist->function_value(x) + c1 * a * inner_prod(dist->calcGradients(x), d)
	while (function_value(x + a * d) > function_value(x))
	{
		std::cout << "new f : " << function_value(x + a * d) << "\n";
		std::cout << "old f : " << function_value(x) << "\n";
		std::cout << "steglängd = " << a << "\n";
		std::cout << "new parameters2 = " << x(0) + a * d(0) << ", " << x(1) + a * d(1) << "\n";
		a = a * 0.5;
	}

	std::cout << "new f : " << function_value(x + a * d) << "\n";
	std::cout << "old f : " << function_value(x) << "\n";
	std::cout << "steglängd = " << a << "\n \n";

	return a;
}

Distribution::~Distribution(void) {}