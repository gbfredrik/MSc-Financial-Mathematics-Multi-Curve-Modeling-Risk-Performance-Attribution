#include "Distribution.h"

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>

#include <iostream>
#include <cmath>

using namespace boost::numeric;

Distribution::Distribution(ublas::matrix<double> series)
	//: time_series(series) {}
{
    ublas::matrix_column<ublas::matrix<double>> x(series, 0);
	time_series = x;
}

ublas::matrix<double> Distribution::calc_num_hessian(ublas::vector<double> const& /*x*/) {
    ublas::matrix<double> v(0, 0);

	return v;
}

ublas::vector<double> Distribution::calc_num_gradients(ublas::vector<double> const& x) {
	double increment{ 0.00001 };
    ublas::vector<double> num_gradients(4);
    ublas::vector<double> x_0diff(x);
	x_0diff(0) += increment;
    ublas::vector<double> x_1diff(x);
	x_1diff(1) += increment;

    double f_val{ function_value(x) };
	num_gradients(0) = (function_value(x_0diff) - f_val) / increment;
	num_gradients(1) = (function_value(x_1diff) - f_val) / increment;

	return num_gradients;
}

ublas::vector<double> Distribution::calc_gradients(ublas::vector<double> const& x) {
    ublas::vector<double> gradients(2);
	gradients(0) = 2 * (200 * pow(x(0), 3) - 200 * x(0) * x(1) + x(0) - 1);
	gradients(1) = 200 * (x(1) - pow(x(0), 2));

	return gradients;
}

double Distribution::function_value(ublas::vector<double> const& x) {
	return pow(1 - x(0), 2) + 100 * pow(x(1) - x(0) * x(0), 2);
}

double Distribution::calc_step_size(
    ublas::vector<double> const& x, 
    ublas::vector<double> const& d
) {
    double a{ 1.0 };
    //double c1{ pow(10, -4) };
    //double c2{ 0.9 };

	while (function_value(x + a * d) > function_value(x)) {
		a *= 0.5;
	}

    return a;
}

Distribution::~Distribution(void) {}
