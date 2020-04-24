#pragma once

#include <boost/numeric/ublas/matrix.hpp>


class rvSim {
private:

public:
	static boost::numeric::ublas::matrix<double> gen_test(int rows, int cols);
	static boost::numeric::ublas::vector<double> gen_normal(double m, double s, int N);
	static double gen_gamma(double a);
	static double gen_normal(double m, double s);
	static double gen_uniform(double l, double u);
};
