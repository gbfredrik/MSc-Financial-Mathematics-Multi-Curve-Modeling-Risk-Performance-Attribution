#pragma once

#include "mex.h"
#include <boost/numeric/ublas/matrix.hpp>

#include <iostream>

class rvSim {
public:
	static boost::numeric::ublas::matrix<double> gen_test(int rows, int cols);
	static boost::numeric::ublas::matrix<double> gen_normal(double m, double s, int k, int N);
	static double gen_gamma(double a);
	static double gen_normal(double m, double s);
	static double gen_uniform(double l, double u);
	static boost::numeric::ublas::matrix<double> genEps(
		boost::numeric::ublas::matrix<double> const& V, 
		boost::numeric::ublas::matrix<double> const& E, 
		boost::numeric::ublas::vector<double> const& sigma, 
		std::string const& type
	);
};
