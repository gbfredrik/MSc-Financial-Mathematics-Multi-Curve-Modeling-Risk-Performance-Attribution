#ifndef RVSIM_H
#define RVSIM_H
#include "mex.h"
#include <boost/numeric/ublas/matrix.hpp>
#include <iostream>

class rvSim {
private:

public:
	static boost::numeric::ublas::matrix<double> gen_test(int rows, int cols);
	static boost::numeric::ublas::matrix<double> gen_normal(double m, double s, int k, int N);
	static double gen_gamma(double a);
	static double gen_normal(double m, double s);
	static double gen_uniform(double l, double u);
	static boost::numeric::ublas::matrix<double> gen_eps(boost::numeric::ublas::matrix<double> V, boost::numeric::ublas::vector<double> sigma, std::string type);
};


#endif
