#pragma once

#include "mex.h"
#include <boost/numeric/ublas/matrix.hpp>

#include <iostream>

extern boost::numeric::ublas::matrix<double> hehe;

class rvSim {
public:
	static boost::numeric::ublas::matrix<double> gen_test(int rows, int cols);
<<<<<<< HEAD
	static boost::numeric::ublas::matrix<double> gen_normal(double m, double s, size_t k, size_t N);
	static double gen_gamma(double df);
	static boost::numeric::ublas::matrix<double> genEps(boost::numeric::ublas::matrix<double> V, 
		boost::numeric::ublas::vector<double> mu, boost::numeric::ublas::vector<double> sigma, std::string type,
		boost::numeric::ublas::vector<double> dfM);
	static boost::numeric::ublas::vector<double> genEps(boost::numeric::ublas::vector<double> V,
		boost::numeric::ublas::vector<double> mu, boost::numeric::ublas::vector<double> sigma, std::string type,
		boost::numeric::ublas::vector<double> dfM);
=======
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
>>>>>>> master
};
