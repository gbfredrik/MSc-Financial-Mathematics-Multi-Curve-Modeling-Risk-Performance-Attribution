#pragma once
#include "Distribution.h"

#include <iostream>
#include <numeric>
#include <cmath>

//Boost packages for numeric represenation
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/io.hpp>

//Boost packages for statistics
#include <boost/math/tools/bivariate_statistics.hpp>
#include <boost/math/distributions/normal.hpp>

#include <boost/math/constants/constants.hpp>

using namespace boost::numeric::ublas;

class Gaussian : public Distribution {
public: 
	Gaussian(vector<double> series);
	vector<double> create_GARCH_vec(vector<double> x);
	double function_value(vector<double> x);
	vector<double> calcGradients(vector<double> x);
	void getSeries();

private:
	vector<double> GARCH_vec;
};