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


#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/variance.hpp>
#include <boost/bind.hpp>


#include <boost/qvm/mat_operations.hpp>


using namespace boost::numeric::ublas;

class Gaussian : public Distribution {
public: 
	Gaussian(vector<double> series);
	void update_GARCH_vec(vector<double> x);
	double function_value(vector<double> x);
	matrix<double> calcNumHessian(vector<double> x);
	vector<double> calcNumGradients(vector<double> x);
	vector<double> calcGradients(vector<double> x);
	void getSeries();
	double calcStepSize(vector<double> x, vector<double> d);

private:
	vector<double> m_GARCH_vec;
	double garch0;
	vector<double> derivative_w(vector<double> x);
	vector<double> derivative_a(vector<double> x);
	vector<double> derivative_b(vector<double> x);

};