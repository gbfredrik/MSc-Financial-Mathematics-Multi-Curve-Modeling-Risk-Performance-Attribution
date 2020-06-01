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

using namespace boost::numeric::ublas;

class T_Copula : public Distribution {
public:
	T_Copula(matrix<double> series);
	double function_value(matrix<double> const& x);
	vector<double> calcGradients(vector<double> const& x);
	vector<double> calcNumGradients(vector<double> const& x);
	void getSeries();
	double calcStepSize(vector<double> const& x, vector<double> const& d);

private:
	vector<double> m_GARCH_vec;
	matrix<double> time_series;
	vector<double> derivative_w(vector<double> const& x);
	vector<double> derivative_a(vector<double> const& x);
	vector<double> derivative_b(vector<double> const& x);

};