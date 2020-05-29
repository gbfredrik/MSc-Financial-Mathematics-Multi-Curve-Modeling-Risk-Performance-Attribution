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

class Student_t : public Distribution {
public:
	Student_t(vector<double> series);
	vector<double> create_GARCH_vec(vector<double> x);
	double function_value(vector<double> x);
	vector<double> calcGradients(vector<double> x);
	void getSeries();
	double calcStepSize(vector<double> x, vector<double> d);

private:
	vector<double> GARCH_vec;

	vector<double> derivative_w(vector<double> x, vector<double> GARCH_vec);
	vector<double> derivative_a(vector<double> x, vector<double> GARCH_vec);
	vector<double> derivative_b(vector<double> x, vector<double> GARCH_vec);

};