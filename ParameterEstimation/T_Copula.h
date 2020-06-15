#pragma once
#include "Distribution.h"
#include "../MathLibrary/matrixOperations.h"
#include <boost/math/distributions/students_t.hpp>

#include <Eigen/Dense>

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
#include <boost/qvm/mat_operations.hpp>

#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/variance.hpp>
#include <boost/bind.hpp>

#include <boost/math/special_functions/digamma.hpp>

using namespace boost::numeric::ublas;

class T_Copula : public Distribution {
public:
	T_Copula(matrix<double> series);
	double function_value(vector<double> const& x);
	vector<double> calcGradients(vector<double> const& x);
	vector<double> calcNumGradients(vector<double> const& x);
	void getSeries();
	double calcStepSize(vector<double> const& x, vector<double> const& d);

//private:
	matrix<double> time_series;
	matrix<double> buildP(vector<double> const& x);
	vector<double> matrixToVector(matrix<double> const& matrix);
	matrix<double> vectorToMatrix(vector<double> const& vec);
	vector<double> getElements(matrix<double> const& matrix);
	vector<double> kronOfVectors(vector<double> const& v1, vector<double> const& v2);
	double dGamma(double t);

};