#pragma once


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

using namespace boost::numeric::ublas;

class Distribution {
public:
	Distribution(vector<double> time_series);
	vector<double> time_series;
	virtual matrix<double> calcNumHessian(vector<double> const& x);
	virtual vector<double> calcNumGradients(vector<double> const& x);
	virtual vector<double> calcGradients(vector<double> const& x);
	virtual double function_value(vector<double> const& x);
	virtual void getSeries();
	virtual double calcStepSize(vector<double> const& x, vector<double> const& d);
	~Distribution(void);
	

private:
	
};