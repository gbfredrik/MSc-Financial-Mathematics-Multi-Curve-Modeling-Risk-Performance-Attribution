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
	virtual vector<double> calcGradients(vector<double> x);
	virtual double function_value(vector<double>);
	virtual void getSeries();
	virtual double calcStepSize(vector<double> x, vector<double> d);
	~Distribution(void);
	

private:
	
};