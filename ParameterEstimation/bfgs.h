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

#include <boost/random.hpp>
#include <boost/random/normal_distribution.hpp>


//Rosenbrock function: f(x,y) = (1-x)^2 + 100(y-x^2)^2 
using namespace boost::numeric::ublas;

class bfgs {
public:
	static vector<double> minimize(vector<double> start, matrix<double> H_inv, int max_iter, float epsilon, Distribution* dist);

	static Distribution dist;
	static int function_type;
	static vector<double> calcGradients(vector<double> x);
	static double calcStepSize(vector<double> x, vector<double> d, Distribution* dist);
	static double f(vector<double> x);
	//static double rosenbrock(vector<double> x);
};