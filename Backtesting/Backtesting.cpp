// Backtesting.cpp : Defines the functions for the static library.

#include "pch.h"
#include "framework.h"
#include "backtesting.h"
#include "../MathLibrary/matrixOperations.h"

#include <iostream>
#include <numeric>
#include <cmath>

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/bind.hpp>
#include <boost/ref.hpp >
#include <boost/range/counting_range.hpp>
#include <unsupported/Eigen/MatrixFunctions>
#include <boost/math/distributions/normal.hpp>
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/variance.hpp>

using namespace boost::numeric;
using namespace boost::math;


// Pricing of FX swaps somewhere else


//Likelihood ratio test functions
int Likelihood::likelihoodRatioTest(ublas::vector<double> const& vector1, ublas::vector<double> const& vector2, double confidence_level)
{
	size_t N = vector1.size();

	ublas::vector<double> d_i(N);
	double d;
	double sigma;
	double s;
	int isBetter;

	d_i = Likelihood::d_i(vector1, vector2);
	d = Likelihood::d(d_i, N);
	sigma = Likelihood::sigma(d, d_i, N);
	s = Likelihood::s(sigma, N);
	isBetter = Likelihood::isBetter(d, s, confidence_level);
	return isBetter;
}

//Likelihood ratio test functions with residual
int Likelihood::likelihoodRatioTestResidual(ublas::vector<double> const& vector1, ublas::vector<double> const& vector2, ublas::vector<double> const& prices1, ublas::vector<double> const& prices2, double confidence_level)
{
	int N = vector1.size();
	//size_t columns = vector1.size2();

	ublas::vector<double> d_i;
	double d;
	double sigma;
	double s;
	int isBetter;
	//ublas::vector<double> confidence_level(rows);

	d_i = Likelihood::dResidual_i(vector1, vector2, prices1, prices2);
	d = Likelihood::d(d_i, N);
	sigma = Likelihood::sigma(d, d_i, N);
	s = Likelihood::s(sigma, N);
	isBetter = Likelihood::isBetter(d, s, confidence_level);
	return isBetter;
}

//Calculate logarithm of function values in Likelihood Ratio Test
ublas::vector<double> Likelihood::d_i(ublas::vector<double> const& vector1, ublas::vector<double> const& vector2)
{
	size_t rows = vector1.size();

	ublas::vector<double> vector1_log(rows);
	ublas::vector<double> vector2_log(rows);
	ublas::vector<double> vector_log(rows);

	vector1_log = matrixOperations::vectorLog(vector1);
	vector2_log = matrixOperations::vectorLog(vector2);
	return vector_log = vector1_log - vector2_log;
}

//Calculate d in Likelihood Ratio Test
double Likelihood::d(ublas::vector<double> const& d_i, int N)
{
	size_t rows = d_i.size();

	double d = 0;
	for (size_t i = 0; i < rows; ++i) {
		d += d_i(i);
	}
	return d/N;
}


//Calculate sigma in Likelihood Ratio Test
double Likelihood::sigma(double d, ublas::vector<double> const& d_i, int N)
{
	size_t rows = d_i.size();
	double sigma{0};

	//Calculate sigma^2
	for (size_t i = 0; i < rows; ++i) {
		sigma += pow(d_i(i)-d,2);
	}
	sigma /= (N - 1);
	return sqrt(sigma);
}

//Calculate s in Likelihood Ratio Test
double Likelihood::s(double sigma, int N)
{
	return sigma / sqrt(N);
}

//Calculate residual in Likelihood Ratio Test
ublas::vector<double> Likelihood::dResidual_i(ublas::vector<double> const& values_functions1, ublas::vector<double> const& values_functions2, ublas::vector<double> const& prices1, ublas::vector<double> const& prices2)
{
	size_t rows = values_functions1.size();
	ublas::vector<double> residual1(rows);
	ublas::vector<double> residual2(rows);
	ublas::vector<double> dResidual_i(rows);

	residual1 = prices1 - values_functions1;
	residual2 = prices2 - values_functions2;
	dResidual_i = element_prod(trans(residual1), residual1) - element_prod(trans(residual2), residual2);
	return dResidual_i;
}

//Calculate probability and test if it's better than cdf of confidence level in Likelihood Ratio Test
int Likelihood::isBetter(double d, double s, double confidence_level)
{
	double p;
	int isBetter;
	normal norm;
	double cdfInv = quantile(norm, confidence_level);
	p = d / s;

	if (p > cdfInv) {
		isBetter = 1;
	}
	else {
		isBetter = 0;
	}
		
	std::cout << "p: " << p << std::endl;
	return isBetter;
}

ublas::vector<double> KernelDensity::kde_multi(ublas::matrix<double> x_simulated, ublas::vector<double> x_realized) {
	size_t m = x_realized.size();
	ublas::vector<double> f(m);

	for (size_t j = 0; j < m; ++j) {
		ublas::matrix_column<ublas::matrix<double>> x_simulated_col(x_simulated, j); //Row or column?
		double x_realized_m = x_realized(j);
		f(j) = KernelDensity::kde(x_simulated_col, x_realized_m);
	}
	return f;
}

double KernelDensity::kde(ublas::vector<double> const& x_simulated, double x_realized) {
	using namespace boost::accumulators;

	int n = x_simulated.size();
	double f;
	double variance;
	double sigma;
	double h;
	double sum = 0;
	const double pi = 3.141592653589793238L;
	
	variance = matrixOperations::vector_variance(x_simulated);
	sigma = sqrt(variance);
	h = pow(4.0 / (3.0* n), 1.0/5.0) * sigma;
	
	for (int i = 0; i < n; ++i){
		sum += std::exp(- pow(x_realized - x_simulated(i),2)) / (2 * pow(h, 2));
	}
	f = (1.0 / (sqrt(2.0 * pi) * h * n))*sum;
	return f;
}






