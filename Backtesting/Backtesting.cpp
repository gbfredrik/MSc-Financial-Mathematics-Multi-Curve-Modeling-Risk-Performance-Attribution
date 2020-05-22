// Backtesting.cpp : Defines the functions for the static library.
//

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
#include <boost/range/counting_range.hpp>
#include <unsupported/Eigen/MatrixFunctions>
#include <boost/math/distributions/normal.hpp>

using namespace boost::numeric::ublas;
using namespace boost::math;

// Pricing of FX swaps somewhere else


//Likelihood ratio test functions
vector<double> Likelihood::likelihoodRatioTest(matrix<double> const& matrix1, matrix<double> const& matrix2)
{
	size_t rows = matrix1.size1();
	size_t columns = matrix1.size2();

	matrix<double> d_i(rows, columns);
	vector<double> d(rows);
	vector<double> sigma(rows);
	vector<double> s(rows);
	vector<double> p(rows);
	double N = columns;
	vector<double> confidence_level(rows);

	d_i = Likelihood::d_i(matrix1, matrix2);
	d = Likelihood::d(d_i, N);
	sigma = Likelihood::sigma(d, d_i, N);
	s = Likelihood::s(sigma, N);
	p = Likelihood::probability(d, s, confidence_level);
	return p;
}

//Likelihood ratio test functions with residual
vector<double> Likelihood::likelihoodRatioTestResidual(matrix<double> const& matrix1, matrix<double> const& matrix2, matrix<double> const& prices1, matrix<double> const& prices2)
{
	size_t rows = matrix1.size1();
	size_t columns = matrix1.size2();

	matrix<double> d_i(rows, columns);
	vector<double> d(rows);
	vector<double> sigma(rows);
	vector<double> s(rows);
	vector<double> p(rows);
	double N = columns;
	vector<double> confidence_level(rows);

	d_i = Likelihood::dResidual_i(matrix1, matrix2, prices1, prices2);
	d = Likelihood::d(d_i, N);
	sigma = Likelihood::sigma(d, d_i, N);
	s = Likelihood::s(sigma, N);
	p = Likelihood::probability(d, s, confidence_level);
	return p;
}

//Calculate logarithm of function values in Likelihood Ratio Test
matrix<double> Likelihood::d_i(matrix<double> const& matrix1, matrix<double> const& matrix2)
{
	size_t rows = matrix1.size1();
	size_t columns = matrix1.size2();

	matrix<double> matrix1_log(rows, columns);
	matrix<double> matrix2_log(rows, columns);
	matrix<double> matrix_log(rows,columns);

	matrix1_log = matrixOperations::matrixLog(matrix1);
	matrix2_log = matrixOperations::matrixLog(matrix2);
	return matrix_log = matrix1_log - matrix2_log;
}

//Calculate d in Likelihood Ratio Test
vector<double> Likelihood::d(matrix<double> const& d_i, double const N)
{
	size_t rows = d_i.size1();
	size_t columns = d_i.size2();
	vector<double> d(rows,0);
	for (size_t i = 0; i < rows; ++i) {
		for (size_t j = 0; j < columns; ++j) {
			d(i) += d_i(i, j);
		}
	}
	return d/N;
}


//Calculate sigma in Likelihood Ratio Test
vector<double> Likelihood::sigma(vector<double> const& d, matrix<double> const& d_i, double const N)
{
	size_t rows = d_i.size1();
	size_t columns = d_i.size2();
	vector<double> sigma(rows, 0);

	//Calculate sigma^2
	for (size_t i = 0; i < rows; ++i) {
		for (size_t j = 0; j < columns; ++j) {
			sigma(i) += pow(d_i(i, j)-d(i),2);
		}
	}
	sigma /= (N - 1);

	//Calculate sigma
	for (size_t i = 0; i < rows; ++i) {
		sigma(i) = sqrt(sigma(i));
	}
	return sigma;
}

//Calculate s in Likelihood Ratio Test
vector<double> Likelihood::s(vector<double> const& sigma, double const N)
{
	size_t length = sigma.size();
	vector<double> s(length);
	s = sigma / sqrt(N);
	return s;
}

//Calculate residual in Likelihood Ratio Test
matrix<double> Likelihood::dResidual_i(matrix<double> const& values_functions1, matrix<double> const& values_functions2, matrix<double> const& prices1, matrix<double> const& prices2)
{
	size_t rows = values_functions1.size1();
	size_t columns = values_functions1.size2();
	matrix<double> residual1(rows, columns);
	matrix<double> residual2(rows, columns);
	matrix<double> dResidual_i(rows, columns);

	residual1 = prices1 - values_functions1;
	residual2 = prices2 - values_functions2;
	dResidual_i = element_prod(residual1, residual1) - element_prod(residual2, residual2);
	return dResidual_i;
}

//Calculate probability in Likelihood Ratio Test
vector<double> Likelihood::probability(vector<double> const& d, vector<double> const& s, vector<double> const& confidence_level)
{
	size_t length = d.size();
	vector<double> p(length);
	normal norm;
	p = element_div(d, s);
	//for (size_t i = 0; i < length; ++i) { 
		//p(i) = quantile(norm, confidence_level(i));
	//}
	return p;
}






