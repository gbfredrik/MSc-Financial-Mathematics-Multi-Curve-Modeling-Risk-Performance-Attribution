// Backtesting.cpp : Defines the functions for the static library.
#include "pch.h"
#include "backtesting.h"

#define _USE_MATH_DEFINES

#include "../MathLibrary/matrixOperations.h"
#include "../MathLibrary/statisticsOperations.h"

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/math/distributions/normal.hpp>

#include <numeric>
#include <cmath>

using namespace boost::numeric;

//Likelihood ratio test functions
int Likelihood::likelihood_ratio_test(
    ublas::vector<double> const& vector1, 
    ublas::vector<double> const& vector2, 
    double const confidence_level
) {
    size_t N{ vector1.size() };
    ublas::vector<double> d_i{ Likelihood::d_i(vector1, vector2) };
    double d{ Likelihood::d(d_i, N) };
    double sigma{ Likelihood::sigma(d, d_i, N) };
    double s{ Likelihood::s(sigma, N) };

    return Likelihood::is_better(d, s, confidence_level);
}

//Likelihood ratio test functions with residual
int Likelihood::likelihood_ratio_test_residual(
    ublas::vector<double> const& vector1, 
    ublas::vector<double> const& vector2, 
    ublas::vector<double> const& prices1, 
    ublas::vector<double> const& prices2, 
    double const confidence_level
) {
    size_t N{ vector1.size() };
    ublas::vector<double> d_i{ Likelihood::dResidual_i(vector1, vector2, prices1, prices2) };
    double d{ Likelihood::d(d_i, N) };
    double sigma{ Likelihood::sigma(d, d_i, N) };
    double s{ Likelihood::s(sigma, N) };

    return Likelihood::is_better(d, s, confidence_level);
}

//Calculate logarithm of function values in Likelihood Ratio Test
ublas::vector<double> Likelihood::d_i(
    ublas::vector<double> const& vector1, 
    ublas::vector<double> const& vector2
) {
    return matrixOperations::vectorLog(vector1) - matrixOperations::vectorLog(vector2);
}

//Calculate d in Likelihood Ratio Test
double Likelihood::d(ublas::vector<double> const& d_i, int const N) {
    double d{ 0 };

    for (size_t i{ 0 }, n{ d_i.size() }; i < n; ++i) {
        d += d_i(i);
    }

    return d / N;
}

//Calculate sigma in Likelihood Ratio Test
double Likelihood::sigma(double const d, ublas::vector<double> const& d_i, int const N) {
    double sigma2{ 0 };

    //Calculate sigma^2
    for (size_t i{ 0 }, rows{ d_i.size() }; i < rows; ++i) {
        sigma2 += pow(d_i(i)-d,2);
    }

    sigma2 /= (N - 1.0);

    return sqrt(sigma2);
}

//Calculate s in Likelihood Ratio Test
double Likelihood::s(double const sigma, int const N) {
    return sigma / sqrt(N);
}

//Calculate residual in Likelihood Ratio Test
ublas::vector<double> Likelihood::dResidual_i(
    ublas::vector<double> const& values_functions1, 
    ublas::vector<double> const& values_functions2, 
    ublas::vector<double> const& prices1, 
    ublas::vector<double> const& prices2
) {
    ublas::vector<double> residual1(prices1 - values_functions1);
    ublas::vector<double> residual2(prices2 - values_functions2);

    return element_prod(residual1, residual1) - element_prod(residual2, residual2);
}

//Calculate probability and test if it's better than cdf of confidence level in Likelihood Ratio Test
int Likelihood::is_better(double const d, double const s, double const confidence_level) {
    double cdfInv{ statisticsOperations::invCDFNorm(confidence_level) };
    double p{ d / s };

    return p > cdfInv;
}

ublas::vector<double> KernelDensity::kde_multi(
    ublas::matrix<double> const& x_simulated,
    ublas::vector<double> const& x_realized
) {
    size_t m{ x_realized.size() };
    ublas::vector<double> f(m);

    for (size_t j{ 0 }; j < m; ++j) {
        f(j) = KernelDensity::kde(column(x_simulated, j), x_realized(j));
    }

    return f;
}

double KernelDensity::kde(
    ublas::vector<double> const& x_simulated, 
    double const x_realized
) {
    size_t n{ x_simulated.size() };
    double variance{ matrixOperations::vector_variance(x_simulated) };
    double sigma{ sqrt(variance) };
    double h{ pow(4.0 / (3.0 * n), 1.0 / 5.0) * sigma };
    double sum{ 0 };
    
    for (size_t i{ 0 }; i < n; ++i) {
        sum += std::exp(-pow(x_realized - x_simulated(i), 2) / (2 * pow(h, 2)));
    }
    
    return (1.0 / (sqrt(2.0 * M_PI) * h * n)) * sum;
}
