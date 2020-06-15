#include "pch.h"
#include "Likelihood.h"

#include "../MathLibrary/matrixOperations.h"
#include "../MathLibrary/statisticsOperations.h"

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>

#include <cmath>

using namespace boost::numeric;

// Likelihood ratio test functions
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

// Likelihood ratio test functions with residual
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

// Calculate logarithm of function values in Likelihood Ratio Test
ublas::vector<double> Likelihood::d_i(
    ublas::vector<double> const& vector1,
    ublas::vector<double> const& vector2
) {
    return matrixOperations::vectorLog(vector1) - matrixOperations::vectorLog(vector2);
}

// Calculate d in Likelihood Ratio Test
double Likelihood::d(ublas::vector<double> const& d_i, int const N) {
    double d{ 0 };

    for (size_t i{ 0 }, n{ d_i.size() }; i < n; ++i) {
        d += d_i(i);
    }

    return d / N;
}

// Calculate sigma in Likelihood Ratio Test
double Likelihood::sigma(double const d, ublas::vector<double> const& d_i, int const N) {
    double sigma2{ 0 };

    //Calculate sigma^2
    for (size_t i{ 0 }, rows{ d_i.size() }; i < rows; ++i) {
        sigma2 += pow(d_i(i) - d, 2);
    }

    sigma2 /= (N - 1.0);

    return sqrt(sigma2);
}

// Calculate s in Likelihood Ratio Test
double Likelihood::s(double const sigma, int const N) {
    return sigma / sqrt(N);
}

// Calculate residual in Likelihood Ratio Test
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

// Calculate probability and test if it's better than cdf of confidence level in Likelihood Ratio Test
int Likelihood::is_better(double const d, double const s, double const confidence_level) {
    double cdfInv{ statisticsOperations::invCDFNorm(confidence_level) };
    double p{ d / s };

    return p > cdfInv;
}
