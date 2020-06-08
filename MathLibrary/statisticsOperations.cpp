#include "pch.h"
#include "statisticsOperations.h"

#include "matrixOperations.h"

//#include "mex.h"
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/math/distributions/normal.hpp>
#include <boost/math/distributions/students_t.hpp>
#include <boost/math/distributions/chi_squared.hpp>

#include <numeric>

using namespace boost::numeric;

// Calculates the covariance matrix
ublas::matrix<double> statisticsOperations::covm(ublas::matrix<double> const& input) {
	//size_t m{ input.size1() };
	//size_t n{ input.size2() };
	//ublas::matrix<double> cov(n, n);
	//ublas::matrix<double> A(m, n);
	//double mean{ 0.0 };

	//for (size_t j{ 0 }; j < n; ++j) {
	//	mean = vectorMean(column(input, j));
	//	for (size_t i{ 0 }; i < m; ++i) {
	//		A(i, j) = input(i, j) - mean;
	//	}
	//}

	//cov = prod(trans(A), A) / (m - 1);

	//return cov;
	ublas::matrix<double> m_centered{ matrixOperations::center_matrix(input) };

	return prod(trans(m_centered), m_centered) / (m_centered.size1() - 1);
}

// Todo: Remove
//double statisticsOperations::vectorMean(ublas::vector<double> const& input) {
//	double sum{ std::accumulate(input.begin(), input.end(), 0.0) };
//	double mean{ sum / input.size() };
//
//	return mean;
//}

// Calculate the Pearson correlation matrix, TODO: reduce the amount of calls to pearson_rho()
ublas::matrix<double> statisticsOperations::corrm(ublas::matrix<double> const& input) {
	size_t m{ input.size1() };
	size_t n{ input.size2() };
	ublas::matrix<double> corr(n, n);
	ublas::vector<double> X(m);
	ublas::vector<double> Y(m);

	for (size_t i{ 0 }; i < n; ++i) {
		for (size_t j{ 0 }; j < n; ++j) {
			if (i == j) {
				corr(i, j) = 1;
			} else {
				corr(i, j) = pearson_rho(column(input, i), column(input, j));
			}
		}
	}

	return corr;

}

// Calculate the Pearson correlation coefficient
double statisticsOperations::pearson_rho(
	ublas::vector<double> const& X, 
	ublas::vector<double> const& Y
) {
	double rho{ 0 };
	size_t m{ X.size() };
	double numerator{ 0 };
	double denomenator_a{ 0 };
	double denomenator_b{ 0 };

	//Calculate mean of input vectors
	double X_hat{ matrixOperations::vector_average(X) }; // Changed from vectorMean
	double Y_hat{ matrixOperations::vector_average(Y) }; // Changed from vectorMean

	for (size_t i{ 0 }; i < m; ++i) {
		numerator += (X(i) - X_hat) * (Y(i) - Y_hat);
		denomenator_a += pow(X(i) - X_hat, 2);
		denomenator_b += pow(Y(i) - Y_hat, 2);
	}

	return numerator / (sqrt(denomenator_a * denomenator_b));
}

// Calculates the first garch volatility values with the full dataset
// Check if the GJR-term is needed
ublas::vector<double> statisticsOperations::GARCH(
	ublas::vector<double> const& omega,
	ublas::vector<double> const& alpha,
	ublas::vector<double> const& beta,
	ublas::vector<double> const& gamma,
	ublas::matrix<double> const& E,
	ublas::matrix<double> const& fHist
) {
	size_t m{ fHist.size1() }; // Number of days in fHist
	size_t k{ E.size2() }; // Number of risk factors

	ublas::vector<double> dXi(k);
	ublas::vector<double> sigmaPrevSq(k);
	ublas::vector<double> sigmaSq(k);
	ublas::vector<double> sigma(k);

	dXi = prod(trans(E), trans(row(fHist, 1) - row(fHist, 0)));

	for (size_t i{ 0 }; i < k; ++i) {
		sigmaPrevSq(i) = omega(i) + alpha(i) * pow(dXi(i), 2) + beta(i) * pow(dXi(i), 2);
	}
	
	for (size_t i{ 3 }; i < m; ++i) {
		dXi = prod(trans(E), trans(row(fHist, i - 1) - row(fHist, i - 2)));
		
		for (size_t j{ 0 }; j < k; ++j) {
			sigmaSq(j) = omega(j) + alpha(j) * pow(dXi(j), 2) + beta(j) * sigmaPrevSq(j);
			sigmaPrevSq(j) = sigmaSq(j);
		}
	}

	for (size_t i{ 0 }; i < k; ++i) {
		sigma(i) = sqrt(sigmaSq(i));
	}

	return sigma;
}

// Calculates the updated garch volatility
// Check if the GJR-term is needed
ublas::vector<double> statisticsOperations::GARCH(
	ublas::vector<double> const& omega,
	ublas::vector<double> const& alpha,
	ublas::vector<double> const& beta,
	ublas::vector<double> const& gamma,
	ublas::matrix<double> const& E,
	ublas::vector<double> const& fPrev,
	ublas::vector<double> const& fPrevPrev,
	ublas::vector<double> const& sigmaPrev
) {
	size_t k{ E.size2() }; // Number of risk factors

	ublas::vector<double> dXi(k);
	ublas::vector<double> sigmaSq(k);
	ublas::vector<double> sigma(k);

	dXi = prod(trans(E), trans(fPrev - fPrevPrev));

	for (size_t i{ 0 }; i < k; ++i) {
		sigmaSq(i) = omega(i) + alpha(i) * pow(dXi(i), 2) + beta(i) * pow(sigmaPrev(i), 2);
		sigma(i) = sqrt(sigmaSq(i));
	}

	return sigma;
}

double statisticsOperations::invCDFNorm(double const u, double const mu, double const sigma) {
	boost::math::normal norm(mu, sigma);

	return quantile(norm, u);
}

double statisticsOperations::invCDFT(double const u, double const df) {
	boost::math::students_t t(df);

	return quantile(t, u);
}

double statisticsOperations::invCDFchi2(double const u, double const df) {
    boost::math::chi_squared chi(df);

    return quantile(chi, u);
}
