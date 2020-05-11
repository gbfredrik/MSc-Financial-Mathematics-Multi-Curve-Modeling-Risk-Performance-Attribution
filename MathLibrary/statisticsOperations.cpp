#include "pch.h"
#include "statisticsOperations.h"
#include "mex.h"
#include <numeric>

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/math/distributions/normal.hpp>
#include <boost/math/distributions/students_t.hpp>
using namespace boost::numeric::ublas;


// Calculates the covariance matrix
matrix<double> statisticsOperations::covm(matrix<double> const& input) {

	size_t m = input.size1();
	size_t n = input.size2();

	matrix<double> cov(n, n);
	matrix<double> A(m, n);
	double mean = 0.0;

	for (size_t j = 0; j < n; j++) {
		mean = vectorMean(column(input, j));
		for (size_t i = 0; i < m; i++) {
			A(i, j) = input(i, j) - mean;
		}
	}

	cov = prod(trans(A), A) / (m - 1);

	return cov;
}

double statisticsOperations::vectorMean(vector<double> const& input) {

	double sum = std::accumulate(input.begin(), input.end(), 0.0);
	double mean = sum / input.size();

	return mean;
}


// Calculate the Pearson correlation matrix, TODO: reduce the amount of calls to pearson_rho()
matrix<double> statisticsOperations::corrm(matrix<double> const& input) {
	size_t m = input.size1();
	size_t n = input.size2();
	matrix<double> corr(n, n);
	vector<double> X(m);
	vector<double> Y(m);

	for (size_t i = 0; i < n; i++) {
		for (size_t j = 0; j < n; j++) {
			if (i == j) {
				corr(i, j) = 1;
			}
			else {
				corr(i, j) = pearson_rho(column(input, i), column(input, j));
			}
		}
	}
	return corr;
}
// Calculate the Pearson correlation coefficient
double statisticsOperations::pearson_rho(vector<double> const& X, vector<double> const& Y) {
	double rho = 0.0;
	size_t m = X.size();
	double numerator = 0;
	double denomenator_a = 0;
	double denomenator_b = 0;

	//Calculate mean of input vectors
	double X_hat = vectorMean(X);
	double Y_hat = vectorMean(Y);

	for (size_t i = 0; i < m; i++) {
		numerator = numerator + (X(i) - X_hat) * (Y(i) - Y_hat);
		denomenator_a = denomenator_a + pow(X(i) - X_hat, 2);
		denomenator_b = denomenator_b + pow(Y(i) - Y_hat, 2);
	}

	return numerator / (sqrt(denomenator_a) * sqrt(denomenator_b));
}

// Calculates the first garch volatility values with the full dataset
// Check if the GJR-term is needed
vector<double> statisticsOperations::GARCH(vector<double> omega, vector<double> alpha, vector<double> beta, matrix<double> E, matrix<double> fHist) {

	size_t m = fHist.size1(); // Number of days in fHist
	size_t n = fHist.size2(); // Number of discretization points on the curves
	size_t k = E.size2(); // Number of risk factors

	vector<double> dXi(k);
	vector<double> sigmaPrevSq(k);
	vector<double> sigmaSq(k);
	vector<double> sigma(k);

	dXi = prod(trans(E), trans(row(fHist, 1) - row(fHist, 0)));

	
	for (size_t i = 0; i < k; i++) {
		sigmaPrevSq(i) = omega(i) + alpha(i) * pow(dXi(i), 2) + beta(i) * pow(dXi(i), 2);
	}
	
	for (size_t i = 3; i < m; i++) {
		
		dXi = prod(trans(E), trans(row(fHist, i - 1) - row(fHist, i - 2)));
		
		for (size_t j = 0; j < k; j++) {
			sigmaSq(j) = omega(j) + alpha(j) * pow(dXi(j), 2) + beta(j) * sigmaPrevSq(j);
			sigmaPrevSq(j) = sigmaSq(j);
		}
	}

	for (size_t i = 0; i < k; i++) {
		sigma(i) = sqrt(sigmaSq(i));
	}
	

	return sigma;
}

// Calculates the updated garch volatility
// Check if the GJR-term is needed
vector<double> statisticsOperations::GARCH(vector<double> omega, vector<double> alpha, vector<double> beta, matrix<double> E, vector<double> ft1, vector<double> ft2, vector<double> sigmat1) {
	
	size_t n = ft1.size(); // Number of discretization points on the curves
	size_t k = E.size2(); // Number of risk factors

	vector<double> dXi(k);
	vector<double> sigmaPrevSq(k);
	vector<double> sigmaSq(k);
	vector<double> sigma(k);

	dXi = prod(trans(E), trans(ft1 - ft2));

	for (size_t i = 0; i < k; i++) {
		sigmaSq(k) = omega(i) + alpha(i) * pow(dXi(i), 2) + beta(i) * pow(sigmat1(i), 2);
		sigma(i) = sqrt(sigmaSq(i));
	}



	return sigma;

}

double statisticsOperations::invCDFNorm(double u, double mu, double sigma) {
	
	double q = 0.0;
	boost::math::normal norm(mu, sigma);
	q = quantile(norm, u);

	return q;

}

double statisticsOperations::invCDFT(double u, double mu, double sigma, double df) {
	double q = 0.0;
	boost::math::students_t t(df); 

	mexPrintf("jappsi");

	q = quantile(t, u);

	return q;


	
}