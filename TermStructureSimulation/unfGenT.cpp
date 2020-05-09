#include "pch.h"
#include "mex.h"
#include "unfGenT.h"
#include "../MathLibrary/matrixOperations.h"
#include "../MathLibrary/rvSim.h"
#include "../MathLibrary/statisticsOperations.h"

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/math/distributions/normal.hpp>
#include <boost/math/distributions/students_t.hpp>


using namespace boost::numeric::ublas;
using namespace boost::math;


matrix<double> unfGenT::TC_sim(matrix<double> const& E, int N) {

	size_t m = E.size1();
	size_t k = E.size2();
	// Calculate the correlation matrix
	matrix<double> cov(k, k);
	cov = statisticsOperations::covm(E);

	// Perform Cholesky decomposition
	matrix<double> L(k, k);
	L = matrixOperations::chol(cov);

	// Generate i.i.d. standard normal random variables
	matrix<double> X(k, N);
	X = rvSim::gen_normal(0.0, 1.0, k, N);
	// Generate Correlated Gaussian samples
	matrix<double> Z(k, N);
	Z = prod(trans(L), X);

	// Generate random variables from the gamma distribution
	// Generate sqrt(normalized chi-square r.v.s)
	int df = 3; // Degrees of freedom
	vector<double> g(k);
	vector<double> Xi(k);
	for (int i = 0; i < k; i++) {
		g(i) = rvSim::gen_gamma(df);
		Xi(i) = sqrt(g(i) / df);
		for (int j = 0; j < N; j++) {
			Z(i, j) = Z(i, j) / Xi(i);
		}
	}
	// Return uniformly distributed r.v.s
	matrix<double> U(k, N);
	students_t dist(df);
	for (int i = 0; i < k; i++) {
		for (int j = 0; j < N; j++) {
			U(i, j) = cdf(dist, Z(i, j));
		}
	}
	return U;

}








