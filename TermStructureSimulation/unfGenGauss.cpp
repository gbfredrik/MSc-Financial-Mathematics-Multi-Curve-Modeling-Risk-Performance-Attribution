#include "pch.h"
#include "unfGenGauss.h"
#include "../MathLibrary/matrixOperations.h"
#include "../MathLibrary/rvSim.h"
#include "../MathLibrary/statistics.h"

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/math/tools/bivariate_statistics.hpp>
#include <boost/math/distributions/normal.hpp>


using namespace boost::numeric::ublas;
using namespace boost::math;


matrix<double> unfGenGauss::GC_sim(matrix<double> const& E, int N) {

	size_t m = E.size1();
	size_t n = E.size2();

	// Calculate the correlation matrix
	matrix<double> corr(n, n);
	corr = statistics::corrm(E);

	// Perform Cholesky decomposition
	matrix<double> L(n, n);
	L = matrixOperations::chol(corr);

	// Generate i.i.d. standard normal random variables
	matrix<double> X(N, n);
	for (int i = 0; i < n; i++) {
		column(X, i) = rvSim::gen_normal(0.0, 1.0, N);
	}

	// Generate Correlated Gaussian samples
	matrix<double> Z(N, n);
	Z = prod(X, trans(L));

	// Return correlated uniformy distributed random variables
	matrix<double> U(N, n);
	normal s;
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < n; j++) {
			U(i, j) = cdf(s, Z(i, j));
		}
	}
	return U;
}

