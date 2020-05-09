#include "pch.h"
#include "mex.h"
#include "unfGenGauss.h"
#include "../MathLibrary/matrixOperations.h"
#include "../MathLibrary/rvSim.h"
#include "../MathLibrary/statisticsOperations.h"

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/math/distributions/normal.hpp>


using namespace boost::numeric::ublas;
using namespace boost::math;


matrix<double> unfGenGauss::GC_sim(matrix<double> const& E, int N) {

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
	Z = prod(L, X);
	
	// Return correlated uniformy distributed random variables
	matrix<double> U(k, N);
	normal s;
	for (int i = 0; i < k; i++) {
		for (int j = 0; j < N; j++) {
			U(i, j) = cdf(s, Z(i, j)); //CHANGE
		}
	}
	return U;
}

