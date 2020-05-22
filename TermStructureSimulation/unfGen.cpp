#include "pch.h"
#include "unfGen.h"
#include "mex.h"

#include <numeric>

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/math/distributions/normal.hpp>
#include <boost/math/distributions/students_t.hpp>


#include "../MathLibrary/matrixOperations.h"
#include "../MathLibrary/rvSim.h"
#include "../MathLibrary/statisticsOperations.h"


using namespace boost::numeric::ublas;
using namespace boost::math;


matrix<double> unfGen::genU(matrix<double> const& rho, int N, std::string copula, vector<double> const& df) {
	
	size_t k = rho.size1();
	matrix<double> U(k, N);

	if (copula == "normal") {
		U = GC_sim(rho, N);
	}
	else if (copula == "t") {
		U = TC_sim(rho, N, df);
	}

	return U;
}


matrix<double> unfGen::TC_sim(matrix<double> const& rho, int N, vector<double> const& dfC) {

	size_t k = rho.size1();

	// Perform Cholesky decomposition
	matrix<double> L(k, k);
	L = matrixOperations::chol(rho);

	// Generate i.i.d. standard normal random variables
	matrix<double> X(k, N);
	X = rvSim::gen_normal(0.0, 1.0, k, N);
	// Generate Correlated Gaussian samples
	matrix<double> Z(k, N);
	Z = prod(L, X);

	// Generate random variables from the gamma distribution
	// Generate sqrt(normalized chi-square r.v.s)
	vector<double> g(k);
	vector<double> Xi(k);
	for (size_t i = 0; i < k; i++) {
		g(i) = rvSim::gen_gamma(dfC(i));
		Xi(i) = sqrt(g(i) / dfC(i));
		for (int j = 0; j < N; j++) {
			Z(i, j) = Z(i, j) / Xi(i);
		}
	}
	// Return uniformly distributed r.v.s	
	matrix<double> U(k, N);
	
	for (size_t i = 0; i < k; i++) {
		students_t dist(dfC(i));
		for (int j = 0; j < N; j++) {
			U(i, j) = cdf(dist, Z(i, j));
		}
	}
	return U;
}



matrix<double> unfGen::GC_sim(matrix<double> const& rho, int N) {

	size_t k = rho.size1();

	// Perform Cholesky decomposition
	matrix<double> L(k, k);
	L = matrixOperations::chol(rho);

	// Generate i.i.d. standard normal random variables
	matrix<double> X(k, N);
	X = rvSim::gen_normal(0.0, 1.0, k, N);
	
	// Generate Correlated Gaussian samples
	matrix<double> Z(k, N);
	Z = prod(L, X);

	// Return correlated uniformy distributed random variables
	matrix<double> U(k, N);					
	normal s(0.0, 1.0);
	for (size_t i = 0; i < k; i++) {
		for (int j = 0; j < N; j++) {
			U(i, j) = cdf(s, Z(i, j)); 
		}
	}
	return U;
}




