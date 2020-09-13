#include "pch.h"
#include "mex.h"

#include "FXMultipleYieldSim.h"
#include "unfGen.h"
#include "varRed.h"
#include "../MathLibrary/rvSim.h"
#include "../MathLibrary/statisticsOperations.h"

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/optional/optional.hpp>

using namespace boost::numeric::ublas;

//boost::numeric::ublas::matrix<double> test(6, 2000);

void FXMultipleYieldSim::FXsimMultipleFull(vector<matrix<double>> const& E, vector<matrix<double>> const& rho, vector<vector<double>> const& mu,
	vector<vector<double>> const& omega, vector<vector<double>> const& alpha, vector<vector<double>> const& beta,
	vector<matrix<double>> const& hist, vector<std::string> marginal, vector<std::string> copula,
	vector<std::string> varRedType, size_t d, size_t N, vector<matrix<double>>& fRes,
	boost::optional<vector<vector<double>>> const& gamma, boost::optional<vector<double>> const& kappa,
	boost::optional<vector<vector<double>>> const& xiHat, boost::optional<vector<vector<double>>> const& dfC,
	boost::optional<vector<vector<double>>> const& dfM) {
	/*
		E: vector, eigenvector matrices
		rho: vector, matrices describing risk factor dependencies

		omega: vector, garch parameter
		alpha: vector, garch parameter
		beta: vector, garch parameter
		hist: vector, curve data matrices

		marginal: vector, marginal distributions used
		copula: vector, copulas used
		varRedType: vector, variance reduction type used

		d: size_t, days ahead to simulate
		N: size_t, number of simulations

		fRes: result vector, one matrix for each tenor

		OPTIONAL PARAMETERS:
		gamma: vector, garch parameter

		kappa: vector, mean reversion speed parameter of each tau
		xiHat: matrix, long term averages of each tau and risk factor


		dfC: matrix, degrees of freedom for student's t-copulas
		dfM: matrix, degrees of freedom for student's t-marginals

	*/

	size_t M = E.size(); // Risk-free curve + number of tenor curves
	vector<size_t> k(M); // Number of risk factors of each tenor
	for (size_t i = 0; i < M; i++) {
		k(i) = E(i).size2();
	}
	size_t n = fRes(0).size1(); // Number of discretization points
	size_t m = hist(0).size1(); // Number of historical data points
	vector<matrix<double>> U(M); // Uniform r.v.s
	vector<matrix<double>> V(M); // Variance reduced r.v.s
	vector<matrix<double>> eps(M); // R.v.s used in the simulation
	vector<vector<double>> fStart(n); // Starting vector in the simulation
	matrix<double> pi(n, M - 1); // Starting vector in the simulation
	vector<vector<double>> sigma(M);

	/*
		Set number the number of columns in the individual matrices
		to the number of risk factors of that tenor
	*/
	for (size_t i = 0; i < M; i++) {
		U(i) = matrix<double>(k(i), N);
		V(i) = matrix<double>(k(i), N);
		eps(i) = matrix<double>(k(i), N);
		sigma(i) = vector<double>(k(i));
	}


	/*
		histPrevSim and histPrevPrevSim keeps track of the
		last and second to last simulated curves and is used to calculate the garch volatility.
		Contains fZero and basis spreads (pi).
	*/

	fStart(0) = row(hist(0), m - 1); //Domestic
	fStart(1) = row(hist(1), m - 1); //Foreign
	fStart(2) = row(hist(2), m - 1); //Demand


	// Simulate N curves the first day
	for (size_t j = 0; j < M; j++) {
		U(j) = unfGen::genU(rho(j), N, copula(j), boost::get(dfC)(j)); // Generate uniformly correlated random variables with desired copula
		V(j) = varRed::redVariance(U(j), varRedType(j)); // Apply desired variance reduction technique		
		//test = U(0);
		sigma(j) = statisticsOperations::GARCH(omega(j), alpha(j),
			beta(j), boost::get(gamma)(j), E(j), hist(j)); // Calculate scaling factor
		eps(j) = rvSim::genEps(V(j), mu(j), sigma(j), marginal(j), boost::get(dfM)(j));
	}
	FXsimMultipleDaily(E, fStart, pi, eps, M, N, fRes);
}


/*
	Function to simulate one day ahead N times.
*/
void FXMultipleYieldSim::FXsimMultipleDaily(vector<matrix<double>> const& E, vector<vector<double>> const& fStart,
	matrix<double> const& pi, vector<matrix<double>> const& eps, int M, int N, vector<matrix<double>>& fRes) {

	for (int i = 0; i < N; i++) { // Simulate N times
		for (int k = 0; k < M; k++) { // Iterate over each curve
				column(fRes(k), i) = fStart(k) + prod(E(k), column(eps(k), i)); // +prod(E(k)(), (kappa(k) * (xiHat(k) - ) + eps(k)(i, j)); // Simulate risk-free curve
		}
	}
}

