#include "pch.h"
#include "mex.h"

#include "MultipleYieldSim.h"
#include "unfGen.h"
#include "varRed.h"
#include "../MathLibrary/rvSim.h"
#include "../MathLibrary/statisticsOperations.h"


#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/optional/optional.hpp>
#include <boost/none.hpp>
#include <boost/none_t.hpp>
using namespace boost::numeric::ublas;


void MultipleYieldSim::simMultipleFull(vector<matrix<double>> const& E, vector<matrix<double>> const& rho, vector<vector<double>> const& mu, 
	vector<vector<double>> const& omega, vector<vector<double>> const& alpha, vector<vector<double>> const& beta,
	vector<matrix<double>> const& hist, vector<std::string> marginal, vector<std::string> copula, 
	vector<std::string> varRedType, size_t d, size_t N, vector<matrix<double>> & fRes, 
	boost::optional<vector<vector<double>>> const& gamma, boost::optional<vector<double>> const& kappa,
	boost::optional<vector<vector<double>>> const& xiHat, boost::optional<vector<vector<double>>> const& dfC,
	boost::optional<vector<vector<double>>> const& dfM){


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

		muC: matrix, mean value parameter for gaussian copulas
		muM: matrix, mean value parameter for gaussian marginals

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
	vector<double> fZero(n); // Starting vector in the simulation
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
	vector<matrix<double>> histPrevSim(M, matrix<double>(n, N)); 
	vector<matrix<double>> histPrevPrevSim(M, matrix<double>(n, N));
	fZero = row(hist(0), m - 1);
	for (size_t k = 0; k < N; k++) {
		column(histPrevPrevSim(0), k) = fZero;
	}
	for (size_t k = 1; k < M; k++) {
		column(pi, k - 1) = row(hist(k), m - 1);
		for (size_t l = 0; l < N; l++) {
			column(histPrevPrevSim(k), l) = column(pi, k - 1);
		}
	}

	for (size_t i = 0; i < d; i++) {
		if (i == 0) { // Simulate N curves the first day
			for (size_t j = 0; j < M; j++) {
				U(j) = unfGen::genU(rho(j), N, copula(j), boost::get(dfC)(j)); // Generate uniformly correlated random variables with desired copula
				V(j) = varRed::redVariance(U(j), varRedType(j)); // Apply desired variance reduction			
			
				sigma(j) = statisticsOperations::GARCH(omega(j), alpha(j),
					beta(j), boost::get(gamma)(j), E(j), hist(j)); // Calculate scaling factor
				eps(j) = rvSim::genEps(U(j), mu(j), sigma(j), marginal(j), boost::get(dfM)(j));
			}
			simMultipleDaily(E, fZero, pi, eps, M, N, fRes, histPrevSim);
		}
		if (i != 0) { // Simulate 1 curve for each following day
			for (size_t j = 0; j < M; j++) {
				U(j) = unfGen::genU(rho(j), N, copula(j), boost::get(dfC)(j)); // Generate uniformly correlated random variables with desired copula
				V(j) = varRed::redVariance(U(j), varRedType(j)); // Apply desired variance reduction			

				for (size_t l = 0; l < N; l++) {
					sigma(j) = statisticsOperations::GARCH(omega(j), alpha(j),
						beta(j), boost::get(gamma)(j), E(j), column(histPrevSim(j), l), column(histPrevPrevSim(j), l), sigma(j)); // Calculate scaling factor
					column(eps(j), l) = rvSim::genEps(column(U(j), l), mu(j), sigma(j), marginal(j), boost::get(dfM)(j));
				}
			}
			simSingleMultipleDaily(E, eps, M, N, fRes, histPrevSim, histPrevPrevSim);
		}
	}
}


/*
	Function to simulate one day ahead N times.
*/
void MultipleYieldSim::simMultipleDaily(vector<matrix<double>> const& E, vector<double> const& fZero, 
	matrix<double> const& pi, vector<matrix<double>> const& eps, int M, int N, vector<matrix<double>>& fRes,
	vector<matrix<double>>& histPrevSim) {

	for (int i = 0; i < N; i++) { // Simulate N times
		for (int k = 0; k < M; k++) { // Iterate over each tenor curve
			if(k == 0){
				column(fRes(k), i) = fZero + prod(E(k), column(eps(k), i)); // +prod(E(k)(), (kappa(k) * (xiHat(k) - ) + eps(k)(i, j)); // Simulate risk-free curve
				column(histPrevSim(k), i) = column(fRes(k), i);
			}
			else {
				column(fRes(k), i) = column(fRes(k - 1), i) + column(pi, k -  1) + prod(E(k), column(eps(k), i)); //+ kappa(k) * (xiHat(k) - ) + eps(k)(i, j); // Simulate tenor curves
				column(histPrevSim(k), i) = column(fRes(k), i) - column(fRes(0), i);
			}
		}
	}
}

/*
	Function to simulate d days ahead 1 time each.
*/
void MultipleYieldSim::simSingleMultipleDaily(vector<matrix<double>> const& E, vector<matrix<double>> const& eps, int M, int N,
	vector<matrix<double>>& fRes, vector<matrix<double>>& histPrevSim, vector<matrix<double>>& histPrevPrevSim) {
	
	size_t n = E(0).size1();
	matrix<double> pi(n, M - 1);


	
	for (int i = 0; i < N; i++) { // Simulate N times
		for (int k = 0; k < M; k++) { // Iterate over each tenor curve
			if (k != M - 1) {
				column(pi, k) = column(fRes(k + 1), i) - column(fRes(0), i);
			}
			
			if (k == 0) {
				column(histPrevPrevSim(k), i) = column(fRes(k), i);
				column(fRes(k), i) = column(fRes(k), i) + prod(E(k), column(eps(k), i)); // +prod(E(k)(), (kappa(k) * (xiHat(k) - ) + eps(k)(i, j)); // Simulate risk-free curve
				column(histPrevSim(k), i) = column(fRes(k), i);
			}
			else {
				column(histPrevPrevSim(k), i) = column(fRes(k), i) - column(fRes(0), k);
				column(fRes(k), i) = column(fRes(k - 1), i) + column(pi, k - 1) + prod(E(k), column(eps(k), i)); //+ kappa(k) * (xiHat(k) - ) + eps(k)(i, j); // Simulate tenor curves
				column(histPrevSim(k), i) = column(fRes(k), i) - column(fRes(0), k);
			}
		}
	}
}




