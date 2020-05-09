#include "pch.h"
#include "MultipleYieldSim.h"
#include "unfGenT.h"
#include "unfGenGauss.h"
#include "lhsd.h"
#include "../MathLibrary/rvSim.h"
#include "../MathLibrary/statisticsOperations.h"
#include "mex.h"

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>



using namespace boost::numeric::ublas;


void MultipleYieldSim::simMultipleFull(vector<matrix<double>> const& E, vector<double> const& fZero, 
	matrix<double> const& pi, vector<double> const& kappa, matrix<double> const& xiHat, vector<double> const& omega, 
	vector<double> const& alpha, vector<double> const& beta, vector<matrix<double>> const& hist, int d, int N, vector<matrix<double>>& fRes) {
	
	/*
		E: vector, eigenvector matrices
		fZero: vector, forward rate curve of time t
		pi: matrix, basis spreads
		kappa: vector, mean reversion speed parameter of each tau
		xiHat: vector, long term averages of each tau
		omega: vector, garch parameter								CHANGE TO MATRIX
		alpha: vector, garch parameter                              CHANGE TO MATRIX
		beta: vector, garch parameter								CHANGE TO MATRIX
		hist: vector, curve data matrices
		d: int, days ahead to simulate
		N: int, number of simulations
		fRes: vector, one matrix for each tenor	
	*/
	
	size_t M = E.size(); // Risk-free curve + number of tenor curves
	vector<int> k(M);
	for (int i = 0; i < M; i++) {
		k(i) = E(i).size2(); // Get number of risk factors for each tau
	}
	
	vector<matrix<double>> U(M); 
	vector<matrix<double>> V(M);
	vector<matrix<double>> eps(M);
	

	/*
		Set number the number of columns in the individual matrices
		to the number of risk factors of that tenor
	*/
	
	for (int i = 0; i < M; i++) { 
		U(i) = matrix<double>(k(i), N);
		V(i) = matrix<double>(k(i), N);
		eps(i) = matrix<double>(k(i), N);
	}
	
	vector<double> sigma(k(0)); // FIX FOR ARBITRARY NUMBER OF RISK FACTORS


	for (int i = 0; i < d; i++) {
		for (int j = 0; j < M; j++) { // FIX, CHOOSE COPULA AS INPUT
			U(j) = unfGenGauss::GC_sim(E(j), N); // Generate uniformly correlated random variables
			V(j) = lhsd::lhsd_gen(U(j)); // Apply LHSD
			
			sigma = statisticsOperations::GARCH(omega, alpha, beta, E(j), hist(j)); // Calculate scaling factor
			eps(j) = rvSim::gen_eps(V(j), sigma, "normal");
		}
		
		if (i == 0){ // Simulate N curves the first day
			//simMultipleDaily(E, fZero, pi, eps, kappa, xiHat, M, N, fRes);
			simMultipleDaily(E, fZero, pi, eps, M, N, fRes);
		}
		else { // Simulate 1 curve for each following day
		}
	}
}

/*
	Function to simulate one day ahead N times.
*/
void MultipleYieldSim::simMultipleDaily(vector<matrix<double>> const& E, vector<double> const& fZero, 
	matrix<double> const& pi, vector<matrix<double>> const& eps, int M, int N, vector<matrix<double>>& fRes) {

	for (int i = 0; i < N; i++) { // Simulate N times
		for (int k = 0; k < M; k++) { // Iterate over each tenor curve
			if(k == 0){
				column(fRes(k), i) = fZero + prod(E(k), column(eps(k), i)); // +prod(E(k)(), (kappa(k) * (xiHat(k) - ) + eps(k)(i, j)); // Simulate risk-free curve
			}
			else {
				column(fRes(k), i) = column(fRes(k - 1), i) + column(pi, k -  1) + prod(E(k), column(eps(k), i)); //+ kappa(k) * (xiHat(k) - ) + eps(k)(i, j); // Simulate tenor curves
			}
		}
	}
}

/*
	Function to simulate M days ahead 1 time each.
*/
/*
void simSingleMultipleDaily() {

}
*/



