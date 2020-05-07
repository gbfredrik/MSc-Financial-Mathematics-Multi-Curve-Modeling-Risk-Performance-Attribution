#include "pch.h"
#include "MultipleYieldSim.h"
#include "unfGenT.h"
#include "unfGenGauss.h"
#include "lhsd.h"
#include "../MathLibrary/rvSim.h"
#include "mex.h"

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>



using namespace boost::numeric::ublas;


void MultipleYieldSim::simMultipleFull(vector<matrix<double>> const& E, vector<double> const& fZero, matrix<double> const& pi, vector<double> const& kappa, matrix<double> const& xiHat, int d, int N, vector<matrix<double>>& fRes) {
	
	/*
		E: vector, eigenvector matrices
		fZero: vector, forward rate curve of time t
		pi: matrix, basis spreads
		kappa: vector, mean reversion speed parameter of each tau
		xiHat: vector, long term averages of each tau
		d: int, days ahead to simulate
		N: int, number of simulations
		fRes: vector, one matrix for each tenor	
	*/
	
	size_t M = E.size(); // Risk-free curve + number of tenor curves
	vector<int> k(M);
	for (int i = 0; i < M; i++) {
		k(i) = E(i).size2(); // Get number of risk factors for each tau
	}
	
	/*
	vector<matrix<double>> U; 
	vector<matrix<double>> V;
	vector<matrix<double>> eps;
	
		Set number the number of columns in the individual matrices
		to the number of risk factors of that tenor
	*/
	/*
	for (int i = 0; i < M; i++) { 
		U(i)(N, k(i));
		V(i)(N, k(i));
		eps(i)(N, k(i));
	}
	*/
	vector<matrix<double>> U(M, matrix<double>(6, N));
	vector<matrix<double>> V(M, matrix<double>(6, N));
	vector<matrix<double>> eps(M, matrix<double>(6, N));



	//vector<double> sigma;
	double sigma = 0.001;
	

	for (int i = 0; i < d; i++) {
		for (int j = 0; j < M; j++) {
			U(j) = unfGenT::TC_sim(E(j), N); // Generate uniformly correlated random variables
			V(j) = lhsd::lhsd_gen(U(j)); // Apply LHSD
			//sigma(j) = garch::garch_gen(omega, alpha, gamma, beta); // Calculate scaling factor
			eps(j) = rvSim::gen_eps(U(j), sigma, "normal");
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
//void simMultipleDaily(vector<matrix<double>> const& E, vector<double> const& fZero, matrix<double> const& pi, vector<matrix<double>> const& eps, vector<double> const& kappa, matrix<double> xiHat, int M, int N, vector<matrix<double>> *fRes) {
void MultipleYieldSim::simMultipleDaily(vector<matrix<double>> const& E, vector<double> const& fZero, matrix<double> const& pi, vector<matrix<double>> eps, int M, int N, vector<matrix<double>>& fRes) {

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



