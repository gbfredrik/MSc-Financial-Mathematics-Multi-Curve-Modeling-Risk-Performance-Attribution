#include "pch.h"
#include "MultipleYieldSim.h"
#include "unfGenT.h"
#include "lhsd.h"
#include "mex.h"

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>

using namespace boost::numeric::ublas;


void MultipleYieldSim::simMultipleFull(vector<matrix<double>> const& E, vector<double> const& fZero, matrix<double> const& pi, vector<double> const& kappa, matrix<double> const& xiHat, int d, vector<matrix<double>>& fRes) {
	
	/*
		E: vector with an arbitrary number of eigenvector matrices
		f: forward rate curves of time = t
		d: days ahead to simulate
		fRes: result vector with one matrix for each tenor curve	
	*/
	
	int M = E.size(); // Risk-free curve + number of tenor curves
	int N = 1; // Number of simulations each day 
	//vector<matrix<double>> U;
	//vector<matrix<double>> V;
	//vector<matrix<double>> eps;
	//vector<double> sigma;

	for (int i = 0; i < d; i++) {
		//for (int j = 0; j < M; j++) {
			//U(j) = unfGenT::TC_sim(E(j), N); // Generate uniformly correlated random variables
			//V(j) = lhsd::lhsd_gen(U(j)); // Apply LHSD
			//sigma(j) = garch::garch_gen(omega, alpha, gamma, beta); // Calculate scaling factor
			//double sigma = 0.1;
			//eps(j) = V(j) / sigma; // Scale epsilon with the garch volatility
		//}
		
		if (i == 0){ // Simulate N curves the first day
			//simMultipleDaily(E, fZero, pi, eps, kappa, xiHat, M, N, fRes);
			simMultipleDaily(E, fZero, pi, M, N, fRes);
		}
		else { // Simulate 1 curve for each following day
		}
	}
}
/*
	Function to simulate one day ahead N times.
*/
//void simMultipleDaily(vector<matrix<double>> const& E, vector<double> const& fZero, matrix<double> const& pi, vector<matrix<double>> const& eps, vector<double> const& kappa, matrix<double> xiHat, int M, int N, vector<matrix<double>> *fRes) {
void MultipleYieldSim::simMultipleDaily(vector<matrix<double>> const& E, vector<double> const& fZero, matrix<double> const& pi, int M, int N, vector<matrix<double>>& fRes) {
	for (int i = 0; i < N; i++) { // Simulate N times
		for (int k = 0; k < M; k++) { // Iterate over each tenor curve
			if(k == 0){
				column(fRes(k), i) = fZero; // +prod(E(k)(), (kappa(k) * (xiHat(k) - ) + eps(k)(i, j)); // Simulate risk-free curve
			}
			else {
				column(fRes(k), i) = column(fRes(k - 1), i) + column(pi, k -  1); //+ kappa(k) * (xiHat(k) - ) + eps(k)(i, j); // Simulate tenor curves
			}
		}
	}
	mexPrintf("yo");
}

/*
	Function to simulate M days ahead 1 time each.
*/
/*
void simSingleMultipleDaily() {

}
*/



