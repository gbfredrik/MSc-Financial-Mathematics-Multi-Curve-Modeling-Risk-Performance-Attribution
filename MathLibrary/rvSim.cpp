#include "pch.h"
#include "mex.h"
#include "rvSim.h"
#include "../MathLibrary/statisticsOperations.h"
#include <random>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/math/distributions/normal.hpp>

using namespace boost::numeric::ublas;

boost::numeric::ublas::matrix<double> hehe(3, 2000);
//static std::seed_seq seed{ 1, 2, 3, 4, 5 };
//static std::default_random_engine e2;
static std::normal_distribution<double> distNorm(0.0, 1.0);
static std::uniform_real_distribution<double> distU(0.0, 1.0);
static std::mt19937 e2(0);

/*
matrix<double> rvSim::gen_test(int rows, int cols) {
	std::random_device rd;
	std::mt19937 mt(rd());
	std::uniform_real_distribution<double> dist(1.0, 10.0);
	matrix<double> test(rows, cols);
	for (size_t i = 0; i < test.size1(); i++) {
		for (size_t j = 0; j < test.size2(); j++) {
			test(i, j) = dist(mt);
		}
	}
	return test;
}
*/
// Generate n normal variables with mean m and standard deviation s
matrix<double> rvSim::gen_normal(double m, double s, size_t k, size_t N) {

	//static std::random_device rd{};
	//std::seed_seq seed{ rd(), rd(), rd(), rd(), rd(), rd(), rd(), rd() };
	
	//std::seed_seq seed{ 1, 2, 3, 4, 5 };
	//static std::mt19937 e2(seed);

	//static std::default_random_engine dre;
	//std::normal_distribution<double> dist(0.0, 1.0);
	matrix<double> rand(k, N);

	for (size_t i = 0; i < k; i++) {
		for (size_t j = 0; j < N; j++) {
			rand(i, j) = distNorm(e2);
			hehe(i, j) = rand(i, j);
			//mexPrintf("%g",rand(i, j));
			//mexPrintf(" ");
		}
	}
	



	return rand;
}

// R.v.s from the gamma distribution using the same method as Matlab (Marsaglia, G. and Tsang, W.W. (2000))
double rvSim::gen_gamma(double df) {

	return std::tgamma(df);

}


matrix<double> rvSim::genEps(matrix<double> V, vector<double> mu, vector<double> sigma, std::string type, 
	vector<double> dfM) {
	
	size_t k = V.size1();
	size_t N = V.size2();

	matrix<double> eps(k, N);
	
	if (type == "normal") {
		for (size_t i = 0; i < k; i++) {
			for (size_t j = 0; j < N; j++) {
				eps(i, j) = mu(i) + statisticsOperations::invCDFNorm(V(i, j), 0, 1) * sigma(i); 
			}
		}
	}
	else if (type == "t") {
		for (size_t i = 0; i < k; i++) {
			for (size_t j = 0; j < N; j++) {
				eps(i, j) = mu(i) + statisticsOperations::invCDFT(V(i, j), dfM(i)) * sigma(i);
			}
		}
	}

	
	return eps;
}

vector<double> rvSim::genEps(vector<double> V, vector<double> mu, vector<double> sigma, std::string type,
	vector<double> dfM) {

	size_t k = V.size();
	vector<double> eps(k);

	if (type == "normal") {
		for (size_t i = 0; i < k; i++) {
			eps(i) = mu(i) + statisticsOperations::invCDFNorm(V(i), 0, 1) * sigma(i);
		}
	}
	else if (type == "t") {
		for (size_t i = 0; i < k; i++) {
			eps(i) = mu(i) + statisticsOperations::invCDFT(V(i), dfM(i)) * sigma(i);
		}
	}


	return eps;
}