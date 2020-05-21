#include "pch.h"
#include "rvSim.h"

#include "../MathLibrary/statisticsOperations.h"

#include "mex.h"
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/math/distributions/normal.hpp>

#include <random>
#include <boost/random.hpp>
#include <boost/random/random_device.hpp>
#include <boost/random/normal_distribution.hpp>

using namespace boost::numeric::ublas;

// boost::random_device dev;
// boost::mt19937 gener(dev);
static boost::mt19937 gener(2);
static boost::normal_distribution<> normal(0, 1);
static boost::variate_generator<boost::mt19937&, boost::normal_distribution<> > rng(gener, normal);

// Generate n normal variables with mean m and standard deviation s
matrix<double> rvSim::gen_normal(size_t k, size_t N) {

	matrix<double> rand(k, N);

	for (size_t i = 0; i < k; i++) {
		for (size_t j = 0; j < N; j++) {
			rand(i, j) = rng();
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
				eps(i, j) = mu(i) + statisticsOperations::invCDFNorm(V(i, j), 0.0, 1.0) * sigma(i); 
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
			eps(i) = mu(i) + statisticsOperations::invCDFNorm(V(i), 0.0, 1.0) * sigma(i);
		}
	}
	else if (type == "t") {
		for (size_t i = 0; i < k; i++) {
			eps(i) = mu(i) + statisticsOperations::invCDFT(V(i), dfM(i)) * sigma(i);
		}
	}

	return eps;
}