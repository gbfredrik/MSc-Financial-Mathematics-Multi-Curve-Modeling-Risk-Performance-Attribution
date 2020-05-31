#include "pch.h"
#include "RvSim.h"

#include "../MathLibrary/statisticsOperations.h"

//#include "mex.h"
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/math/distributions/normal.hpp>
#include <boost/random.hpp>
#include <boost/random/random_device.hpp>
#include <boost/random/normal_distribution.hpp>

#include <random>

using namespace boost::numeric;

// boost::random_device dev;
// boost::mt19937 gener(dev);
static boost::mt19937 gener(2);
static boost::normal_distribution<> normal(0, 1);
static boost::variate_generator<boost::mt19937&, boost::normal_distribution<>> rng(gener, normal);

// Generate k times n normal variables
ublas::matrix<double> rvSim::gen_normal(size_t const k, size_t const N) {
	ublas::matrix<double> rand(k, N);

	for (size_t i{ 0 }; i < k; ++i) {
		for (size_t j{ 0 }; j < N; ++j) {
			rand(i, j) = rng();
		}
	}
	
	return rand;
}

double rvSim::gen_gamma(double const df) {
	return std::tgamma(df);
}

ublas::matrix<double> rvSim::genEps(
	ublas::matrix<double> const& V,
	ublas::vector<double> const& mu,
	ublas::vector<double> const& sigma,
	std::string const& type,
	ublas::vector<double> const& dfM
) {
	size_t k{ V.size1() };
	size_t N{ V.size2() };

	ublas::matrix<double> eps(k, N);
	
	// Todo: Kolla om nedan kan ers�ttas av:
	// for (size_t i{ 0 }; i < k; ++i) {
	//	   row(eps, i) = genEps(row(V, i), mu, sigma, type, dfM);
	// }

	if (type == "normal") {
		for (size_t i{ 0 }; i < k; ++i) {
			for (size_t j{ 0 }; j < N; ++j) {
				eps(i, j) = mu(i) + statisticsOperations::invCDFNorm(V(i, j), 0, 1) * sigma(i); 
			}
		}
	} else if (type == "t") {
		for (size_t i{ 0 }; i < k; ++i) {
			for (size_t j{ 0 }; j < N; ++j) {
				eps(i, j) = mu(i) + statisticsOperations::invCDFT(V(i, j), dfM(i)) * sigma(i);
			}
		}
	}

	return eps;
}

ublas::vector<double> rvSim::genEps(
	ublas::vector<double> const& V,
	ublas::vector<double> const& mu,
	ublas::vector<double> const& sigma,
	std::string const& type,
	ublas::vector<double> const& dfM
) {
	size_t k{ V.size() };
	ublas::vector<double> eps(k);

	if (type == "normal") {
		for (size_t i{ 0 }; i < k; ++i) {
			eps(i) = mu(i) + statisticsOperations::invCDFNorm(V(i), 0, 1) * sigma(i);
		}
	} else if (type == "t") {
		for (size_t i{ 0 }; i < k; ++i) {
			eps(i) = mu(i) + statisticsOperations::invCDFT(V(i), dfM(i)) * sigma(i);
		}
	}

	return eps;
}
