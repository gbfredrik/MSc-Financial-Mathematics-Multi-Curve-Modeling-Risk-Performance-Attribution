#include "pch.h"
#include "rvSim.h"

#include "../MathLibrary/statisticsOperations.h"

#include "mex.h"
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/math/distributions/normal.hpp>

#include <random>

using namespace boost::numeric::ublas;

matrix<double> rvSim::gen_test(int rows, int cols) {
	std::random_device rd;
	std::mt19937 mt(rd());
	std::uniform_real_distribution<double> dist(1.0, 10.0);
	matrix<double> test(rows, cols);
	for (size_t i = 0; i < test.size1(); ++i) {
		for (size_t j = 0; j < test.size2(); ++j) {
			test(i, j) = dist(mt);
		}
	}

	return test;
}

// Generate n normal variables with mean m and standard deviation s
matrix<double> rvSim::gen_normal(double m, double s, int k, int N) {
	static std::random_device rd;
	static std::mt19937 e2(rd());
	std::normal_distribution<> dist(m, s);
	matrix<double> rand(k, N);

	for (int i = 0; i < k; ++i) {
		for (int j = 0; j < N; ++j) {
			rand(i, j) = dist(e2);
		}
	}
	
	return rand;
}

// R.v.s from the gamma distribution using the same method as Matlab (Marsaglia, G. and Tsang, W.W. (2000))
double rvSim::gen_gamma(double a) {
	double d, c, x, v, u;
	d = a - 1.0 / 3.0;
	c = 1.0 / sqrt(9.0 * d);

	while (true) {
		do {
			x = gen_normal(0.0, 1.0);
			v = pow(1.0 + c * x, 3);
		} while (v <= 0);

		u = gen_uniform(0.0, 1.0);
		if (u < 1 - 0.0331 * pow(x, 4)) || (log(u) < 0.5 * pow(x, 2) + d * (1.0 - v + log(v))) {
			return d * v;
		}
	}
}

double rvSim::gen_normal(double m, double s) {
	std::random_device rd;
	std::default_random_engine generator(rd());
	std::normal_distribution<double> distribution(m, s);

	return distribution(generator);
}

double rvSim::gen_uniform(double l, double u) {
	std::random_device rd;
	std::default_random_engine generator(rd());
	std::uniform_real_distribution<double> distribution(l, u);

	return distribution(generator);
}

matrix<double> rvSim::genEps(
	matrix<double> const& V, 
	matrix<double> const& E, 
	vector<double> const& sigma, 
	std::string const& type
) {
	size_t m = V.size1();
	size_t n = V.size2();

	matrix<double> cov(m, m);
	cov = statisticsOperations::covm(E);

	matrix<double> eps(m, n);

	if (type == "normal") {
		for (size_t i = 0; i < m; ++i) {
			for (size_t j = 0; j < n; ++j) {
				eps(i, j) = statisticsOperations::invCDFNorm(V(i, j), 0, sqrt(cov(i, i))) * sigma(i);
			}
		}
	} else if (type == "t") {
		for (size_t i = 0; i < m; ++i) {
			for (size_t j = 0; j < n; ++j) {
				eps(i, j) = statisticsOperations::invCDFT(V(i, j), 0, sqrt(cov(i, i)), 5) * sigma(i);
			}
		}
	}

	return eps;
}
