#include "pch.h"

#include "rvSim.h"

#include <random>
#include <boost/numeric/ublas/matrix.hpp>

using namespace boost::numeric::ublas;


matrix<double> rvSim::gen_test(int rows, int cols) {
	std::random_device rd;
	std::mt19937 mt(rd());
	std::uniform_real_distribution<double> dist(1.0, 10.0);
	matrix<double> test(rows, cols);
	for (int i = 0; i < test.size1(); i++) {
		for (int j = 0; j < test.size2(); j++) {
			test(i, j) = dist(mt);
		}
	}
	return test;
}

// Generate n normal variables with mean m and standard deviation s
vector<double> rvSim::gen_normal(double m, double s, int N) {
	std::random_device rd;
	std::default_random_engine generator(rd());
	std::normal_distribution<double> distribution(m, s);
	vector<double> rand(N);

	for (int i = 0; i < N; ++i) {
		double number = distribution(generator);
		rand[i] = number;
	}
	return rand;
}

// R.v.s from the gamma distribution using the same method as Matlab (Marsaglia, G. and Tsang, W.W. (2000))
float rvSim::gen_gamma(float a) {
	float d, c, x, v, u;
	d = a - 1. / 3.;
	c = 1. / sqrt(9. * d);

	for (;;) {
		do {
			x = gen_normal(0.0, 1.0);
			v = pow(1. + c * x, 3);
		} while (v <= 0);

		u = gen_uniform(0.0, 1.0);
		if (u < 1 - 0.0331 * pow(x, 4)) {
			return d * v;
		}

		if (log(u) < 0.5 * pow(x, 2) + d * (1. - v + log(v))) {
			return (d * v);
		}
	}
}

float rvSim::gen_normal(double m, double s) {
	std::random_device rd;
	std::default_random_engine generator(rd());
	std::normal_distribution<double> distribution(m, s);
	return distribution(generator);
}

float rvSim::gen_uniform(float l, float u) {
	std::random_device rd;
	std::default_random_engine generator(rd());
	std::uniform_real_distribution<double> distribution(l, u);
	return distribution(generator);
}

