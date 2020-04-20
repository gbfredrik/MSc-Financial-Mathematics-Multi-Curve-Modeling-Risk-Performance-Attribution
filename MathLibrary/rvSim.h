#ifndef RVSIM_H
#define RVSIM_H

#include <boost/numeric/ublas/matrix.hpp>


class rvSim {
private:

public:
	static boost::numeric::ublas::matrix<double> gen_test(int rows, int cols);
	static boost::numeric::ublas::vector<double> gen_normal(double m, double s, int N);
	static float gen_gamma(float a);
	static float gen_normal(double m, double s);
	static float gen_uniform(float l, float u);
};


#endif
