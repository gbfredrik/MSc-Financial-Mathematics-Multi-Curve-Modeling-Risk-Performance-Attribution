#include "algo4.h"
#include "algo5.h"

#include "../MathFunctions/rvSim.h"

#include <iostream>
#include <boost/numeric/ublas/matrix.hpp>

using namespace boost::numeric::ublas;

int main() {
	
	int N = 2000;
	int m = 11000;
	int n = 6;

	// Generate a test matrix
	matrix<double> test(m, n);
	test = rvSim::gen_test(m, n);

	algo4::GC_sim(test, N);
	algo5::TC_sim(test, N);


}
