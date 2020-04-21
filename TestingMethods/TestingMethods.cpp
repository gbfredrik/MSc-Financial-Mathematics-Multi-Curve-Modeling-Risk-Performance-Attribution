#include "../TermStructureSimulation/unfGenGauss.h"
#include "../TermStructureSimulation/unfGenT.h"
#include "../TermStructureSimulation/lhsd.h"
#include "../MathLibrary/rvSim.h"

#include <iostream>
#include <numeric>

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/range/counting_range.hpp>


using namespace boost::numeric::ublas;
using namespace boost::range;
void test_algo4_5();

int main() {

	test_algo4_5();

}



void test_algo4_5() {
	int N = 2000;
	int m = 11000;
	int n = 6;

	
	// Generate a test matrix
	matrix<double> test(m, n);
	test = rvSim::gen_test(m, n);
	
	// Generate uniformly distributed variables using the Gaussian and t-copula
	matrix<double> U_Gauss(N, n);
	U_Gauss = unfGenGauss::GC_sim(test, N);
	matrix<double> U_T(N, n);
	U_T = unfGenT::TC_sim(test, N);
	

	// Perform LHSD
	matrix<double> V(N, n);
	V = lhsd::lhsd_gen(U_T);


	size_t f;
	f = column(V, 1).size();
	double average = 0.0f;
	if (f != 0) {
		average = accumulate(column(V, 1).begin(), column(V, 1).end(), 0.0) / f;
	}
	
}






