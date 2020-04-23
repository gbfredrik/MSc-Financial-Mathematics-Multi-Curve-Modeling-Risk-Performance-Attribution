#include "../TermStructureSimulation/unfGenGauss.h"
#include "../TermStructureSimulation/unfGenT.h"
#include "../TermStructureSimulation/lhsd.h"
#include "../MathLibrary/rvSim.h"
#include "../MathLibrary/matrixOperations.h"
#include "../RiskFactorCalculation/arnoldi.h"

#include <iostream>
#include <numeric>
#include <Eigen/Core>
#include <tuple>

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/range/counting_range.hpp>


using namespace boost::numeric::ublas;
using namespace boost::range;
using namespace Eigen;

//void test_algo4_5();
//void test_matrix();
//void test_arnoldi();

int main() {

	//test_algo4_5();
	//test_matrix();
	//test_arnoldi();

}

/*

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

*/

/*
void test_matrix() {

	size_t m = 5;
	size_t n = 5;

	matrix<double> u_test(m, n);
	MatrixXd x_test(m, n);

	MatrixXd A = MatrixXd::Random(m, n);
	matrix<double> B(m, n);
	B = rvSim::gen_test(m, n);

	x_test = matrixOperations::ublasToMatrixXd(B);
	u_test = matrixOperations::matrixXdToUblas(A);

	std::cout << B << std::endl;
	std::cout << x_test << std::endl;

	std::cout << A << std::endl;
	std::cout << u_test << std::endl;


}
*/

/*
void test_arnoldi() {

	size_t m = 20;
	size_t n = 20;

	int k = 6;

	matrix<double> D(m, n);
	matrix<double> E(m, n);
	vector<double> lambda(k);
	D = rvSim::gen_test(m, n);

	tie(E, lambda) = arnoldi::iram(D, k);

	std::cout << E << std::endl;
	std::cout << lambda << std::endl;
}
*/



