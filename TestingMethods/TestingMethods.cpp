//#include "../TermStructureSimulation/unfGenGauss.h"
//#include "../TermStructureSimulation/unfGenT.h"
//#include "../TermStructureSimulation/lhsd.h"
//#include "../MathLibrary/rvSim.h"
#include "../MathLibrary/matrixOperations.h"
#include "../Backtesting/backtesting.h"
//#include "../RiskFactorCalculation/FactorCalculation.h"

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
void test_likelihood();
void test_likelihood2();

int main() {

	//test_algo4_5();

	//test_matrix();
	test_likelihood();
	test_likelihood2();


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
void test_likelihood() {
	matrix<double> functions1(2, 3);
	matrix<double> functions2(2, 3);
	matrix<double> d_i(2, 3);
	vector<double> d(2);
	vector<double> sigma(2);
	vector<double> s(2);
	vector<double> p(2);
	vector<double> confidence_level(2);
	matrix<double> prices1(2,3);
	matrix<double> prices2(2,3);
	matrix<double> dR_i(2, 2);
	functions2(0, 0) = 1;
	functions2(0, 1) = 2;
	functions2(0, 2) = 3;
	functions2(1, 0) = 4;
	functions2(1, 1) = 5;
	functions2(1, 2) = 6;
	functions1(0, 0) = 7;
	functions1(0, 1) = 8;
	functions1(0, 2) = 9;
	functions1(1, 0) = 10;
	functions1(1, 1) = 11;
	functions1(1, 2) = 12;

	prices2(0, 0) = 1.1;
	prices2(0, 1) = 2.1;
	prices2(0, 2) = 3.5;
	prices2(1, 0) = 3.9;
	prices2(1, 1) = 5.1;
	prices2(1, 2) = 6.5;
	prices1(0, 0) = 7.1;
	prices1(0, 1) = 8.2;
	prices1(0, 2) = 9.3;
	prices1(1, 0) = 10.2;
	prices1(1, 1) = 10.9;
	prices1(1, 2) = 11.8;

	double N = 3.0;
	confidence_level(0)= 0.95;
	confidence_level(1) = 0.95;

	std::cout << "functions1: " << functions1 << std::endl;
	std::cout << "prices1: " << prices1 << std::endl;
	d_i = Likelihood::d_i(functions1, functions2);
	std::cout << "d_i: " << d_i << std::endl;
	d = Likelihood::d(d_i, N);
	std::cout << "d: " << d << std::endl;
	sigma = Likelihood::sigma(d,d_i, N);
	std::cout << "sigma: " << sigma << std::endl;
	s = Likelihood::s(sigma, N);
	std::cout << "s: " << s << std::endl;
	dR_i = Likelihood::dResidual_i(functions1, functions2, prices1, prices2);
	std::cout << "dR_i: " << dR_i << std::endl;
	p = Likelihood::probability(d, s, confidence_level);
	std::cout << "p: " << p << std::endl;
}

void test_likelihood2() {
	matrix<double> functions1(2, 3);
	matrix<double> functions2(2, 3);
	matrix<double> prices1(2, 3);
	matrix<double> prices2(2, 3);
	functions2(0, 0) = 1;
	functions2(0, 1) = 2;
	functions2(0, 2) = 3;
	functions2(1, 0) = 4;
	functions2(1, 1) = 5;
	functions2(1, 2) = 6;
	functions1(0, 0) = 7;
	functions1(0, 1) = 8;
	functions1(0, 2) = 9;
	functions1(1, 0) = 10;
	functions1(1, 1) = 11;
	functions1(1, 2) = 12;

	prices2(0, 0) = 1.1;
	prices2(0, 1) = 2.1;
	prices2(0, 2) = 3.5;
	prices2(1, 0) = 3.9;
	prices2(1, 1) = 5.1;
	prices2(1, 2) = 6.5;
	prices1(0, 0) = 7.1;
	prices1(0, 1) = 8.2;
	prices1(0, 2) = 9.3;
	prices1(1, 0) = 10.2;
	prices1(1, 1) = 10.9;
	prices1(1, 2) = 11.8;

	vector<double> p(2);
	vector<double> p2(2);

	p = Likelihood::likelihoodRatioTest(functions1, functions2);
	std::cout << "p2: " << p << std::endl;
	p2 = Likelihood::likelihoodRatioTestResidual(functions1, functions2, prices1, prices2);
	std::cout << "p3: " << p2 << std::endl;

}

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
