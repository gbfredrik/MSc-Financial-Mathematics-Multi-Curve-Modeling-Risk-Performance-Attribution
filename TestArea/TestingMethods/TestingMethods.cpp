#include "../../CurveLibrary/sample_handler.h"
#include "../MathLibrary/rvSim.h"
#include "../MathLibrary/matrixOperations.h"
#include "../RiskMeasurement/backtesting.h"

#include <iostream>
//#include <numeric>
//#include <Eigen/Core>
//#include <tuple>

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/range/counting_range.hpp>

using namespace boost::numeric::ublas;
using namespace boost::range;

void test_likelihood();
void test_likelihood2();
void test_read();

int main() {
	//test_likelihood();
	//test_likelihood2();
	test_read();
}

void test_read() {
	matrix<double> matrix1;
	matrix1 = read_csv_matrix("FX_SEK_base.csv");
	if (matrix1.size1() == 0){
		std::cout << "Failed reading file! "  << std::endl;
	}
	else {
		std::cout << "Succesfull reading! " << std::endl;
		std::cout << "Matrix: " << matrix1 << std::endl;
	}

	matrix1.resize(matrix1.size1(), 5);
	std::cout << "\n Matrix: " << matrix1 << std::endl;
}

/*
void test_likelihood() {
	vector<double> functions1(2);
	vector<double> functions2(2);
	vector<double> d_i(2);
	double d;
	double sigma;
	double s;
	double p;
	//vector<double> confidence_level(2);
	vector<double> prices1(2);
	vector<double> prices2(2);
	vector<double> dR_i(2);
	functions2(0) = 1;
	functions2(1) = 6;
	functions1(0) = 7;
	functions1(1) = 12;

	prices2(0) = 1.1;
	prices2(1) = 6.5;
	prices1(0) = 7.1;
	prices1(1) = 11.8;

	int N = 2;
	double confidence_level = 0.95;
	//confidence_level(1) = 0.95;

	std::cout << "functions1: " << functions1 << std::endl;
	std::cout << "prices1: " << prices1 << std::endl;
	d_i = Likelihood::d_i(functions1, functions2);
	std::cout << "d_i: " << d_i << std::endl;
	d = Likelihood::d(d_i, N);
	std::cout << "d: " << d << std::endl;
	sigma = Likelihood::sigma(d,d_i, N);
	std::cout << "sigma: " << sigma << std::endl;
	s = Likelihood::s(sigma, N);
	std::cout << "s: " << sigma / sqrt(N) << std::endl;
	dR_i = Likelihood::dResidual_i(functions1, functions2, prices1, prices2);
	std::cout << "dR_i: " << dR_i << std::endl;
	p = Likelihood::is_better(d, s, confidence_level);
	std::cout << "is_better: " << p << std::endl;
}

void test_likelihood2() {
	vector<double> functions1(2);
	vector<double> functions2(2);
	vector<double> prices1(2);
	vector<double> prices2(2);
	functions2(0) = 1;
	functions2(1) = 2;
	functions1(0) = 4;
	functions1(1) = 5;

	prices2(0) = 1.1;
	prices2(1) = 2.5;
	prices1(0) = 4.1;
	prices1(1) = 5.3;

	int isBetter;
	int isBetter2;
	
	double confidence_level = 0.95;
	std::cout << "------------------ " << std::endl;
	isBetter = Likelihood::likelihood_ratio_test(functions1, functions2, confidence_level);
	std::cout << "is_better2: " << isBetter << std::endl;
	std::cout << "------------------ " << std::endl;
	isBetter2 = Likelihood::likelihood_ratio_test_residual(functions1, functions2, prices1, prices2, confidence_level);
	std::cout << "is_better3: " << isBetter2 << std::endl;

	vector<double> x_sim(3);
	matrix<double> x_simM(3,3);
	vector<double> x_realizedV(3);
	vector<double> fV(3);
	double f;
	double x_realized = 0.011;
	x_sim(0) = 0.01;
	x_sim(1) = 0.012;
	x_sim(2) = 0.009;

	x_simM(0,0) = 0.01;
	x_simM(1,0) = 0.012;
	x_simM(2,0) = 0.009;
	x_simM(0, 1) = 0.023;
	x_simM(1, 1) = 0.022;
	x_simM(2, 1) = 0.021;
	x_simM(0, 2) = 0.033;
	x_simM(1, 2) = 0.032;
	x_simM(2, 2) = 0.031;

	x_realizedV(0) = 0.011;
	x_realizedV(1) = 0.019;
	x_realizedV(2) = 0.031;
	std::cout << "------------------ " << std::endl;
	f = KernelDensity::kde(x_sim, x_realized);
	std::cout << "f: " << f << std::endl;

	//fV = KernelDensity::kde_multi(x_simM, x_realizedV);
	//std::cout << "fV: " << fV << std::endl;

	//vector<double> x_simMultiTest = read_csv_vector("KernelRnd1.csv");
    //std::cout << "x_simMultiTest: " << x_simMultiTest << std::endl;
	//vector<double> x = read_csv_vector("KernelX.csv");
	//std::cout << "x: " << x << std::endl;
	//vector<double> fMultiTest = KernelDensity::kde_multi(x_simMultiTest, x);
	//std::cout << "fMultiTest: " << fMultiTest << std::endl;
}
*/