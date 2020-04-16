#include <iostream>
#include <numeric>
#include <cmath>

//Boost packages for numeric represenation
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/io.hpp>

//Boost packages for statistics
#include <boost/math/tools/bivariate_statistics.hpp>
#include <boost/math/distributions/normal.hpp>

//Package to generate random variables
#include <random>

using namespace boost::numeric::ublas;
using namespace boost::math;

void GC_sim();
vector<double> gen_normal(double m, double s, int n);
matrix<double> gen_test(int rows, int cols);
double pearson_rho(vector<double> X, vector<double> Y);
matrix<double> corrm(matrix<double> input);
matrix<double> chol(matrix<double> input);

int main() {
	GC_sim();
}

// Controller method
void GC_sim() {

	int N = 2000;
	int m = 11000;
	int n = 6;
	// Generate a test matrix
	matrix<double> test(m, n);
	test = gen_test(m, n);

	// Calculate the correlation matrix
	matrix<double> corr(n, n);
	corr = corrm(test);

	// Perform Cholesky decomposition
	matrix<double> L(n, n);
	L = chol(corr);

	// Generate i.i.d. standard normal random variables
	matrix<double> X(N, n);
	for (int i = 0; i < n; i++) {
		column(X, i) = gen_normal(0.0, 1.0, N);
	}

	// Generate Correlated Gaussian samples
	matrix<double> Z(N, n);
	Z = prod(X, trans(L));

	// Return correlated uniformy distributed random variables
	matrix<double> U(N, n);
	normal s;
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < n; j++) {
			U(i, j) = cdf(s, Z(i, j));
		}
	}
	//std::cout << U << std::endl;
}

// Cholesky decomposition
matrix<double> chol(matrix<double> input) {
	size_t n = input.size1();
	matrix<double> L(n, n, 0);
	// Decomposing a matrix into lower triangular 
	for (int i = 0; i < n; i++) {
		for (int k = 0; k < (i + 1); k++) {
			double sum = 0;
			for (int j = 0; j < k; j++) {
				sum += L(i, j) * L(k, j);
			}
			L(i, k) = (i == k) ? sqrt(input(i, i) - sum) : (1.0 / L(k, k) * (input(i, k) - sum));
		}
	}
	return L;
}

// Calculate the Pearson correlation matrix, TODO: reduce the amount of calls to pearson_rho()
matrix<double> corrm(matrix<double> input) {
	size_t m = input.size1();
	size_t n = input.size2();
	matrix<double> corr(n, n);
	vector<double> X(m);
	vector<double> Y(m);

	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			if (i == j) {
				corr(i, j) = 1;
			}
			else {
				corr(i, j) = pearson_rho(column(input, i), column(input, j));
			}
		}
	}
	return corr;
}
// Calculate the Pearson correlation coefficient
double pearson_rho(vector<double> X, vector<double> Y) {
	double rho = 0.0;
	size_t m = X.size();
	double numerator = 0;
	double denomenator_a = 0;
	double denomenator_b = 0;

	//Calculate mean of input vectors
	double sum = std::accumulate(X.begin(), X.end(), 0.0);
	double X_hat = sum / X.size();
	sum = std::accumulate(Y.begin(), Y.end(), 0.0);
	double Y_hat = sum / Y.size();

	for (int i = 0; i < m; i++) {
		numerator = numerator + (X(i) - X_hat) * (Y(i) - Y_hat);
		denomenator_a = denomenator_a + pow(X(i) - X_hat, 2);
		denomenator_b = denomenator_b + pow(Y(i) - Y_hat, 2);
	}

	return numerator / (sqrt(denomenator_a) * sqrt(denomenator_b));
}

matrix<double> gen_test(int rows, int cols) {
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
vector<double> gen_normal(double m, double s, int N) {
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
