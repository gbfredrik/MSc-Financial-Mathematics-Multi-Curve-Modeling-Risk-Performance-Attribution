#include "statistics.h"

#include <numeric>

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>

using namespace boost::numeric::ublas;


// Calculate the Pearson correlation matrix, TODO: reduce the amount of calls to pearson_rho()
matrix<double> statistics::corrm(matrix<double> const& input) {
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
double statistics::pearson_rho(vector<double> const& X, vector<double> const& Y) {
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