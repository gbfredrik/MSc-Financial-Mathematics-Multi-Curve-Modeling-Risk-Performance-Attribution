#include "pch.h"
#include "matrixOperations.h"

#include <boost/numeric/ublas/matrix.hpp>

using namespace boost::numeric::ublas;

matrix<double> matrixOperations::chol(matrix<double> const& input) {
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

