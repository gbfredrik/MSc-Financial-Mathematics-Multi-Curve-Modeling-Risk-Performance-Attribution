#include "pch.h"
#include "matrixOperations.h"

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>

using namespace boost::numeric::ublas;

matrix<double> matrixOperations::chol(matrix<double> const& input) {
	size_t n = input.size1();
	matrix<double> L(n, n, 0);
	// Decomposing a matrix into lower triangular 
	for (size_t i = 0; i < n; i++) {
		for (size_t k = 0; k < (i + 1); k++) {
			double sum = 0;
			for (size_t j = 0; j < k; j++) {
				sum += L(i, j) * L(k, j);
			}
			L(i, k) = (i == k) ? sqrt(input(i, i) - sum) : (1.0 / L(k, k) * (input(i, k) - sum));
		}
	}
	return L;
}

matrix<double> matrixOperations::diff_matrix(matrix<double>& m_curves) { // Can this be "streamlined"? E.g. by directly returning a subset - a subset?
	matrix_range<matrix<double>> m_1(m_curves, range(0, m_curves.size1() - 1), range(0, m_curves.size2()));
	matrix_range<matrix<double>> m_2(m_curves, range(1, m_curves.size1()), range(0, m_curves.size2()));

	//matrix<double> m_diff(m_2 - m_1);
	return m_2 - m_1;
}
