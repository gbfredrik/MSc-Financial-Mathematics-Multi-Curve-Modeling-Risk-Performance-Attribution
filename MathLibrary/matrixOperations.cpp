#include "pch.h"
#include "matrixOperations.h"
#include <Eigen/Core>
#include <boost/numeric/ublas/matrix.hpp>

using namespace boost::numeric::ublas;
using namespace Eigen;


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

matrix<double> matrixOperations::matrixXdToUblas(MatrixXd const& xdMatrix) {

	size_t m = xdMatrix.innerSize();
	size_t n = xdMatrix.outerSize();

	matrix<double> uMatrix(m, n);


	for (size_t i = 0; i < m; i++) {
		for (size_t j = 0; j < n; j++) {
			uMatrix(i, j) = xdMatrix(i, j);
		}
	}

	return uMatrix;

}

MatrixXd matrixOperations::ublasToMatrixXd(matrix<double> const& uMatrix) {

	size_t m = uMatrix.size1();
	size_t n = uMatrix.size2();

	MatrixXd xdMatrix(m, n);

	for (size_t i = 0; i < m; i++) {
		for (size_t j = 0; j < n; j++) {
			xdMatrix(i, j) = uMatrix(i, j);
		}
	}

	return xdMatrix;
}

//Logarithm of a matrix
matrix<double> matrixOperations::matrixLog(matrix<double> const& input) {

	size_t rows = input.size1();
	size_t columns = input.size2();

	matrix<double> logMatrix(rows, columns);

	std::cout << "rows: " << rows << std::endl;
	std::cout << "columns: " << columns << std::endl;

	for (size_t i = 0; i < rows; ++i) {
		for (size_t j = 0; j < columns; ++j) {
			logMatrix(i, j) = log(input(i, j));
			std::cout << "logMatrix loop: " << logMatrix(i,j) << std::endl;
		}
	}
	std::cout << "logMatrix: " << logMatrix(1,1) << std::endl;
	return logMatrix;
}