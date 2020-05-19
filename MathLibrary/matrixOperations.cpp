#include "pch.h"
#include "matrixOperations.h"

//#include "mex.h"
#include <Eigen/Core>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/io.hpp>// Temp: remove later

#include <fstream> 

using namespace boost::numeric::ublas;
using namespace Eigen;

matrix<double> matrixOperations::chol(matrix<double> const& input) {
	size_t n{ input.size1() };
	matrix<double> L(n, n);

	// Decomposing a matrix into lower triangular 
	for (size_t i{ 0 }; i < n; ++i) {
		for (size_t k{ 0 }; k < (i + 1); ++k) {
			double sum{ 0 };
			for (size_t j{ 0 }; j < k; ++j) {
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


	for (size_t i{ 0 }; i < m; ++i) {
		for (size_t j{ 0 }; j < n; ++j) {
			uMatrix(i, j) = xdMatrix(i, j);
		}
	}

	return uMatrix;
}

MatrixXd matrixOperations::ublasToMatrixXd(matrix<double> const& uMatrix) {
	size_t m{ uMatrix.size1() };
	size_t n{ uMatrix.size2() };
	MatrixXd xdMatrix(m, n);

	for (size_t i = 0; i < m; i++) {
		for (size_t j = 0; j < n; j++) {
			xdMatrix(i, j) = uMatrix(i, j);
		}
	}

	return xdMatrix;
}

vector<double> matrixOperations::vectorXdToUblas(VectorXd const& xdVector) {
	size_t length = xdVector.innerSize();
	vector<double> uVector(length);

	for (size_t i = 0; i < length; i++) {
		uVector(i) = xdVector(i);
	}

	return uVector;
}

VectorXd matrixOperations::ublasToVectorXd(vector<double> const& uVector) {
	size_t length = uVector.size();
	VectorXd xdVector(length);

	for (size_t i = 0; i < length; i++) {
		xdVector(i) = uVector(i);
	}

	return xdVector;
}

matrix<double> matrixOperations::diff_matrix(matrix<double>& m_curves) { // Can this be "streamlined"? E.g. by directly returning a subset - a subset?
	matrix_range<matrix<double>> m_1(m_curves, range(0, m_curves.size1() - 1), range(0, m_curves.size2()));
	matrix_range<matrix<double>> m_2(m_curves, range(1, m_curves.size1()), range(0, m_curves.size2()));

	//matrix<double> m_diff(m_2 - m_1);
	return m_2 - m_1;
}

double matrixOperations::column_average(vector<double> const& vec) {
	return sum(vec) / vec.size();
}

matrix<double> matrixOperations::center_matrix(matrix<double> const& diff_matrix) {
	matrix<double> centered_matrix(diff_matrix.size1(), diff_matrix.size2());
	vector<double> column_averages(diff_matrix.size2());

	for (size_t i = 0; i < diff_matrix.size2(); ++i) {
		column_averages(i) = column_average(column(diff_matrix, i));
	}

	for (size_t i = 0; i < diff_matrix.size1(); ++i) {
		row(centered_matrix, i) = row(diff_matrix, i) - column_averages;
	}

	return centered_matrix;
}
