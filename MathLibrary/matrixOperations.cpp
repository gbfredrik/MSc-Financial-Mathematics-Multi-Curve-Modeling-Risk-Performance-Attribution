#include "pch.h"
#include "matrixOperations.h"

//#include "mex.h"
#include <Eigen/Core>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>

#include <fstream> 

using namespace boost::numeric;

ublas::matrix<double> matrixOperations::chol(ublas::matrix<double> const& input) {
	size_t n{ input.size1() };
	ublas::matrix<double> L(n, n);

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

ublas::matrix<double> matrixOperations::matrixXdToUblas(Eigen::MatrixXd const& xdMatrix) {
	size_t m = xdMatrix.innerSize(); // Todo: Kontrollera typ
	size_t n = xdMatrix.outerSize(); // Dito
	ublas::matrix<double> uMatrix(m, n);

	for (size_t i{ 0 }; i < m; ++i) {
		for (size_t j{ 0 }; j < n; ++j) {
			uMatrix(i, j) = xdMatrix(i, j);
		}
	}

	return uMatrix;
}

Eigen::MatrixXd matrixOperations::ublasToMatrixXd(ublas::matrix<double> const& uMatrix) {
	size_t m{ uMatrix.size1() };
	size_t n{ uMatrix.size2() };
	Eigen::MatrixXd xdMatrix(m, n);

	for (size_t i{ 0 }; i < m; ++i) {
		for (size_t j{ 0 }; j < n; ++j) {
			xdMatrix(i, j) = uMatrix(i, j);
		}
	}

	return xdMatrix;
}

ublas::vector<double> matrixOperations::vectorXdToUblas(Eigen::VectorXd const& xdVector) {
	size_t length = xdVector.innerSize(); // Todo: Kontrollera typ
	ublas::vector<double> uVector(length);

	for (size_t i{ 0 }; i < length; ++i) {
		uVector(i) = xdVector(i);
	}

	return uVector;
}

Eigen::VectorXd matrixOperations::ublasToVectorXd(ublas::vector<double> const& uVector) {
	size_t length{ uVector.size() };
	Eigen::VectorXd xdVector(length);

	for (size_t i{ 0 }; i < length; ++i) {
		xdVector(i) = uVector(i);
	}

	return xdVector;
}

ublas::matrix<double> matrixOperations::diff_matrix(ublas::matrix<double>& m_curves) { // Can this be "streamlined"? E.g. by directly returning a subset - a subset?
	ublas::matrix_range<ublas::matrix<double>> m_1(
		m_curves, 
		ublas::range(0, m_curves.size1() - 1), 
		ublas::range(0, m_curves.size2())
	);
	ublas::matrix_range<ublas::matrix<double>> m_2(
		m_curves, 
		ublas::range(1, m_curves.size1()), 
		ublas::range(0, m_curves.size2())
	);

	return m_2 - m_1;
}

ublas::matrix<double> matrixOperations::center_matrix(ublas::matrix<double> const& diff_matrix) {
	ublas::matrix<double> centered_matrix(diff_matrix.size1(), diff_matrix.size2());
	ublas::vector<double> column_averages(diff_matrix.size2());

	for (size_t i{ 0 }, cols{ diff_matrix.size2() }; i < cols; ++i) {
		column_averages(i) = vector_average(column(diff_matrix, i));
	}

	for (size_t i{ 0 }, rows{ diff_matrix.size1() }; i < rows; ++i) {
		row(centered_matrix, i) = row(diff_matrix, i) - column_averages;
	}

	return centered_matrix;
}

double matrixOperations::vector_average(ublas::vector<double> const& vec) {
	return sum(vec) / vec.size();
}

//Logarithm of a matrix
ublas::matrix<double> matrixOperations::matrixLog(ublas::matrix<double> const& input) {
	size_t rows{ input.size1() };
	size_t columns{ input.size2() };

	ublas::matrix<double> logMatrix(rows, columns);

	//std::cout << "rows: " << rows << std::endl;
	//std::cout << "columns: " << columns << std::endl;

	for (size_t i{ 0 }; i < rows; ++i) {
		for (size_t j{ 0 }; j < columns; ++j) {
			logMatrix(i, j) = log(input(i, j));
			//std::cout << "logMatrix loop: " << logMatrix(i,j) << std::endl;
		}
	}
	//std::cout << "logMatrix: " << logMatrix(1,1) << std::endl;

	return logMatrix;
}
