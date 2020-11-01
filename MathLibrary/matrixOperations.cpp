#include "pch.h"
#include "matrixOperations.h"

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>

#include <Eigen/Core>
#include <Eigen/Dense>

using namespace boost::numeric;

ublas::matrix<double> matrixOperations::chol(ublas::matrix<double> const& input) {
	size_t n{ input.size1() };
	ublas::matrix<double> L(n, n);

	for (size_t i{ 0 }; i < n; ++i) {
		for (size_t k{ 0 }; k < (i + 1); ++k) {
			double sum{ 0.0 };

			for (size_t j{ 0 }; j < k; ++j) {
				sum += L(i, j) * L(k, j);
			}
			L(i, k) = (i == k) ? sqrt(input(i, i) - sum) : (1.0 / L(k, k) * (input(i, k) - sum));
		}
	}

	return L;
}

ublas::matrix<double> matrixOperations::matrixXdToUblas(Eigen::MatrixXd const& xdMatrix) {
    size_t m{ static_cast<size_t>(xdMatrix.innerSize()) };
    size_t n{ static_cast<size_t>(xdMatrix.outerSize()) };
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
    size_t length{ static_cast<size_t>(xdVector.innerSize()) };
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

ublas::vector<double> matrixOperations::matrixToVector(ublas::matrix<double> const& mat) {
    size_t n{ mat.size1() };
    ublas::vector<double> resVec(n * n);
    size_t counter{ 0 };

    for (size_t j{ 0 }; j < n; ++j) {
        for (size_t i{ 0 }; i < n; ++i) {
            resVec(counter) = mat(i, j);
            ++counter;
        }
    }

    return resVec;
}

ublas::matrix<double> matrixOperations::vectorToMatrix(
    ublas::vector<double> const& vec, 
    size_t const n
) {
    //size_t n{ time_series.size2() };
    ublas::matrix<double> resMatrix(n, n);
    size_t counter{ 0 };

    for (size_t j{ 0 }; j < n; ++j) {
        for (size_t i{ 0 }; i < n; ++i) {
            resMatrix(i, j) = vec(counter);
            ++counter;
        }
    }

    return resMatrix;
}

ublas::matrix<double> matrixOperations::diff_matrix(ublas::matrix<double>& m_curves) {
    ublas::matrix_range<ublas::matrix<double>> m_1{
        m_curves,
        ublas::range(0, m_curves.size1() - 1),
        ublas::range(0, m_curves.size2())
    };
    ublas::matrix_range<ublas::matrix<double>> m_2{
        m_curves,
        ublas::range(1, m_curves.size1()),
        ublas::range(0, m_curves.size2())
    };

	return m_2 - m_1;
}

ublas::matrix<double> matrixOperations::center_matrix(ublas::matrix<double> const& diff_matrix) {
    size_t rows{ diff_matrix.size1() };
    size_t cols{ diff_matrix.size2() };
    ublas::matrix<double> centered_matrix(rows, cols);
	ublas::vector<double> column_averages(cols);

	for (size_t i{ 0 }; i < cols; ++i) {
		column_averages(i) = vector_average(column(diff_matrix, i));
	}

	for (size_t i{ 0 }; i < rows; ++i) {
		row(centered_matrix, i) = row(diff_matrix, i) - column_averages;
	}

	return centered_matrix;
}

double matrixOperations::vector_average(ublas::vector<double> const& vec) {
	return sum(vec) / vec.size();
}

double matrixOperations::vector_variance(ublas::vector<double> const& vec) {
    size_t length{ vec.size() };
    double variance{ 0.0 };
    double mean{ matrixOperations::vector_average(vec) };

    for (size_t i{ 0 }; i < length; ++i) {
		double var = pow(vec(i) - mean, 2);
		variance += pow(vec(i) - mean, 2);
	}

	return variance / (length - 1);
}

ublas::matrix<double> matrixOperations::matrixLog(ublas::matrix<double> const& input) {
	size_t rows{ input.size1() };
	size_t columns{ input.size2() };
	ublas::matrix<double> logMatrix(rows, columns);

	for (size_t i{ 0 }; i < rows; ++i) {
		for (size_t j{ 0 }; j < columns; ++j) {
			logMatrix(i, j) = log(input(i, j));
		}
	}

	return logMatrix;
}

ublas::vector<double> matrixOperations::vectorLog(ublas::vector<double> const& input) {
	size_t rows{ input.size() };
	ublas::vector<double> logVector(rows);

	for (size_t i{ 0 }; i < rows; ++i) {
		logVector(i) = log(input(i));
	}

	return logVector;
}

ublas::matrix<double> matrixOperations::matrix_inv(ublas::matrix<double> const& input) {
    return matrixXdToUblas(ublasToMatrixXd(input).inverse());
}

double matrixOperations::matrix_det(ublas::matrix<double> const& input) {
    return ublasToMatrixXd(input).determinant();
}

ublas::vector<double> matrixOperations::kron_prod_vec(ublas::vector<double> const& v1, ublas::vector<double> const& v2) {
    ublas::vector<double> vKron(v1.size() * v2.size());
    size_t counter{ 0 };
#pragma omp parallel for
    for (size_t i{ 0 }; i < v1.size(); ++i) {
        for (size_t j{ 0 }; j < v2.size(); ++j) {
            vKron(counter) = v1(i) * v2(j);
            ++counter;
        }
    }

    return vKron;
}

ublas::matrix<double> matrixOperations::abs(ublas::matrix<double> const& m) {
	ublas::matrix<double> out(m.size1(), m.size2());

	for (size_t r{ 0 }, rows{ m.size1() }; r < rows; ++r) {
		for (size_t c{ 0 }, cols{ m.size2() }; c < cols; ++c) {
			out(r, c) = std::abs(m(r, c));
		}
	}

	return out;
}

