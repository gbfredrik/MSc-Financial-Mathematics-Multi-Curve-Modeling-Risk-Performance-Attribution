#include "FactorCalculation.h"

#include "../MathLibrary/matrixOperations.h"

#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/QR>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <Spectra/MatOp/DenseSymMatProd.h>
#include <Spectra/SymEigsSolver.h>

#include <iostream>
#include <Windows.h>
#include <WinUser.h>

using namespace boost::numeric;

bool FactorCalculation::iram(
	ublas::matrix<double> const& input, 
	int const k, 
	ublas::matrix<double>& m_E, 
	ublas::vector<double>& v_Lambda
) {
	// Utilizes the SymEigsSolver (IRAM) from the Spectra library
	using namespace Spectra;

	// Transform ublas matrix into Eigen::matrix
	Eigen::MatrixXd D{ matrixOperations::ublasToMatrixXd(input) };
	
	if (D.rows() < D.cols()) {
		Eigen::HouseholderQR<Eigen::MatrixXd> qr(D.transpose());

		// Define the positive definite RTR matrix
		Eigen::MatrixXd RTR{ qr.matrixQR().transpose() * qr.matrixQR() };
		Eigen::MatrixXd thinQ(Eigen::MatrixXd::Identity(D.cols(), D.rows()));
		thinQ = qr.householderQ() * thinQ;

		// Construct matrix operation objects
		DenseSymMatProd<double> op(RTR);

		// Construct eigen solver object, requesting the largest k eigenvalues in magnitude
		SymEigsSolver<double, LARGEST_ALGE, DenseSymMatProd<double>> eigs(&op, k, 2 * k + 1);
		
		// Initialize and compute
		eigs.init();
		int nconv{ eigs.compute() };

		// Retrieve results
		if (eigs.info() == SUCCESSFUL) {
			// Return eigenpairs
			m_E = matrixOperations::matrixXdToUblas(thinQ * eigs.eigenvectors());
			v_Lambda = matrixOperations::vectorXdToUblas(eigs.eigenvalues()); // TODO: Scale by n-1
		}
	} else {
		// Define the positive definite C matrix
		Eigen::MatrixXd C{ D.transpose() * D };

		// Construct matrix operation objects
		DenseSymMatProd<double> op(C);

		// Construct eigen solver object, requesting the largest k eigenvalues in magnitude
		SymEigsSolver<double, LARGEST_ALGE, DenseSymMatProd<double>> eigs(&op, k, 2 * k + 1);

		eigs.init();
		int nconv{ eigs.compute() };

		// Retrieve results
		if (eigs.info() == SUCCESSFUL) {
			// Return eigenpairs
			m_E = matrixOperations::matrixXdToUblas(eigs.eigenvectors());
			v_Lambda = matrixOperations::vectorXdToUblas(eigs.eigenvalues()); // TODO: Scale by n-1
		}
	}
	
	return v_Lambda.size() > 0;
}

bool FactorCalculation::eigen_bdcsvd(
	ublas::matrix<double> const& input, 
	int const k, 
	ublas::matrix<double>& m_E, 
	ublas::vector<double>& v_Lambda
) {
	// Utilizes BDCSVD from Eigen library
	int svd_opt{ Eigen::ComputeThinU | Eigen::ComputeThinV };

	Eigen::MatrixXd H{ matrixOperations::ublasToMatrixXd(input) };

	if (H.rows() < H.cols()) {
		Eigen::HouseholderQR<Eigen::MatrixXd> qr(H.transpose()); // FullPivHouseholderQR<Matrix<double, Dynamic, Size>> fpqr(A.rows(), A.cols());
		Eigen::MatrixXd thinQ(qr.householderQ() * Eigen::MatrixXd::Identity(H.cols(), H.rows()));

		//Eigen::MatrixXd thinQ(Eigen::MatrixXd::Identity(H.cols(), H.rows()));
		//thinQ = qr.householderQ() * thinQ;
		Eigen::MatrixXd RTR{ qr.matrixQR().transpose() * qr.matrixQR() };

		Eigen::BDCSVD<Eigen::MatrixXd> bdcsvd(RTR, svd_opt);

		// Return eigenpairs
		m_E = matrixOperations::matrixXdToUblas(thinQ * bdcsvd.matrixV().block(0, 0, bdcsvd.matrixV().rows(), k));
		v_Lambda = matrixOperations::vectorXdToUblas(bdcsvd.singularValues().head(k).array().square()); // TODO: Control square method
		// See: https://stackoverflow.com/questions/34373757/piece-wise-square-of-vector-piece-wise-product-of-two-vectors-in-c-eigen
	} else {
		Eigen::BDCSVD<Eigen::MatrixXd> bdcsvd(H.transpose() * H, svd_opt);
		m_E = matrixOperations::matrixXdToUblas(bdcsvd.matrixV().block(0, 0, bdcsvd.matrixV().rows(), k));
		v_Lambda = matrixOperations::vectorXdToUblas(bdcsvd.singularValues().head(k).array().square()); // TODO: Control square method
	}
		
	return v_Lambda.size() > 0;
}

ublas::matrix<double> FactorCalculation::compute_risk_factors(
	ublas::matrix<double> const& m_E_k, 
	ublas::matrix<double> const& m_delta_f
) {
	return prod(trans(m_E_k), trans(m_delta_f));
}

double FactorCalculation::smallest_eigval(ublas::matrix<double> const& input) {
	// Utilizes the SymEigsSolver (IRAM) from the Spectra library
	using namespace Spectra;
	int k{ 1 };
	// Transform ublas matrix into Eigen::matrix
	Eigen::MatrixXd C{ matrixOperations::ublasToMatrixXd(input) };

	// Construct matrix operation objects
	DenseSymMatProd<double> op(C);

	// Construct eigen solver object, requesting the largest k eigenvalues in magnitude
	SymEigsSolver<double, SMALLEST_ALGE, DenseSymMatProd<double>> eigs(&op, k, 2 * k + 1);

	eigs.init();
	int nconv{ eigs.compute() };

	// Retrieve results
	if (eigs.info() == SUCCESSFUL) {
		// Return smallest eigenvalue
		return eigs.eigenvalues()[0];
	}
	
    MessageBoxA(
        NULL,
        "Failed to find smallest eigenvalue.",
        "Status",
        MB_OK);
	return -1.0;
}

double FactorCalculation::eig_norm_error(
	ublas::matrix<double> const& m_A, 
	ublas::vector<double> const& v_x, 
	double const lambda
) {
	return norm_2(prod(m_A, v_x) - lambda * v_x);
}

ublas::vector<double> FactorCalculation::eig_all_norm_errors(
	ublas::matrix<double> const& m_A,
	ublas::matrix<double> const& m_x,
	ublas::vector<double> const& v_Lambda
) {
	size_t len{ v_Lambda.size() };
	ublas::vector<double> v_errors(len);

	for (size_t i{ 0 }; i < len; ++i) {
		v_errors(i) = eig_norm_error(m_A, column(m_x, i), v_Lambda(i));
	}

	return v_errors;
}

ublas::matrix<double> FactorCalculation::clean_data(ublas::matrix<double> const& m) {
	// TODO: Implement

	return ublas::matrix<double>();
}
