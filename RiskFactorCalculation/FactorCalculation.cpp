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
	Eigen::MatrixXd D = matrixOperations::ublasToMatrixXd(input);
	
	if (D.rows() < D.cols()) {
		Eigen::HouseholderQR<Eigen::MatrixXd> qr(D.transpose());

		// Define the positive definite RTR matrix
		Eigen::MatrixXd RTR = qr.matrixQR().transpose() * qr.matrixQR();
		Eigen::MatrixXd thinQ(Eigen::MatrixXd::Identity(D.cols(), D.rows()));
		thinQ = qr.householderQ() * thinQ;

		// Construct matrix operation objects
		DenseSymMatProd<double> op(RTR);

		// Construct eigen solver object, requesting the largest k eigenvalues in magnitude
		SymEigsSolver<double, LARGEST_ALGE, DenseSymMatProd<double>> eigs(&op, k, 2.0 * k + 1.0);

		// Initialize and compute
		eigs.init();
		int nconv = eigs.compute();

		// Retrieve results
		if (eigs.info() == SUCCESSFUL) {
			// Return eigenpairs
			m_E = matrixOperations::matrixXdToUblas(thinQ * eigs.eigenvectors());
			v_Lambda = matrixOperations::vectorXdToUblas(eigs.eigenvalues()); // TODO: Scale by n-1
		}
	} else {
		// Define the positive definite C matrix
		Eigen::MatrixXd C = D.transpose() * D;

		// Construct matrix operation objects
		DenseSymMatProd<double> op(C);

		// Construct eigen solver object, requesting the largest k eigenvalues in magnitude
		SymEigsSolver<double, LARGEST_ALGE, DenseSymMatProd<double>> eigs(&op, k, 2.0 * k + 1.0);

		eigs.init();
		int nconv = eigs.compute();

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
	boost::numeric::ublas::vector<double>& v_Lambda
) {
	// Utilizes BDCSVD from Eigen library
	int svd_opt = Eigen::ComputeThinU | Eigen::ComputeThinV;

	Eigen::MatrixXd H = matrixOperations::ublasToMatrixXd(input);	

	if (H.rows() < H.cols()) {
		Eigen::HouseholderQR<Eigen::MatrixXd> qr(H.transpose()); // FullPivHouseholderQR<Matrix<double, Dynamic, Size>> fpqr(A.rows(), A.cols());
		Eigen::MatrixXd thinQ(Eigen::MatrixXd::Identity(H.cols(), H.rows()));
		thinQ = qr.householderQ() * thinQ;
		Eigen::MatrixXd RTR = qr.matrixQR().transpose() * qr.matrixQR();

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
	return ublas::prod(ublas::trans(m_E_k), ublas::trans(m_delta_f));
}
