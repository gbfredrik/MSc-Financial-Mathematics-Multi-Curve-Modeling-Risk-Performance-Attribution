#include "FactorCalculation.h"

#include "../MathLibrary/matrixOperations.h"

#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/QR>
#include <boost/numeric/ublas/matrix.hpp>
#include <Spectra/MatOp/SparseGenMatProd.h>
#include <Spectra/MatOp/DenseCholesky.h>
#include <Spectra/SymGEigsSolver.h>
#include <Spectra/MatOp/DenseSymMatProd.h>

#include <iostream>

using namespace boost::numeric::ublas;
using namespace Spectra;
using namespace Eigen;

bool FactorCalculation::iram(matrix<double> const& input, int k, boost::numeric::ublas::matrix<double>& m_E, boost::numeric::ublas::vector<double>& v_Lambda) {
	// Transform ublas matrix into Eigen::matrix
	Eigen::MatrixXd D = matrixOperations::ublasToMatrixXd(input);
	
	HouseholderQR<MatrixXd> qr(D.transpose());

	// Define the positive definite C matrix
	MatrixXd RTR = qr.matrixQR().transpose() * qr.matrixQR();
	MatrixXd thinQ(MatrixXd::Identity(D.cols(), D.rows()));
	thinQ = qr.householderQ() * thinQ;


	// Construct matrix operation objects
	DenseSymMatProd<double> op(D);
	DenseCholesky<double> Bop(RTR);

	// Construct eigen solver object, requesting the largest k eigenvalues in magnitude
	SymGEigsSolver<double, LARGEST_MAGN, DenseSymMatProd<double>, DenseCholesky<double>, GEIGS_CHOLESKY> geigs(&op, &Bop, k, 2.0 * k + 1.0);

	// Initialize and compute
	geigs.init();
	int nconv = geigs.compute();
	
	// Retrieve results
	if (geigs.info() == SUCCESSFUL){
		// Return eigenpairs
		m_E = matrixOperations::matrixXdToUblas(thinQ * geigs.eigenvectors());
		v_Lambda = matrixOperations::vectorXdToUblas(geigs.eigenvalues()); // TODO: Scale by n-1
		return true;
	}
	
	return true;
}

bool FactorCalculation::eigen_bdcsvd(boost::numeric::ublas::matrix<double> const& input, boost::numeric::ublas::matrix<double>& m_E, boost::numeric::ublas::vector<double>& v_Lambda) {
	// Utilizes BDCSVD from Eigen library
	int svd_opt = ComputeThinU | ComputeThinV;

	MatrixXd H = matrixOperations::ublasToMatrixXd(input);

	//FullPivHouseholderQR<Matrix<double, Dynamic, Size>> fpqr(A.rows(), A.cols());
	HouseholderQR<MatrixXd> qr(H.transpose());

	MatrixXd thinQ(MatrixXd::Identity(H.cols(), H.rows()));
	thinQ = qr.householderQ() * thinQ;
	MatrixXd RTR = qr.matrixQR().transpose() * qr.matrixQR();
	BDCSVD<MatrixXd> bdcsvd(RTR, svd_opt);

	// Return eigenpairs
	m_E = matrixOperations::matrixXdToUblas(thinQ * bdcsvd.matrixV()); // OK?
	v_Lambda = matrixOperations::vectorXdToUblas(bdcsvd.singularValues().array().square()); // TODO: Control square method
	// See: https://stackoverflow.com/questions/34373757/piece-wise-square-of-vector-piece-wise-product-of-two-vectors-in-c-eigen
	return v_Lambda.size() > 0;
}
