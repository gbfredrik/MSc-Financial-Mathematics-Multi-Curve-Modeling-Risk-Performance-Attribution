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
	
	// Create return vector and matrix
	Eigen::VectorXd lambda;
	Eigen::MatrixXd E;
	
	// Transform ublas matrix into Eigen::matrix
	Eigen::MatrixXd D = matrixOperations::ublasToMatrixXd(input);
	
	// Define the positive definite C matrix
	Eigen::MatrixXd C = D.transpose() * D;

	// Construct matrix operation objects
	DenseSymMatProd<double> op(D);
	DenseCholesky<double> Bop(C);

	// Construct eigen solver object, requesting the largest k eigenvalues in magnitude
	//SymGEigsSolver<double, LARGEST_MAGN, DenseSymMatProd<double>, DenseCholesky<double>, GEIGS_CHOLESKY> geigs(&op, &Bop, k, 2.0 * k + 1.0);

	// Initialize and compute
	//geigs.init();
	//int nconv = geigs.compute();
	
	// Retrieve results
	//if (geigs.info() == SUCCESSFUL){
	//	lambda = geigs.eigenvalues();
	//	E = geigs.eigenvectors();
	//	// Return k eigenpairs as a tuple
	//	m_E = matrixOperations::matrixXdToUblas(E);
	//	v_Lambda = matrixOperations::vectorXdToUblas(lambda);
	//	return true;
	//}
	
	return true;
}

bool FactorCalculation::eigen_bdcsvd(boost::numeric::ublas::matrix<double> const& input, boost::numeric::ublas::matrix<double>& m_E, boost::numeric::ublas::vector<double>& v_Lambda) {
	// Utilizes BDCSVD from Eigen library
	int svd_opt = ComputeThinU | ComputeThinV;

	MatrixXd H = matrixOperations::ublasToMatrixXd(input);

	//FullPivHouseholderQR<Matrix<double, Dynamic, Size>> fpqr(A.rows(), A.cols());

	//HouseholderQR<MatrixXd> qr(H.rows(), H.cols());
	HouseholderQR<MatrixXd> qr(H.transpose());

	MatrixXd thinQ(MatrixXd::Identity(H.cols(), H.rows()));
	//MatrixXd Q;
	//Q = qr.householderQ();
	thinQ = qr.householderQ() * thinQ;
	MatrixXd RR = qr.matrixQR().transpose() * qr.matrixQR();
	BDCSVD<MatrixXd> bdcsvd(RR, svd_opt);

	// Return eigenpairs
	m_E = matrixOperations::matrixXdToUblas(thinQ * bdcsvd.matrixV()); // OK?
	v_Lambda = matrixOperations::vectorXdToUblas(bdcsvd.singularValues().array().square()); // TODO: Control square method
	// See: https://stackoverflow.com/questions/34373757/piece-wise-square-of-vector-piece-wise-product-of-two-vectors-in-c-eigen
	return v_Lambda.size() > 0;
}
