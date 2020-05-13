#include "FactorCalculation.h"

#include "../MathLibrary/matrixOperations.h"

#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/QR>
#include <boost/numeric/ublas/matrix.hpp>
#include <Spectra/MatOp/DenseSymMatProd.h>
//#include <Spectra/MatOp/SparseGenMatProd.h>
//#include <Spectra/MatOp/DenseCholesky.h>
//#include <Spectra/SymGEigsSolver.h>
#include <Spectra/SymEigsSolver.h>

#include <iostream>

using namespace boost::numeric::ublas;

bool FactorCalculation::iram(matrix<double> const& input, int const k, boost::numeric::ublas::matrix<double>& m_E, boost::numeric::ublas::vector<double>& v_Lambda) {
	using namespace Spectra;
	
	// Transform ublas matrix into Eigen::matrix
	Eigen::MatrixXd D = matrixOperations::ublasToMatrixXd(input);
	
	Eigen::HouseholderQR<Eigen::MatrixXd> qr(D.transpose());

	// Define the positive definite C matrix
	Eigen::MatrixXd RTR = qr.matrixQR().transpose() * qr.matrixQR();
	Eigen::MatrixXd thinQ(Eigen::MatrixXd::Identity(D.cols(), D.rows()));
	thinQ = qr.householderQ() * thinQ;


	// Construct matrix operation objects
	DenseSymMatProd<double> op(RTR);
	//DenseCholesky<double> Bop(RTR);

	// Construct eigen solver object, requesting the largest k eigenvalues in magnitude
	//SymGEigsSolver<double, LARGEST_MAGN, DenseSymMatProd<double>, DenseCholesky<double>, GEIGS_CHOLESKY> geigs(&op, &Bop, k, 2.0 * k + 1.0);
	SymEigsSolver<double, LARGEST_ALGE, DenseSymMatProd<double> > eigs(&op, k, 2.0 * k + 1.0);

	// Initialize and compute
	eigs.init();
	int nconv = eigs.compute();
	
	// Retrieve results
	if (eigs.info() == SUCCESSFUL){
		// Return eigenpairs
		m_E = matrixOperations::matrixXdToUblas(thinQ * eigs.eigenvectors());
		v_Lambda = matrixOperations::vectorXdToUblas(eigs.eigenvalues()); // TODO: Scale by n-1
		return true;
	}
	
	return true;
}

bool FactorCalculation::eigen_bdcsvd(boost::numeric::ublas::matrix<double> const& input, int const k, boost::numeric::ublas::matrix<double>& m_E, boost::numeric::ublas::vector<double>& v_Lambda) {
	// Utilizes BDCSVD from Eigen library
	int svd_opt = Eigen::ComputeThinU | Eigen::ComputeThinV;

	Eigen::MatrixXd H = matrixOperations::ublasToMatrixXd(input);	
	if (H.rows() < H.cols()) {
		Eigen::HouseholderQR<Eigen::MatrixXd> qr(H.transpose());//FullPivHouseholderQR<Matrix<double, Dynamic, Size>> fpqr(A.rows(), A.cols());
	
		Eigen::MatrixXd thinQ(Eigen::MatrixXd::Identity(H.cols(), H.rows()));
		thinQ = qr.householderQ() * thinQ;
		Eigen::MatrixXd RTR = qr.matrixQR().transpose() * qr.matrixQR();

		Eigen::BDCSVD<Eigen::MatrixXd> bdcsvd(RTR, svd_opt);
		// Return eigenpairs
		m_E = matrixOperations::matrixXdToUblas(thinQ * bdcsvd.matrixV()); // OK?
		v_Lambda = matrixOperations::vectorXdToUblas(bdcsvd.singularValues().array().square()); // TODO: Control square method
		// See: https://stackoverflow.com/questions/34373757/piece-wise-square-of-vector-piece-wise-product-of-two-vectors-in-c-eigen
	} else {
		Eigen::BDCSVD<Eigen::MatrixXd> bdcsvd(H.transpose() * H, svd_opt);
		m_E = matrixOperations::matrixXdToUblas(bdcsvd.matrixV()); // OK?
		v_Lambda = matrixOperations::vectorXdToUblas(bdcsvd.singularValues().array().square()); // TODO: Control square method
	}
		
	return v_Lambda.size() > 0; // Ändra så bara k=6 egenpar returneras!
}
