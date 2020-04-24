#include "FactorCalculation.h"
#include <Eigen/Core>
#include <boost/numeric/ublas/matrix.hpp>
#include "../MathLibrary/matrixOperations.h"
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
	SymGEigsSolver<double, LARGEST_MAGN, DenseSymMatProd<double>, DenseCholesky<double>, GEIGS_CHOLESKY> geigs(&op, &Bop, k, 2.0 * k + 1.0);

	// Initialize and compute
	geigs.init();
	int nconv = geigs.compute();
	
	// Retrieve results
	if (geigs.info() == SUCCESSFUL){
		lambda = geigs.eigenvalues();
		E = geigs.eigenvectors();
		// Return k eigenpairs as a tuple
		m_E = matrixOperations::matrixXdToUblas(E);
		v_Lambda = matrixOperations::vectorXdToUblas(lambda);
		return true;
	}
	
	return false;
}


