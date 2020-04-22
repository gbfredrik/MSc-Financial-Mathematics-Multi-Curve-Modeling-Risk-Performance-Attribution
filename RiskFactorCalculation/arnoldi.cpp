#include "arnoldi.h"
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

std::tuple<matrix<double>, vector<double>> arnoldi::iram(matrix<double> const& input, int k) {
	
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
	SymGEigsSolver<double, SMALLEST_MAGN, DenseSymMatProd<double>, DenseCholesky<double>, GEIGS_CHOLESKY> geigs(&op, &Bop, k, 2.0 * k + 1.0);

	// Initialize and compute
	geigs.init();
	int nconv = geigs.compute();
	
	// Retrieve results
	if (geigs.info() == SUCCESSFUL){
		lambda = geigs.eigenvalues();
		E = geigs.eigenvectors();
	}
	
	// Return k eigenpairs as a tuple
	return std::make_tuple(matrixOperations::matrixXdToUblas(E), matrixOperations::vectorXdToUblas(lambda));
	
}


