#pragma once

//#include "mex.h"
#include <boost/numeric/ublas/matrix.hpp>
#include <Eigen/Core>

class matrixOperations {
public:
    static boost::numeric::ublas::matrix<double> chol(boost::numeric::ublas::matrix<double> const& input);

	static boost::numeric::ublas::matrix<double> matrixXdToUblas(Eigen::MatrixXd const& xdMatrix);
	static Eigen::MatrixXd ublasToMatrixXd(boost::numeric::ublas::matrix<double> const& uMatrix);

	static boost::numeric::ublas::vector<double> vectorXdToUblas(Eigen::VectorXd const& xdVector);
	static Eigen::VectorXd ublasToVectorXd(boost::numeric::ublas::vector<double> const& uVector);
    
	static boost::numeric::ublas::matrix<double> diff_matrix(boost::numeric::ublas::matrix<double>& m_curves);
	static boost::numeric::ublas::matrix<double> center_matrix(boost::numeric::ublas::matrix<double> const& diff_matrix);
	static double vector_average(boost::numeric::ublas::vector<double> const& vec);
  
	static boost::numeric::ublas::matrix<double> matrixLog(boost::numeric::ublas::matrix<double> const& input);
};
