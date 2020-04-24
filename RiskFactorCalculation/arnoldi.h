#ifndef ARNOLDI
#define ARNOLDI

#include <boost/numeric/ublas/matrix.hpp>
#include <Eigen/Core>


class arnoldi {
private:

public:
	static boost::numeric::ublas::matrix<double> matrixXdToUblas(Eigen::MatrixXd xdMatrix);
	static Eigen::MatrixXd ublasToMatrixXd(boost::numeric::ublas::matrix<double> uMatrix);
};

#endif