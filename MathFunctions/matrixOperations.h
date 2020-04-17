#ifndef MATRIXOPERATIONS_H
#define MATRIXOPERATIONS_H

#include <boost/numeric/ublas/matrix.hpp>


class matrixOperations {

public:
	static boost::numeric::ublas::matrix<double> chol(boost::numeric::ublas::matrix<double> const& input);
};


#endif