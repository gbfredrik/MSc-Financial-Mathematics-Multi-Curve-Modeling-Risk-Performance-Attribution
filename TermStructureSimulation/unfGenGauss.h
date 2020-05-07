#ifndef UNFGENGAUSS
#define UNFGENGAUSS
#include "mex.h"
#include <boost/numeric/ublas/matrix.hpp>

class unfGenGauss {
private:

public:
	static boost::numeric::ublas::matrix<double> GC_sim(boost::numeric::ublas::matrix<double> const& E, int N);
};

#endif