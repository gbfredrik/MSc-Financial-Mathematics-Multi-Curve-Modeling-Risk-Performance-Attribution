#ifndef UNFGENGAUSS
#define UNFGENGAUSS

#include <boost/numeric/ublas/matrix.hpp>

using namespace boost::numeric::ublas;

class unfGenGauss {
private:

public:
	static matrix<double> GC_sim(matrix<double> const& E, int N);
};

#endif