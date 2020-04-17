#ifndef ALGO4
#define ALGO4

#include <boost/numeric/ublas/matrix.hpp>

using namespace boost::numeric::ublas;

class algo4 {
private:

public:
	static matrix<double> GC_sim(matrix<double> const& E, int N);
};

#endif