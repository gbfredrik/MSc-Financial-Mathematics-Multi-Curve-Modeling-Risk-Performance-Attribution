#ifndef ALGO5
#define ALGO5

#include <boost/numeric/ublas/matrix.hpp>

using namespace boost::numeric::ublas;

class algo5 {
private:

public:
	static matrix<double> TC_sim(matrix<double> const& E, int N);
};

#endif
