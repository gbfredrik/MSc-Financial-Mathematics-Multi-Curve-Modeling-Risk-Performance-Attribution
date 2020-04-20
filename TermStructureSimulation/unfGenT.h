#ifndef UNFGENT
#define UNFGENT

#include <boost/numeric/ublas/matrix.hpp>

using namespace boost::numeric::ublas;

class unfGenT {
private:

public:
	static matrix<double> TC_sim(matrix<double> const& E, int N);
};

#endif

