#ifndef UNFGENT
#define UNFGENT

#include <boost/numeric/ublas/matrix.hpp>

class unfGenT {
private:

public:
	static boost::numeric::ublas::matrix<double> TC_sim(boost::numeric::ublas::matrix<double> const& E, int N);
};

#endif

