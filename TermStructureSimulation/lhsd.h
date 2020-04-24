#ifndef LHSD
#define LHSD

#include <boost/numeric/ublas/matrix.hpp>

class lhsd {
private:
	static boost::numeric::ublas::vector<double> rank(boost::numeric::ublas::vector<double> const& U);
	static boost::numeric::ublas::vector<double> sort(boost::numeric::ublas::vector<double> const& U);
public:
	static boost::numeric::ublas::matrix<double> lhsd_gen(boost::numeric::ublas::matrix<double> const& U);
};

#endif
