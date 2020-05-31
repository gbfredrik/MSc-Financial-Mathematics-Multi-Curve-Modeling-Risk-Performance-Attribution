#pragma once

#include "mex.h"
#include <boost/numeric/ublas/matrix.hpp>
//extern boost::numeric::ublas::matrix<double> hehe;
class varRed {
private:
	static boost::numeric::ublas::matrix<double> lhsd_gen(boost::numeric::ublas::matrix<double> const& U);
	static boost::numeric::ublas::vector<double> rank(boost::numeric::ublas::vector<double> const& U);
	static boost::numeric::ublas::vector<double> sort(boost::numeric::ublas::vector<double> const& U);

public:
	static boost::numeric::ublas::matrix<double> redVariance(boost::numeric::ublas::matrix<double> const& U, std::string varRedType);
};
