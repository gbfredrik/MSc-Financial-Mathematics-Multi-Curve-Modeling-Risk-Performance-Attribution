#pragma once

#include <boost/numeric/ublas/matrix.hpp>


class statistics {
private:

public:
	static boost::numeric::ublas::matrix<double> corrm(boost::numeric::ublas::matrix<double> const& input);
	static double pearson_rho(boost::numeric::ublas::vector<double> const& X, boost::numeric::ublas::vector<double> const& Y);

};
