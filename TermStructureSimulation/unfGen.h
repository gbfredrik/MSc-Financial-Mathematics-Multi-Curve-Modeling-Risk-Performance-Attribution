#pragma once

#include "mex.h"
#include <boost/numeric/ublas/matrix.hpp>

class unfGen {
private:
	static boost::numeric::ublas::matrix<double> GC_sim(boost::numeric::ublas::matrix<double> const& rho, int N);
	static boost::numeric::ublas::matrix<double> TC_sim(boost::numeric::ublas::matrix<double> const& rho, int N,
		boost::numeric::ublas::vector<double> const& df);

public:
	static boost::numeric::ublas::matrix<double> genU(boost::numeric::ublas::matrix<double> const& rho, int N, 
		std::string copula, boost::numeric::ublas::vector<double> const& df);
};
