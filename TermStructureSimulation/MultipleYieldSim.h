#ifndef MULTIPLEYIELDSIM
#define MULTIPLEYIELDSIM
#include "mex.h"
#include <tuple>
#include <boost/numeric/ublas/matrix.hpp>

class MultipleYieldSim {
private:

public:
	static std::tuple<boost::numeric::ublas::matrix<double>, boost::numeric::ublas::matrix<double>> simMultiple(boost::numeric::ublas::vector<double> eps, boost::numeric::ublas::matrix<double> Ezero, boost::numeric::ublas::matrix<double> Etau, boost::numeric::ublas::vector<double> fZeroCurr, boost::numeric::ublas::vector<double> fTauCurr, int N);
	static void testSim(double x, double* y, double* z, mwSize n);
};

#endif
