#ifndef MULTIPLEYIELDSIM
#define MULTIPLEYIELDSIM
#include "mex.h"
#include <boost/numeric/ublas/matrix.hpp>

class MultipleYieldSim {
private:
	static void simMultipleDaily(boost::numeric::ublas::vector<boost::numeric::ublas::matrix<double>> const& E, 
		boost::numeric::ublas::vector<double> const& fZero, boost::numeric::ublas::matrix<double> const& pi, 
		boost::numeric::ublas::vector<boost::numeric::ublas::matrix<double>> const& eps,
		int M, int N, boost::numeric::ublas::vector<boost::numeric::ublas::matrix<double>>& fRes);
public:
	static void simMultipleFull(boost::numeric::ublas::vector<boost::numeric::ublas::matrix<double>> const& E, 
		boost::numeric::ublas::vector<double> const& fZero, boost::numeric::ublas::matrix<double> const& pi, 
		boost::numeric::ublas::vector<double> const& kappa, boost::numeric::ublas::matrix<double> const& xiHat, 
		boost::numeric::ublas::vector<double> const& omega, boost::numeric::ublas::vector<double> const& alpha,
		boost::numeric::ublas::vector<double> const& beta, boost::numeric::ublas::vector<boost::numeric::ublas::matrix<double>> const& hist,
		boost::numeric::ublas::vector<std::string> marginal, boost::numeric::ublas::vector<std::string> copula,
		boost::numeric::ublas::vector<std::string> varRedType, boost::numeric::ublas::matrix<double> const& mu,
		boost::numeric::ublas::matrix<double> const& df,
		int d, int N, boost::numeric::ublas::vector<boost::numeric::ublas::matrix<double>>& fRes);
};

#endif
