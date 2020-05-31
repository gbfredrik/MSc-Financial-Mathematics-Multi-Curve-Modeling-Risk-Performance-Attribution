#pragma once

#include "mex.h"
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/optional/optional.hpp>
#include <boost/none_t.hpp>

class MultipleYieldSim {
private:
	static void simMultipleDaily(boost::numeric::ublas::vector<boost::numeric::ublas::matrix<double>> const& E, 
		boost::numeric::ublas::vector<double> const& fZero, boost::numeric::ublas::matrix<double> const& pi, 
		boost::numeric::ublas::vector<boost::numeric::ublas::matrix<double>> const& eps,
		int M, int N, boost::numeric::ublas::vector<boost::numeric::ublas::matrix<double>>& fRes,
		boost::numeric::ublas::vector<boost::numeric::ublas::matrix<double>>& histPrevSim);
	static void simSingleMultipleDaily(boost::numeric::ublas::vector<boost::numeric::ublas::matrix<double>> const& E,
		boost::numeric::ublas::vector<boost::numeric::ublas::matrix<double>> const& eps,
		int M, int N, boost::numeric::ublas::vector<boost::numeric::ublas::matrix<double>>& fRes,
		boost::numeric::ublas::vector<boost::numeric::ublas::matrix<double>>& histPrevSim,
		boost::numeric::ublas::vector<boost::numeric::ublas::matrix<double>>& histPrevPrevSim);
public:
	static void simMultipleFull(boost::numeric::ublas::vector<boost::numeric::ublas::matrix<double>> const& E, 
		boost::numeric::ublas::vector<boost::numeric::ublas::matrix<double>> const& rho, boost::numeric::ublas::vector<boost::numeric::ublas::vector<double>> const& mu,
		boost::numeric::ublas::vector<boost::numeric::ublas::vector<double>> const& omega, boost::numeric::ublas::vector<boost::numeric::ublas::vector<double>> const& alpha,
		boost::numeric::ublas::vector<boost::numeric::ublas::vector<double>> const& beta, boost::numeric::ublas::vector<boost::numeric::ublas::matrix<double>> const& hist,
		boost::numeric::ublas::vector<std::string> marginal, boost::numeric::ublas::vector<std::string> copula, boost::numeric::ublas::vector<std::string> varRedType, size_t d, size_t N, 
		boost::numeric::ublas::vector<boost::numeric::ublas::matrix<double>>& fRes, boost::optional<boost::numeric::ublas::vector<boost::numeric::ublas::vector<double>>> const& gamma,
		boost::optional<boost::numeric::ublas::vector<double>> const& kappa, boost::optional<boost::numeric::ublas::vector<boost::numeric::ublas::vector<double>>> const& xiHat,
		boost::optional<boost::numeric::ublas::vector<boost::numeric::ublas::vector<double>>> const& dfC, boost::optional<boost::numeric::ublas::vector<boost::numeric::ublas::vector<double>>> const& dfM);


};
