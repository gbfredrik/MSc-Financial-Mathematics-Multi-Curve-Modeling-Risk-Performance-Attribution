#pragma once
#ifndef MULTIPLEYIELDSIM
#define MULTIPLEYIELDSIM
#include "mex.h"
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/optional/optional.hpp>

//extern boost::numeric::ublas::matrix<double> test;

class FXMultipleYieldSim {
private:
	static void FXsimMultipleDaily(boost::numeric::ublas::vector<boost::numeric::ublas::matrix<double>> const& E,
		boost::numeric::ublas::vector<boost::numeric::ublas::vector<double>> const& fStart, boost::numeric::ublas::matrix<double> const& pi,
		boost::numeric::ublas::vector<boost::numeric::ublas::matrix<double>> const& eps,
		int M, int N, boost::numeric::ublas::vector<boost::numeric::ublas::matrix<double>>& fRes);
public:
	static void FXsimMultipleFull(boost::numeric::ublas::vector<boost::numeric::ublas::matrix<double>> const& E,
		boost::numeric::ublas::vector<boost::numeric::ublas::matrix<double>> const& rho, boost::numeric::ublas::vector<boost::numeric::ublas::vector<double>> const& mu,
		boost::numeric::ublas::vector<boost::numeric::ublas::vector<double>> const& omega, boost::numeric::ublas::vector<boost::numeric::ublas::vector<double>> const& alpha,
		boost::numeric::ublas::vector<boost::numeric::ublas::vector<double>> const& beta, boost::numeric::ublas::vector<boost::numeric::ublas::matrix<double>> const& hist,
		boost::numeric::ublas::vector<std::string> marginal, boost::numeric::ublas::vector<std::string> copula, boost::numeric::ublas::vector<std::string> varRedType, size_t d, size_t N,
		boost::numeric::ublas::vector<boost::numeric::ublas::matrix<double>>& fRes, boost::optional<boost::numeric::ublas::vector<boost::numeric::ublas::vector<double>>> const& gamma,
		boost::optional<boost::numeric::ublas::vector<double>> const& kappa, boost::optional<boost::numeric::ublas::vector<boost::numeric::ublas::vector<double>>> const& xiHat,
		boost::optional<boost::numeric::ublas::vector<boost::numeric::ublas::vector<double>>> const& dfC, boost::optional<boost::numeric::ublas::vector<boost::numeric::ublas::vector<double>>> const& dfM);


};

#endif
