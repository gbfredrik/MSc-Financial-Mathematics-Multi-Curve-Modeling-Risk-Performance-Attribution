#pragma once

//#include "mex.h"
#include <boost/numeric/ublas/matrix.hpp>

class statisticsOperations {
public:
	static boost::numeric::ublas::matrix<double> covm(boost::numeric::ublas::matrix<double> const& input);
	static boost::numeric::ublas::matrix<double> corrm(boost::numeric::ublas::matrix<double> const& input);
	static double pearson_rho(
		boost::numeric::ublas::vector<double> const& X,
		boost::numeric::ublas::vector<double> const& Y
	);
	static boost::numeric::ublas::vector<double> GARCH(
		boost::numeric::ublas::vector<double> const& omega, 
		boost::numeric::ublas::vector<double> const& alpha, 
		boost::numeric::ublas::vector<double> const& beta, 
		boost::numeric::ublas::vector<double> const& gamma, 
		boost::numeric::ublas::matrix<double> const& E,
		boost::numeric::ublas::matrix<double> const& fHist
	);
	static boost::numeric::ublas::vector<double> GARCH(
		boost::numeric::ublas::vector<double> const& omega,
		boost::numeric::ublas::vector<double> const& alpha,
		boost::numeric::ublas::vector<double> const& beta,
		boost::numeric::ublas::vector<double> const& gamma,
		boost::numeric::ublas::matrix<double> const& E,
		boost::numeric::ublas::vector<double> const& fPrev,
		boost::numeric::ublas::vector<double> const& fPrevPrev,
		boost::numeric::ublas::vector<double> const& sigmaPrev
	);
	static double invCDFNorm(double u, double mu, double sigma);
	static double invCDFT(double u, double df);

private:
	static double vectorMean(boost::numeric::ublas::vector<double> const& input);
};
