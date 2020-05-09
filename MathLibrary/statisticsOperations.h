#ifndef STATISTICSOPERATIONS
#define STATISTICSOPERATIONS

#include "mex.h"
#include <boost/numeric/ublas/matrix.hpp>


class statisticsOperations {
private:
	static double vectorMean(boost::numeric::ublas::vector<double> const& input);
public:
	static boost::numeric::ublas::matrix<double> covm(boost::numeric::ublas::matrix<double> const& input);
	static boost::numeric::ublas::matrix<double> corrm(boost::numeric::ublas::matrix<double> const& input);
	static double pearson_rho(boost::numeric::ublas::vector<double> const& X, boost::numeric::ublas::vector<double> const& Y);
	static boost::numeric::ublas::vector<double> GARCH(boost::numeric::ublas::vector<double> omega, boost::numeric::ublas::vector<double> alpha, boost::numeric::ublas::vector<double> beta, boost::numeric::ublas::matrix<double> E, boost::numeric::ublas::matrix<double> fHist);
	static boost::numeric::ublas::vector<double> GARCH(boost::numeric::ublas::vector<double> omega, boost::numeric::ublas::vector<double> alpha, boost::numeric::ublas::vector<double> beta, boost::numeric::ublas::matrix<double> E, boost::numeric::ublas::vector<double> ft1, boost::numeric::ublas::vector<double> ft2, boost::numeric::ublas::vector<double> sigmat1);

};


#endif