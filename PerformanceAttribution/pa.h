#ifndef PA
#define PA
#include "mex.h"
#include <boost/numeric/ublas/matrix.hpp>


class pa {
private:

public:
	static void performanceAttribution(boost::numeric::ublas::vector<double> const& N, boost::numeric::ublas::vector<double> const& y,
		boost::numeric::ublas::vector<boost::numeric::ublas::matrix<double>> const& E, boost::numeric::ublas::vector<size_t> k,
		boost::numeric::ublas::vector<boost::numeric::ublas::vector<double>> const& floatCashFlows, 
		boost::numeric::ublas::vector<boost::numeric::ublas::vector<double>> const& fixCashFlows, 
		boost::numeric::ublas::vector<boost::numeric::ublas::matrix<double>> const& curveData,
		boost::numeric::ublas::vector<double> const& times, int startDate, int endDate, boost::numeric::ublas::vector<double>& NPV,
		boost::numeric::ublas::vector<double>& carry, boost::numeric::ublas::vector<double>& sumRiskFactors,
		boost::numeric::ublas::vector<double>& epsI, boost::numeric::ublas::vector<double>& epsA, 
		boost::numeric::ublas::vector<double>& epsP);
};

#endif
