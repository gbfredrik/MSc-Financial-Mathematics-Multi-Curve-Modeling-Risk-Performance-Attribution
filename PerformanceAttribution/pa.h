#ifndef PA
#define PA
#include "mex.h"
#include <boost/numeric/ublas/matrix.hpp>


class pa {
private:
	static void intMatrix(boost::numeric::ublas::matrix<double>& A);
	static void seta(boost::numeric::ublas::matrix<double> const& A, boost::numeric::ublas::vector<boost::numeric::ublas::matrix<double>> const& E, 
		int kZero, int kPi, boost::numeric::ublas::matrix<double>& a, std::string type);
	static double irsPriceRiskFactor(double N, double y, boost::numeric::ublas::vector<int> const& floatCashFlows,
		boost::numeric::ublas::vector<int> const& fixCashFlows, boost::numeric::ublas::matrix<double> const& aZero,
		boost::numeric::ublas::matrix<double> aPi, boost::numeric::ublas::vector<double> const& deltaTj, boost::numeric::ublas::vector<double> const& XiZero,
		boost::numeric::ublas::vector<double> XiPi);
	static double irsPrice(double N, double y, boost::numeric::ublas::vector<int> const& floatCashFlows,
		boost::numeric::ublas::vector<int> const& fixCashFlows, boost::numeric::ublas::vector<double> const& deltaTj,
		boost::numeric::ublas::vector<double> const& r, boost::numeric::ublas::vector<double> const& pi);
	static void setXiBar(boost::numeric::ublas::vector<boost::numeric::ublas::matrix<double>> const& E, int n,
		boost::numeric::ublas::vector<double> const& f, boost::numeric::ublas::vector<double> const& pi,
		boost::numeric::ublas::vector<double>& XiBar, std::string type);
	static void setdXi(boost::numeric::ublas::vector<boost::numeric::ublas::matrix<double>> const& E, int kZero, int kPi,
		boost::numeric::ublas::vector<double> const& f, boost::numeric::ublas::vector<double> const& pi,
		boost::numeric::ublas::vector<double> const& fPrev, boost::numeric::ublas::vector<double> const& piPrev,
		boost::numeric::ublas::vector<double> & XiBar, std::string type);
	static boost::numeric::ublas::matrix<double> hess(double N, double y, boost::numeric::ublas::vector<int> floatCashFlows,
		boost::numeric::ublas::vector<int> fixCashFlows, boost::numeric::ublas::matrix<double> aZero,
		boost::numeric::ublas::matrix<double> aPi, boost::numeric::ublas::vector<double> deltaTj, boost::numeric::ublas::vector<double> r,
		boost::numeric::ublas::vector<double> pi);
	static boost::numeric::ublas::vector<double> grad(double N, double y, boost::numeric::ublas::vector<int> floatCashFlows,
		boost::numeric::ublas::vector<int> fixCashFlows, boost::numeric::ublas::matrix<double> aZero,
		boost::numeric::ublas::matrix<double> aPi, boost::numeric::ublas::vector<double> deltaTj,
		boost::numeric::ublas::vector<double> r, boost::numeric::ublas::vector<double> pi);
public:
	static void performanceAttribution(boost::numeric::ublas::vector<double> const& N, boost::numeric::ublas::vector<double> const& y,
		boost::numeric::ublas::vector<boost::numeric::ublas::matrix<double>> const& E, boost::numeric::ublas::vector<int> k,
		boost::numeric::ublas::vector<boost::numeric::ublas::vector<int>> floatCashFlows, 
		boost::numeric::ublas::vector<boost::numeric::ublas::vector<int>> fixCashFlows, 
		boost::numeric::ublas::vector<boost::numeric::ublas::matrix<double>> const& curveData,
		boost::numeric::ublas::vector<int> const& times, boost::numeric::ublas::vector<int> const& startDate,
		boost::numeric::ublas::vector<int> const& endDate, boost::numeric::ublas::vector<boost::numeric::ublas::vector<double>>& NPV,
		boost::numeric::ublas::vector<boost::numeric::ublas::vector<double>>& carry, boost::numeric::ublas::vector<boost::numeric::ublas::vector<double>>& sumRiskFactors,
		boost::numeric::ublas::vector<boost::numeric::ublas::vector<double>>& epsI, boost::numeric::ublas::vector<boost::numeric::ublas::vector<double>>& epsA,
		boost::numeric::ublas::vector<boost::numeric::ublas::vector<double>>& epsP);
};
#endif
