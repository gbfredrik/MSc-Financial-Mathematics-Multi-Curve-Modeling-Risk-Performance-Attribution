#ifndef BACKTESTING
#define BACKTESTING

#include <boost/numeric/ublas/matrix.hpp>

class Likelihood {
	//private:
public:
	static boost::numeric::ublas::vector<double> d_i(boost::numeric::ublas::vector<double> const& vector1, boost::numeric::ublas::vector<double> const& vector2);
	static double d(boost::numeric::ublas::vector<double> const& d_i, int N);
	static double s(double sigma, int N);
	static double sigma(double d, boost::numeric::ublas::vector<double> const& d_i, int N);
	static int isBetter(double d, double s, double confidence_level);
	static boost::numeric::ublas::vector<double> dResidual_i(boost::numeric::ublas::vector<double> const& values_functions1, boost::numeric::ublas::vector<double> const& values_functions2, boost::numeric::ublas::vector<double> const& prices1, boost::numeric::ublas::vector<double> const& prices2);
	static int likelihoodRatioTest(boost::numeric::ublas::vector<double> const& values_functions1, boost::numeric::ublas::vector<double> const& values_functions2, double confidence_level);
	static int likelihoodRatioTestResidual(boost::numeric::ublas::vector<double> const& vector1, boost::numeric::ublas::vector<double> const& vector2, boost::numeric::ublas::vector<double> const& prices1, boost::numeric::ublas::vector<double> const& prices2, double confidence_level);
};								  

class KernelDensity {
	//private:
public:
	static boost::numeric::ublas::vector<double> kde_multi(boost::numeric::ublas::matrix<double> x_simulated, boost::numeric::ublas::vector<double> x_realized);
	static double kde(boost::numeric::ublas::vector<double> const& x_simulated, double x_realized);
};
#endif
