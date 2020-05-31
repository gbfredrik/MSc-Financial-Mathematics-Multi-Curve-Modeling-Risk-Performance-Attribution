#ifndef BACKTESTING
#define BACKTESTING

#include <boost/numeric/ublas/matrix.hpp>

class Likelihood {
	//private:
public:
	static boost::numeric::ublas::matrix<double> d_i(boost::numeric::ublas::matrix<double> const& matrix1, boost::numeric::ublas::matrix<double> const& matrix2);
	static boost::numeric::ublas::vector<double> d(boost::numeric::ublas::matrix<double> const& d_i, double const N);
	static boost::numeric::ublas::vector<double> s(boost::numeric::ublas::vector<double> const& sigma, double const N);
	static boost::numeric::ublas::vector<double> sigma(boost::numeric::ublas::vector<double> const& d, boost::numeric::ublas::matrix<double> const& d_i, double const N);
	static boost::numeric::ublas::vector<double> probability(boost::numeric::ublas::vector<double> const& d, boost::numeric::ublas::vector<double> const& s, boost::numeric::ublas::vector<double> const& confidence_level);
	static boost::numeric::ublas::matrix<double> dResidual_i(boost::numeric::ublas::matrix<double> const& values_functions1, boost::numeric::ublas::matrix<double> const& values_functions2, boost::numeric::ublas::matrix<double> const& prices1, boost::numeric::ublas::matrix<double> const& prices2);
	static boost::numeric::ublas::vector<double> likelihoodRatioTest(boost::numeric::ublas::matrix<double> const& values_functions1, boost::numeric::ublas::matrix<double> const& values_functions2);
	static boost::numeric::ublas::vector<double> likelihoodRatioTestResidual(boost::numeric::ublas::matrix<double> const& matrix1, boost::numeric::ublas::matrix<double> const& matrix2, boost::numeric::ublas::matrix<double> const& prices1, boost::numeric::ublas::matrix<double> const& prices2);
};								  

#endif
