#pragma once

#include <boost/numeric/ublas/matrix.hpp>

class RiskMeasures {
public:
	static double VaR(boost::numeric::ublas::vector<double> const& outcomes, double const c);
	static boost::numeric::ublas::vector<double> VaR_series(
		boost::numeric::ublas::vector<double>& outcomes, 
		double const c, 
		int const window
	);

	static bool VaR_hypothesis_test(
		boost::numeric::ublas::vector<double> const& VaR,
		boost::numeric::ublas::vector<double> const& PnL,
		double const c,
		double const alpha
	);
	//static bool VaR_Christoffersen_test();

	static double ES(boost::numeric::ublas::vector<double> const& outcomes, double const c);

	//static bool ES_Acerbi_test();
	// Kernel density estimator - include

private:
	static int VaR_index(double const c, size_t const n);

	// VaR hypothesis backtesting - Helper functions
	static int test_statistic_X_T(
		boost::numeric::ublas::vector<double> const& VaR,
		boost::numeric::ublas::vector<double> const& PnL
	);
	static double test_statistic_Z(int const X_T, int const T, double const p);

	// VaR Christoffersen's test - Helper functions


};