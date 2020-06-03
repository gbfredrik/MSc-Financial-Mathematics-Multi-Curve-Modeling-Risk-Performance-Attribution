#pragma once

#include <boost/numeric/ublas/matrix.hpp>

class RiskMeasures {
public:
	static double VAR(boost::numeric::ublas::vector<double> const& losses, double const conf_level);
	static double ES(boost::numeric::ublas::vector<double> const& losses, double const conf_level);


	//static bool VAR_hypothesis_test();
	//static bool VAR_Christoffersen_test();

	//static bool ES_Acerbi_test();
	// Kernel density estimator - include
};