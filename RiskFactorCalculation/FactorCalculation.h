#pragma once

#include <boost/numeric/ublas/matrix.hpp>
#include <Eigen/Core>


class FactorCalculation {
public:
	static bool iram(
		boost::numeric::ublas::matrix<double> const& input, 
		int const k, 
		boost::numeric::ublas::matrix<double>& m_E, 
		boost::numeric::ublas::vector<double>& v_Lambda
	);
	static bool eigen_bdcsvd(
		boost::numeric::ublas::matrix<double> const& input, 
		int const k, 
		boost::numeric::ublas::matrix<double>& m_E, 
		boost::numeric::ublas::vector<double>& v_Lambda
	);

	static boost::numeric::ublas::matrix<double> compute_risk_factors(
		boost::numeric::ublas::matrix<double> const& m_E_k, 
		boost::numeric::ublas::matrix<double> const& m_delta_f
	);
};
