#pragma once

#include <boost/numeric/ublas/matrix.hpp>
#include <Eigen/Core>
#include <Spectra/Util/SelectionRule.h>


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

    static double smallest_eigval(
		boost::numeric::ublas::matrix<double> const& input
	);

	static double eig_norm_error(
		boost::numeric::ublas::matrix<double> const& m_A,
		boost::numeric::ublas::vector<double> const& v_x,
		double const lambda
	);
	static boost::numeric::ublas::vector<double> eig_all_norm_errors(
		boost::numeric::ublas::matrix<double> const& m_A,
		boost::numeric::ublas::matrix<double> const& m_x,
		boost::numeric::ublas::vector<double> const& v_lambda
	);

	static boost::numeric::ublas::matrix<double> clean_data(
		boost::numeric::ublas::matrix<double> const& m
	);
};
