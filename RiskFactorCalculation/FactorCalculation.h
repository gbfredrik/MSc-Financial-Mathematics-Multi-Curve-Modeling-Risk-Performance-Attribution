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
		boost::numeric::ublas::vector<double>& v_Lambda,
		double& approximation_error,
		boost::numeric::ublas::vector<double>& v_norm_errors
	);
	static bool eigen_bdcsvd(
		boost::numeric::ublas::matrix<double> const& input,
		int const k,
		boost::numeric::ublas::matrix<double>& m_E,
		boost::numeric::ublas::vector<double>& v_Lambda,
		double& approximation_error,
		boost::numeric::ublas::vector<double>& v_norm_errors
	);
	static bool eigen_rsvd(
		boost::numeric::ublas::matrix<double> const& input,
		int const k,
		boost::numeric::ublas::matrix<double>& m_E,
		boost::numeric::ublas::vector<double>& v_Lambda,
		double& approximation_error,
		boost::numeric::ublas::vector<double>& v_norm_errors
	);

	static boost::numeric::ublas::matrix<double> compute_risk_factors(
		boost::numeric::ublas::matrix<double> const& m_E_k, 
		boost::numeric::ublas::matrix<double> const& m_delta_f
	);

    static double smallest_eigval(
		boost::numeric::ublas::matrix<double> const& input
	);

	static Eigen::MatrixXd svd_approximation(
		Eigen::MatrixXd const& m_U,
		Eigen::MatrixXd const& v_D,
		Eigen::MatrixXd const& m_V
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
	static double relativeFrobeniusNormError(
		Eigen::MatrixXd const& m_original,
		Eigen::MatrixXd const& m_approximation
	);

	static boost::numeric::ublas::matrix<double> clean_data(
		boost::numeric::ublas::matrix<double> const& m
	);
};
