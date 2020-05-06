#pragma once

#include <boost/numeric/ublas/matrix.hpp>
#include <Eigen/Core>


class FactorCalculation {
private:

public:
	static bool iram(boost::numeric::ublas::matrix<double> const& input, int k, boost::numeric::ublas::matrix<double>& m_E, boost::numeric::ublas::vector<double>& v_Lambda);
	static bool eigen_bdcsvd(boost::numeric::ublas::matrix<double> const& input, boost::numeric::ublas::matrix<double>& m_E, boost::numeric::ublas::vector<double>& v_Lambda);


};
