#pragma once
#define _USE_MATH_DEFINES
#include "Distribution.h"

#include <boost/numeric/ublas/matrix.hpp>

class Gaussian : public Distribution {
public: 
    Gaussian(boost::numeric::ublas::matrix<double> series);

    void update_GARCH_vec(boost::numeric::ublas::vector<double> const& x);
    double function_value(boost::numeric::ublas::vector<double> const& x);
    boost::numeric::ublas::vector<double> calcGradients(boost::numeric::ublas::vector<double> const& x);
    boost::numeric::ublas::vector<double> calcNumGradients(boost::numeric::ublas::vector<double> const& x);
    boost::numeric::ublas::matrix<double> calcNumHessian(boost::numeric::ublas::vector<double> const& x);
    double calcStepSize(
        boost::numeric::ublas::vector<double> const& x, 
        boost::numeric::ublas::vector<double> const& d
    );

private:
    boost::numeric::ublas::vector<double> time_series;
    boost::numeric::ublas::vector<double> m_GARCH_vec;
    double garch0;

    boost::numeric::ublas::vector<double> derivative_w(boost::numeric::ublas::vector<double> const& x);
    boost::numeric::ublas::vector<double> derivative_a(boost::numeric::ublas::vector<double> const& x);
    boost::numeric::ublas::vector<double> derivative_b(boost::numeric::ublas::vector<double> const& x);
};
