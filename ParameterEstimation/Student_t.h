#pragma once
#define _USE_MATH_DEFINES
#include "Distribution.h"

#include <boost/numeric/ublas/matrix.hpp>

class Student_t : public Distribution {
public:
    Student_t(boost::numeric::ublas::matrix<double> series);

    void update_GARCH_vec(boost::numeric::ublas::vector<double> const& x);
    double function_value(boost::numeric::ublas::vector<double> const& x);
    boost::numeric::ublas::vector<double> calcGradients(boost::numeric::ublas::vector<double> const& x);
    boost::numeric::ublas::vector<double> calcNumGradients(boost::numeric::ublas::vector<double> const& x);
    double calcStepSize(
        boost::numeric::ublas::vector<double> const& x, 
        boost::numeric::ublas::vector<double> const& d
    );

private:
    boost::numeric::ublas::vector<double> m_GARCH_vec;
    boost::numeric::ublas::vector<double> time_series;

    boost::numeric::ublas::vector<double> derivative_w(boost::numeric::ublas::vector<double> const& x);
    boost::numeric::ublas::vector<double> derivative_a(boost::numeric::ublas::vector<double> const& x);
    boost::numeric::ublas::vector<double> derivative_b(boost::numeric::ublas::vector<double> const& x);
};
