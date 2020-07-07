#pragma once
#include "Distribution.h"

#include <boost/numeric/ublas/matrix.hpp>

class T_Copula : public Distribution {
public:
    boost::numeric::ublas::matrix<double> time_series;

    T_Copula(boost::numeric::ublas::matrix<double> series);

    double function_value(boost::numeric::ublas::vector<double> const& x);
    boost::numeric::ublas::vector<double> calc_gradients(boost::numeric::ublas::vector<double> const& x);
    boost::numeric::ublas::vector<double> calc_num_gradients(boost::numeric::ublas::vector<double> const& x);
    double calc_step_size(
        boost::numeric::ublas::vector<double> const& x,
        boost::numeric::ublas::vector<double> const& d
    );

private:
    boost::numeric::ublas::matrix<double> build_p(boost::numeric::ublas::vector<double> const& x);
    boost::numeric::ublas::vector<double> get_elements(boost::numeric::ublas::matrix<double> const& matrix);

    double dGamma(double const t);
};
