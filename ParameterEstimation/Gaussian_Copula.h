#pragma once
#include "Distribution.h"

#include <boost/numeric/ublas/matrix.hpp>

class Gaussian_Copula : public Distribution {
public:
    boost::numeric::ublas::matrix<double> time_series;

    Gaussian_Copula(boost::numeric::ublas::matrix<double> series);

    double function_value(boost::numeric::ublas::vector<double> const& x);
    boost::numeric::ublas::vector<double> calcGradients(boost::numeric::ublas::vector<double> const& x);
    boost::numeric::ublas::vector<double> calcNumGradients(boost::numeric::ublas::vector<double> const& x);
    double calcStepSize(boost::numeric::ublas::vector<double> const& x, boost::numeric::ublas::vector<double> const& d);

private:
    boost::numeric::ublas::matrix<double> buildP(boost::numeric::ublas::vector<double> const& x);
    boost::numeric::ublas::vector<double> getElements(boost::numeric::ublas::matrix<double> const& matrix);
};
