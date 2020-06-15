#pragma once

#define _USE_MATH_DEFINES
#include <cmath>

#include <boost/numeric/ublas/matrix.hpp>

class RiskMeasures {
public:
    // Value-at-risk
    static double VaR(boost::numeric::ublas::vector<double> const& outcomes, double const c);
    //static boost::numeric::ublas::vector<double> VaR_series(
    //    boost::numeric::ublas::vector<double>& outcomes, 
    //    double const c, 
    //    int const window
    //);

    // Expected shortfall
    static double ES(boost::numeric::ublas::vector<double> const& outcomes, double const c);
    //static boost::numeric::ublas::vector<double> ES_series(
    //    boost::numeric::ublas::vector<double>& outcomes, 
    //    double const c,
    //    int const window
    //);

    // Kernel density estimator
    static double kde(
        boost::numeric::ublas::vector<double> const& x_simulated,
        double const x_realized
    );
    static boost::numeric::ublas::vector<double> kde_multi(
        boost::numeric::ublas::matrix<double> const& x_simulated,
        boost::numeric::ublas::vector<double> const& x_realized
    );

private:
    static int VaR_index(double const c, size_t const n);
};
