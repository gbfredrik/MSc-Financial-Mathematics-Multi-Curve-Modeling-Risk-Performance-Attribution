#pragma once

#include <boost/numeric/ublas/matrix.hpp>

class RiskMeasures {
public:
    static double VaR(boost::numeric::ublas::vector<double> const& outcomes, double const c);
    static boost::numeric::ublas::vector<double> VaR_series(
        boost::numeric::ublas::vector<double>& outcomes, 
        double const c, 
        int const window
    );

    static bool VaR_hypothesis_test(
        boost::numeric::ublas::vector<double> const& VaR,
        boost::numeric::ublas::vector<double> const& PnL,
        double const c,
        double const alpha
    );
    static bool VaR_Christoffersen_test(
        boost::numeric::ublas::vector<double> const& VaR,
        boost::numeric::ublas::vector<double> const& PnL,
        double const alpha
    );

    static double ES(boost::numeric::ublas::vector<double> const& outcomes, double const c);
    static boost::numeric::ublas::vector<double> ES_series(
        boost::numeric::ublas::vector<double>& outcomes, 
        double const c,
        int const window
    );

    static bool ES_Acerbi_Szekely(
        boost::numeric::ublas::vector<double> const& VaR,
        boost::numeric::ublas::vector<double> const& ES,
        boost::numeric::ublas::vector<double> const& PnL,
        double const c
    );
    // Kernel density estimator - include

//private:
    static int VaR_index(double const c, size_t const n);

    // VaR hypothesis backtesting - Helper functions
    static int test_statistic_X_T(
        boost::numeric::ublas::vector<double> const& VaR,
        boost::numeric::ublas::vector<double> const& PnL
    );
    static double test_statistic_Z(int const X_T, int const T, double const p);

    // VaR Christoffersen's test - Helper functions
    static void number_of_periods(
        boost::numeric::ublas::vector<double> const& VaR,
        boost::numeric::ublas::vector<double> const& PnL,
        double& u_00,
        double& u_01,
        double& u_10,
        double& u_11
    );
    static double test_statistic_christoffersen(
        double const pi,
        double const pi_01,
        double const pi_11,
        double const u_00,
        double const u_01,
        double const u_10,
        double const u_11
    );

};