#pragma once

#include <boost/numeric/ublas/matrix.hpp>

class Backtesting {
public:
    // Value-at-risk
    static bool VaR_hypothesis_test(
        boost::numeric::ublas::vector<double> const& VaRs,
        boost::numeric::ublas::vector<double> const& PnLs,
        double const c,
        double const alpha
    );
    static bool VaR_Christoffersen_test(
        boost::numeric::ublas::vector<double> const& VaRs,
        boost::numeric::ublas::vector<double> const& PnLs,
        double const alpha
    );

    // Expected shortfall
    static bool ES_Acerbi_Szekely(
        boost::numeric::ublas::matrix<double> const& X,
        boost::numeric::ublas::vector<double> const& VaRs,
        boost::numeric::ublas::vector<double> const& ESs,
        boost::numeric::ublas::vector<double> const& PnLs,
        double const phi
    );

private:
    // VaR hypothesis backtesting - Helper functions
    static int test_statistic_X_T(
        boost::numeric::ublas::vector<double> const& VaRs,
        boost::numeric::ublas::vector<double> const& PnLs
    );
    static double test_statistic_Z(int const X_T, int const T, double const p);

    // VaR Christoffersen's test - Helper functions
    static void number_of_periods(
        boost::numeric::ublas::vector<double> const& VaRs,
        boost::numeric::ublas::vector<double> const& PnLs,
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

    static double test_statistic_Z1(
        boost::numeric::ublas::vector<double> const& VaRs,
        boost::numeric::ublas::vector<double> const& ESs,
        boost::numeric::ublas::vector<double> const& PnLs
    );

    // Helper functions
    static boost::numeric::ublas::vector<double> VaR_breaches(
        boost::numeric::ublas::vector<double> const& VaRs,
        boost::numeric::ublas::vector<double> const& PnLs
    );
};