#include "RiskMeasures.h"

#include "../MathLibrary/matrixOperations.h"
#include "../MathLibrary/statisticsOperations.h"

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/operation.hpp>
#include <boost/math/distributions/inverse_chi_squared.hpp>

using namespace boost::numeric;

// --- VaR ---
double RiskMeasures::VaR(
    ublas::vector<double> const& outcomes,
    double const c
) {
    ublas::vector<double> sorted_outcomes(-outcomes);
    std::sort(sorted_outcomes.begin(), sorted_outcomes.end());

    int index_VaR{ VaR_index(c, sorted_outcomes.size()) };

    return sorted_outcomes(index_VaR);
}

//ublas::vector<double> RiskMeasures::VaR_series(
//    ublas::vector<double>& outcomes,
//    double const c,
//    int const window
//) {
//    size_t n{ outcomes.size() };
//    ublas::vector<double> VaR_measures(n - window);
//
//    for (size_t i{ 0 }; i < n - window; ++i) {
//        VaR_measures(i) = VaR(
//            ublas::vector_range<ublas::vector<double>> (
//                outcomes,
//                ublas::range(i, i + window)
//            ), 
//            c
//        );
//    }
//    
//    return VaR_measures;
//}

// --- VaR Backtesting---
// VaR hypothesis backtesting
bool RiskMeasures::VaR_hypothesis_test(
    ublas::vector<double> const& VaRs, 
    ublas::vector<double> const& PnLs, 
    double const c,
    double const alpha
) {
    int X_T{ test_statistic_X_T(VaRs, PnLs) };
    int T{ static_cast<int>(VaRs.size()) };
    double p{ 1 - c };
    double Z{ test_statistic_Z(X_T, T, p) };
    double tilde_m{ statisticsOperations::invCDFNorm(1 - alpha) }; // = invCDFNorm( 1 - alpha, 0.0, 1.0);
    
    return Z > tilde_m;
}

int RiskMeasures::test_statistic_X_T(
    ublas::vector<double> const& VaRs, 
    ublas::vector<double> const& PnLs
) {
    return sum(VaR_breaches(VaRs, PnLs));
}

double RiskMeasures::test_statistic_Z(int const X_T, int const T, double const p) {
    return (X_T - T * p) / sqrt(T * p * (1.0 - p));
}

// VaR Christoffersen's test
bool RiskMeasures::VaR_Christoffersen_test(
    ublas::vector<double> const& VaRs, 
    ublas::vector<double> const& PnLs, 
    double const alpha
) {
    double u_00{ 0.0 };
    double u_01{ 0.0 };
    double u_10{ 0.0 };
    double u_11{ 0.0 };
    number_of_periods(VaRs, PnLs, u_00, u_01, u_10, u_11);
    
    double pi{ (u_01 + u_11) / (u_00 + u_01 + u_10 + u_11) };
    double pi_01{ (u_01) / (u_00 + u_01) };
    double pi_11{ (u_11) / (u_10 + u_11) };

    double christoffersen_statistic{ 
        test_statistic_christoffersen(
            pi,
            pi_01,
            pi_11,
            u_00,
            u_01,
            u_10,
            u_11
        )
    };

    double z{ statisticsOperations::invCDFchi2(1 - alpha) };

    return christoffersen_statistic > z;
}

void RiskMeasures::number_of_periods(
    ublas::vector<double> const& VaRs,
    ublas::vector<double> const& PnLs,
    double& u_00,
    double& u_01,
    double& u_10,
    double& u_11
) {
    ublas::vector<int> indicator(VaR_breaches(VaRs, PnLs));
    
    for (size_t i{ 0 }, n{ VaRs.size() - 1}; i < n; ++i) {
        bool exceed_current{ indicator(i) == 1};
        bool exceed_next{ indicator(i + 1) == 1};

        if (!exceed_current && !exceed_next) {
            ++u_00;
        } else if (!exceed_current && exceed_next) {
            ++u_01;
        } else if (exceed_current && !exceed_next) {
            ++u_10;
        } else if (exceed_current && exceed_next) {
            ++u_11;
        }
    }
}

double RiskMeasures::test_statistic_christoffersen(
    double const pi, 
    double const pi_01, 
    double const pi_11, 
    double const u_00,
    double const u_01,
    double const u_10,
    double const u_11
) {
    return -2.0 * ((u_00 + u_10) * log(1.0 - pi) + (u_01 + u_11) * log(pi))
        + 2.0 * (u_00 * log(1.0 - pi_01) + u_01 * log(pi_01) + u_10 * log(1.0 - pi_11) + u_11 * log(pi_11));
}

// --- ES ---
double RiskMeasures::ES(
    ublas::vector<double> const& outcomes,
    double const c
) {
    ublas::vector<double> sorted_outcomes(-outcomes);
    std::sort(sorted_outcomes.begin(), sorted_outcomes.end());

    int index_VaR{ VaR_index(c, sorted_outcomes.size()) };

    ublas::vector_range<ublas::vector<double>> exceedances(
        sorted_outcomes,
        ublas::range(index_VaR, sorted_outcomes.size())
    );

    return matrixOperations::vector_average(exceedances);
}

//ublas::vector<double> RiskMeasures::ES_series(
//    ublas::vector<double>& PnLs,
//    double const c,
//    int const window
//) {
//    size_t n{ PnLs.size() };
//    ublas::vector<double> ES_measures(n - window);
//
//    for (size_t i{ 0 }; i < n - window; ++i) {
//        ES_measures(i) = ES(
//            ublas::vector_range<ublas::vector<double>>(
//                PnLs,
//                ublas::range(i, i + window)
//            ),
//            c
//        );
//    }
//
//    return ES_measures;
//}

// --- ES Backtesting ---
// Test 1: Acerbi and Szekely

bool RiskMeasures::ES_Acerbi_Szekely(
    ublas::matrix<double> const& X,
    ublas::vector<double> const& VaRs, 
    ublas::vector<double> const& ESs, 
    ublas::vector<double> const& PnLs, 
    double const phi
) {
    size_t M{ X.size2() }; // columns
    ublas::vector<double> Z(M);
    double Z_X_actual{ test_statistic_Z1(VaRs, ESs, PnLs) };
    int n{ 0 };

    for (size_t i{ 0 }; i < M; ++i) {
        Z(i) = test_statistic_Z1(VaRs, ESs, column(X, i)); // Todo: Control
        if (Z(i) < Z_X_actual) {
            ++n;
        }
    }

    double p{ n / static_cast<double>(M) };

    return p < phi;
}

double RiskMeasures::test_statistic_Z1(
    ublas::vector<double> const& VaRs, 
    ublas::vector<double> const& ESs, 
    ublas::vector<double> const& PnLs
) {
    ublas::vector<int> indicator(VaR_breaches(VaRs, PnLs));

    int N_T{ sum(indicator) };
    if (N_T == 0) { // Exit if no breaches to avoid div. by zero in Z_X.
        return 0;
    }

    double Z_X{ sum(element_div(element_prod(PnLs, indicator), ESs)) / N_T + 1.0 };
    
    return Z_X;
}

// Kernel Density Estimator
// Todo: Include

// --- Helper functions ---
int RiskMeasures::VaR_index(double const c, size_t const n) {
    return static_cast<int>(c * n);
}

ublas::vector<double> RiskMeasures::VaR_breaches(
    ublas::vector<double> const& VaRs, 
    ublas::vector<double> const& PnLs
) {
    size_t T{ VaRs.size() };
    ublas::vector<double> indicator(T);

    for (size_t i{ 0 }; i < T; ++i) {
        indicator(i) = -VaRs(i) > PnLs(i);
    }
    
    return indicator;
}
