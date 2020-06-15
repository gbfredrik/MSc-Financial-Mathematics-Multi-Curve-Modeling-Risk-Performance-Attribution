#include "RiskMeasures.h"

#include "../MathLibrary/matrixOperations.h"
#include "../MathLibrary/statisticsOperations.h"

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>

using namespace boost::numeric;

// Value-at-risk
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

// Expected shortfall
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

// Kernel density estimator
double RiskMeasures::kde(
    ublas::vector<double> const& x_simulated,
    double const x_realized
) {
    size_t n{ x_simulated.size() };
    double variance{ matrixOperations::vector_variance(x_simulated) };
    double sigma{ sqrt(variance) };
    double h{ pow(4.0 / (3.0 * n), 1.0 / 5.0) * sigma };
    double sum{ 0 };

    for (size_t i{ 0 }; i < n; ++i) {
        sum += std::exp(-pow(x_realized - x_simulated(i), 2) / (2 * pow(h, 2)));
    }

    return (1.0 / (sqrt(2.0 * M_PI) * h * n)) * sum;
}

ublas::vector<double> RiskMeasures::kde_multi(
    ublas::matrix<double> const& x_simulated,
    ublas::vector<double> const& x_realized
) {
    size_t m{ x_realized.size() };
    ublas::vector<double> f(m);

    for (size_t j{ 0 }; j < m; ++j) {
        f(j) = kde(column(x_simulated, j), x_realized(j));
    }

    return f;
}

// Helper functions
int RiskMeasures::VaR_index(double const c, size_t const n) {
    return static_cast<int>(c * n);
}
