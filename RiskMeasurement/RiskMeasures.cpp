#include "RiskMeasures.h"

#include "../MathLibrary/matrixOperations.h"
#include "../MathLibrary/statisticsOperations.h"

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/math/distributions/inverse_chi_squared.hpp>

using namespace boost::numeric;

// --- VaR ---
double RiskMeasures::VaR(
	ublas::vector<double> const& outcomes,
	double const c
) {
	ublas::vector<double> sorted_losses(-outcomes);
	std::sort(sorted_losses.begin(), sorted_losses.end());

	int index_VaR{ VaR_index(c, sorted_losses.size()) };
	//std::cout << index_VaR << "\n";

	return sorted_losses(index_VaR);
}

ublas::vector<double> RiskMeasures::VaR_series(
	ublas::vector<double>& outcomes, 
	double const c,
	int const window
) {
	size_t n{ outcomes.size() };
	ublas::vector<double> VaR_measures(n - window);

	for (size_t i{ 0 }; i < n - window; ++i) {
		VaR_measures(i) = VaR(
			ublas::vector_range<ublas::vector<double>> (
				outcomes, 
				ublas::range(i, i + window)
			), 
			c
		);
	}
	
	return VaR_measures;
}

// --- VaR Backtesting---
// VaR hypothesis backtesting
bool RiskMeasures::VaR_hypothesis_test(
	ublas::vector<double> const& VaR, 
	ublas::vector<double> const& PnL, 
	double const c,
	double const alpha
) {
	int X_T{ test_statistic_X_T(VaR, PnL) };
	int T{ static_cast<int>(VaR.size()) };
	double p{ 1 - c };
	double Z{ test_statistic_Z(X_T, T, p) };
	double tilde_m{ statisticsOperations::invCDFNorm(1 - alpha) }; // = invCDFNorm( 1 - alpha, 0.0, 1.0);
	
	std::cout << "Z: " << Z << std::endl;
	std::cout << "tilde_m: " << tilde_m << std::endl;
	return Z > tilde_m;
}

int RiskMeasures::test_statistic_X_T(ublas::vector<double> const& VaR, ublas::vector<double> const& PnL) {
	int X_T{ 0 };
	for (size_t i{ 0 }, n{ VaR.size() }; i < n; ++i) {
		if (-VaR(i) > PnL(i)) {
			++X_T;
		}
	}
	
	return X_T;
}

double RiskMeasures::test_statistic_Z(int const X_T, int const T, double const p) {
	return (X_T - T * p) / sqrt(T * p * (1.0 - p));
}

// VaR Christoffersen's test


// --- ES ---
double RiskMeasures::ES(
	ublas::vector<double> const& outcomes,
	double const c
) {
	ublas::vector<double> sorted_losses(-outcomes);
	std::sort(sorted_losses.begin(), sorted_losses.end());

	int index_VaR{ VaR_index(c, sorted_losses.size()) };

	ublas::vector_range<ublas::vector<double>> exceedances(
		sorted_losses,
		ublas::range(index_VaR, sorted_losses.size())
	);

	return matrixOperations::vector_average(exceedances);
}

// --- ES Backtesting ---
// Test 1: Acerbi and Szekely


// Kernel Density Estimator

// --- Helper functions ---
int RiskMeasures::VaR_index(double const c, size_t const n) {
	return static_cast<int>(c * n);
}