#include "RiskMeasures.h"

#include "../MathLibrary/matrixOperations.h"

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>

using namespace boost::numeric;

double RiskMeasures::VAR(
	ublas::vector<double> const& outcomes,
	double const conf_level
) {
	ublas::vector<double> sorted_losses(-outcomes);
	std::sort(sorted_losses.begin(), sorted_losses.end());

	int index_VAR{ VAR_index(conf_level, sorted_losses.size()) };
	std::cout << index_VAR << "\n";

	return sorted_losses(index_VAR);
}

double RiskMeasures::ES(
	boost::numeric::ublas::vector<double> const& outcomes,
	double const conf_level
) {
	ublas::vector<double> sorted_losses(-outcomes);
	std::sort(sorted_losses.begin(), sorted_losses.end());

	int index_VAR{ VAR_index(conf_level, sorted_losses.size()) };

	ublas::vector_range<ublas::vector<double>> exceedances(
		sorted_losses,
		ublas::range(index_VAR, sorted_losses.size())
	);

	return matrixOperations::vector_average(exceedances);
}

int RiskMeasures::VAR_index(double const conf_level, size_t const n) {
	return static_cast<int>(conf_level * n);
}
