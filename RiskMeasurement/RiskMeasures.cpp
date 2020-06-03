#include "RiskMeasures.h"

#include <boost/numeric/ublas/matrix.hpp>

using namespace boost::numeric;

double RiskMeasures::VAR(
	ublas::vector<double> const& losses, 
	double const conf_level
) {
	ublas::vector<double> sorted_losses(losses);
	std::sort(sorted_losses.begin(), sorted_losses.end());

	int index_VAR{ static_cast<int>(ceil(conf_level * sorted_losses.size())) };
	std::cout << index_VAR << "\n";

	return sorted_losses(index_VAR);
}

double RiskMeasures::ES(
	boost::numeric::ublas::vector<double> const& losses,
	double const conf_level
) {


	return 0.0;
}
