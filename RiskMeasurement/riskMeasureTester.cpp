#include "RiskMeasures.h"

#include "../MathLibrary/rvSim.h"

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/io.hpp>

#include <iostream>

using namespace boost::numeric;

int main() {
	ublas::vector<double> test_losses(2000);
	
	test_losses = row(rvSim::gen_normal(1, 2000), 0);
	std::cout << test_losses << "\n";

	std::cout << RiskMeasures::VAR(test_losses, 0.95) << "\n";
	std::cout << RiskMeasures::ES(test_losses, 0.95) << "\n";
}