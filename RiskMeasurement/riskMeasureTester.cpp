#include "RiskMeasures.h"

#include "../MathLibrary/rvSim.h"

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>

#include <iostream>

using namespace boost::numeric;

int main() {
	ublas::vector<double> test_losses(10);
	
	test_losses = row(rvSim::gen_normal(1, 10), 0);

	//std::cout << RiskMeasures::VAR(test_losses, 0.95);
	
}