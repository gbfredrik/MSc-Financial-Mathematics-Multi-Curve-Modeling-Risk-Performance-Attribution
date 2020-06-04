#include "RiskMeasures.h"

#include "../MathLibrary/rvSim.h"
#include "../CurveLibrary/sample_handler.h"

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/io.hpp>

#include <iostream>

using namespace boost::numeric;

int main() {
	//ublas::vector<double> test_losses(2000);
	//test_losses = row(rvSim::gen_normal(1, 2000), 0);
	//write_csv_vector(test_losses, "test_losses.csv");
	double c{ 0.95 };
	int window{ 500 };

	ublas::vector<double> returns( read_csv_vector("returns.csv") );

	//std::cout << test_losses << "\n";

	std::cout << RiskMeasures::VaR(returns, c) << std::endl;
	//std::cout << RiskMeasures::ES(returns, c) << std::endl;

	ublas::vector<double> VaRs(RiskMeasures::VaR_series(returns, c, window));
	//std::cout << VaRs << std::endl;
	ublas::vector_range<ublas::vector<double>> PnLs(returns, ublas::range(window, returns.size()));

	std::cout << "Returns: " << returns.size() << std::endl;
	std::cout << "VaRs: " << VaRs.size() << std::endl;
	std::cout << "PnLs: " << PnLs.size() << std::endl;
	
	std::cout << "Reject H_0 for H_1: " << RiskMeasures::VaR_hypothesis_test(
		VaRs,
		ublas::vector_range<ublas::vector<double>> (
			returns,
			ublas::range(window, returns.size())
		),
		c,
		0.05
	) << std::endl;

	write_csv_vector(VaRs, "VaRs.csv");
	write_csv_vector(PnLs, "PnLs.csv");
}
