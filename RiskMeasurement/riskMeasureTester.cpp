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
    int window{ 1000 };

    ublas::vector<double> returns( read_csv_vector("returns.csv") );

    //std::cout << test_losses << "\n";

    std::cout << RiskMeasures::VaR(returns, c) << std::endl;
    //std::cout << RiskMeasures::ES(returns, c) << std::endl;

    ublas::vector<double> VaRs(RiskMeasures::VaR_series(returns, c, window));
    ublas::vector<double> ESs(RiskMeasures::ES_series(returns, c, window));
    //std::cout << VaRs << std::endl;
    ublas::vector_range<ublas::vector<double>> PnLs(
        returns, 
        ublas::range(window, returns.size())
    );

    std::cout << "Returns: " << returns.size() << std::endl;
    std::cout << "VaRs: " << VaRs.size() << std::endl;
    std::cout << "ESs: " << ESs.size() << std::endl;
    std::cout << "PnLs: " << PnLs.size() << std::endl << std::endl;
    
    std::cout << "VaR Backtest:" << std::endl;
    std::cout << "Reject H_0 for H_1: " << RiskMeasures::VaR_hypothesis_test(
        VaRs,
        PnLs,
        c,
        0.05
    ) << std::endl << std::endl;
    
    std::cout << "Christoffersen: " << RiskMeasures::VaR_Christoffersen_test(
        VaRs,
        ublas::vector_range<ublas::vector<double>>(
            returns,
            ublas::range(window, returns.size())
            ),
        0.05
    ) << std::endl;

    std::cout << "ES Backtest:" << std::endl;
    std::cout << "Acerbi & Szekely, Reject H_0 for H_1: " << RiskMeasures::ES_Acerbi_Szekely(VaRs, ESs, PnLs, 0.05) << std::endl;

    write_csv_vector(VaRs, "VaRs.csv");
    write_csv_vector(ESs, "ESs.csv");
    write_csv_vector(PnLs, "PnLs.csv");
}

 //std::cout << RiskMeasures::test_statistic_christoffersen(
    //    0.037260825780463,
    //    0.036610878661088,
    //    0.054054054054054,
    //    921,
    //    35,
    //    35,
    //    2
    //) << std::endl; // Correct!