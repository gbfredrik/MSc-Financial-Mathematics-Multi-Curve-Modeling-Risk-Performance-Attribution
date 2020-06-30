#include "../../CurveLibrary/sample_handler.h"
#include "../../MathLibrary/matrixOperations.h"
#include "../../ParameterEstimation/bfgs.h"
#include "../../ParameterEstimation/Distribution.h"
#include "../../ParameterEstimation/T_Copula.h"
#include "../../ParameterEstimation/Gaussian_Copula.h"


#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/io.hpp>

#include <iostream>
#include <vector>

#include <chrono>

typedef std::chrono::high_resolution_clock Clock;
using namespace boost::numeric;

//ublas::matrix<double> gen_start_params(size_t const n, std::string dist);
//ublas::vector<double> run_bfgs_gaussian(size_t const n, ublas::vector<double> const& series, int const max_iter, double const epsilon);
//ublas::vector<double> run_bfgs_t(size_t const n, ublas::vector<double> const& series, int const max_iter, double const epsilon);
//ublas::vector<double> get_uniform_timeseries(ublas::vector<double> const& series, ublas::vector<double> const& params);
//ublas::vector<double> garch_vec(ublas::vector<double> const& time_series, ublas::vector<double> const& x);
//ublas::matrix<double> gen_copula_params(size_t const n, int const nRiskFactors, std::string dist);
//ublas::vector<double> run_bfgs_t_copula(size_t const n, T_Copula* t, int const max_iter, double const epsilon);
//ublas::vector<double> run_bfgs_gaussian_copula(size_t const n, Gaussian_Copula* norm, int const max_iter, double const epsilon);

int main() {
    ublas::vector<double> time_series;
    ublas::matrix<double> hist_rf;
    ublas::matrix<double> E_rf;

    try {
        hist_rf = read_csv_matrix("fHist.csv"); // Senaste kurvan sist
    } catch (std::exception& ex) {
        std::cout << "Error:" << ex.what() << std::endl;
        return 1;
    }

    //Läs in egenmatrisen för riskfria kurvan
    try {
        E_rf = read_csv_matrix("EZero.csv"); // Senaste kurvan sist
    } catch (std::exception& ex) {
        std::cout << "Error:" << ex.what() << std::endl;
        return 1;
    }

    auto t1{ Clock::now() };
    // Beräkna riskfaktorer på riskfria kurvan.
    ublas::matrix<double> delta_f{ matrixOperations::diff_matrix(hist_rf) };
    ublas::matrix<double> hist_riskfactors{ prod(trans(E_rf), trans(delta_f)) };
    hist_riskfactors = trans(hist_riskfactors);

    //Save uniform variables in U matrix for copula estimation
    ublas::matrix<double> U(hist_riskfactors.size1(), hist_riskfactors.size2());
    size_t nRiskfactors{ U.size2() };

    //Save optimal parameters for all riskfactors
    std::vector<ublas::vector<double>> OptParamsAll(0);

    // Optimize parameters for all riskfactors
    int nSolutions{ 100 }; // 2
    int max_iter{ 100 };
    double epsilon{ pow(10, -7) };

    for (size_t i{ 0 }; i < nRiskfactors; ++i) {
        ublas::matrix_column<ublas::matrix<double>> xi{ hist_riskfactors, i }; // Get risk factor nr i.

        ublas::vector<double> optGaussian_xi{ bfgs::run_bfgs_gaussian(nSolutions, xi, max_iter, epsilon) };
        ublas::vector<double> optStudent_xi{ bfgs::run_bfgs_t(nSolutions, xi, max_iter, epsilon) };

        //Get uniform Timeseries
        if (optGaussian_xi(4) < optStudent_xi(5)) {    // Choose Gaussian marginal dist
            ublas::matrix_column<ublas::matrix<double>> U_column(U, i);
            U_column = bfgs::get_uniform_timeseries(xi, optGaussian_xi);
            OptParamsAll.push_back(optGaussian_xi);
        } else {                                        // Choose Student t marginals
            ublas::matrix_column<ublas::matrix<double>> U_column(U, i);
            U_column = bfgs::get_uniform_timeseries(xi, optStudent_xi);
            OptParamsAll.push_back(optStudent_xi);
        }
    }

    //Copula estimation
    T_Copula dist(U);
    T_Copula* TC{ &dist };

    Gaussian_Copula distG(U);
    Gaussian_Copula* gaussianC{ &distG };

    ublas::vector<double> t_copula_results{ bfgs::run_bfgs_t_copula(10, TC, max_iter, epsilon) };
    ublas::vector<double> norm_copula_results{ bfgs::run_bfgs_gaussian_copula(10, gaussianC, max_iter, epsilon) };

    std::cout << "OptParams first eigenvector = " << OptParamsAll[0] << std::endl;
    std::cout << "OptParams second eigenvector = " << OptParamsAll[1] << std::endl;
    std::cout << "OptParams third eigenvector = " << OptParamsAll[2] << std::endl;
    std::cout << "FV Students t copula = " << t_copula_results(t_copula_results.size() - 1) << std::endl;
    std::cout << "FV Gaussian copula = " << norm_copula_results(norm_copula_results.size() - 1) << std::endl;
    
    auto t2{ Clock::now() };
    std::cout << "Delta t2-t1: "
        << std::chrono::duration_cast<std::chrono::seconds>(t2 - t1).count()
        << " seconds" << std::endl;
}
