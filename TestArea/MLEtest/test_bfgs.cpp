#include "../../CurveLibrary/sample_handler.h"
#include "../../MathLibrary/matrixOperations.h"
#include "../../ParameterEstimation/bfgs.h"
#include "../../ParameterEstimation/Distribution.h"
#include "../../ParameterEstimation/Gaussian.h"
#include "../../ParameterEstimation/Student_t.h"
#include "../../ParameterEstimation/T_Copula.h"
#include "../../ParameterEstimation/Gaussian_Copula.h"

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/math/distributions/normal.hpp>
#include <boost/math/distributions/students_t.hpp>
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/variance.hpp>
#include <boost/bind.hpp>

#include <iostream>
#include <vector>
#include <algorithm>
#include <string>
#include <random>
#include <chrono>

typedef std::chrono::high_resolution_clock Clock;
using namespace boost::numeric;

ublas::matrix<double> gen_start_params(size_t const n, std::string dist);
ublas::vector<double> runBFGS_Gaussian(size_t const n, ublas::vector<double> const& series, int const max_iter, double const epsilon);
ublas::vector<double> runBFGS_T(size_t const n, ublas::vector<double> const& series, int const max_iter, double const epsilon);
ublas::vector<double> getUniformTimeseries(ublas::vector<double> const& series, ublas::vector<double> const& params);
ublas::vector<double> GARCH_vec(ublas::vector<double> const& time_series, ublas::vector<double> const& x);
ublas::matrix<double> gen_copula_params(size_t const n, int const nRiskFactors, std::string dist);
ublas::vector<double> runBFGS_TCopula(size_t const n, T_Copula* t, int const max_iter, double const epsilon);
ublas::vector<double> runBFGS_normCopula(size_t const n, Gaussian_Copula* norm, int const max_iter, double const epsilon);

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
    int nSolutions{ 2 };
    int max_iter{ 100 };
    double epsilon{ pow(10, -7) };

    for (size_t i{ 0 }; i < nRiskfactors; ++i) {
        ublas::matrix_column<ublas::matrix<double>> xi{ hist_riskfactors, i }; // Get risk factor nr i.

        ublas::vector<double> optGaussian_xi{ runBFGS_Gaussian(nSolutions, xi, max_iter, epsilon) };
        ublas::vector<double> optStudent_xi{ runBFGS_T(nSolutions, xi, max_iter, epsilon) };

        //Get uniform Timeseries
        if (optGaussian_xi(4) < optStudent_xi(5)) {    // Choose Gaussian marginal dist
            ublas::matrix_column<ublas::matrix<double>> U_column(U, i);
            U_column = getUniformTimeseries(xi, optGaussian_xi);
            OptParamsAll.push_back(optGaussian_xi);
        } else {                                        // Choose Student t marginals
            ublas::matrix_column<ublas::matrix<double>> U_column(U, i);
            U_column = getUniformTimeseries(xi, optStudent_xi);
            OptParamsAll.push_back(optStudent_xi);
        }
    }

    //Copula estimation
    T_Copula dist(U);
    T_Copula* TC{ &dist };

    Gaussian_Copula distG(U);
    Gaussian_Copula* gaussianC{ &distG };

    ublas::vector<double> t_copula_results{ runBFGS_TCopula(10, TC, max_iter, epsilon) };
    ublas::vector<double> norm_copula_results{ runBFGS_normCopula(10, gaussianC, max_iter, epsilon) };

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

ublas::vector<double> getUniformTimeseries(ublas::vector<double> const& series, ublas::vector<double> const& params) {
    ublas::vector<double> garchVec{ GARCH_vec(series, params) };
    double mu{ params(3) };
    ublas::vector<double> U(series.size());

    if (params.size() == 5) {
        for (size_t i{ 0 }, n{ series.size() }; i < n; ++i) {
            boost::math::normal norm{ boost::math::normal::normal_distribution(0, 1) };
            U(i) = cdf(norm, (series(i) - mu) / sqrt(garchVec(i)));
            if (U(i) == 1) {
                U(i) = U(i - 1);
            }
        }
    } else {
        boost::math::students_t stud{ boost::math::students_t::students_t_distribution(params(4)) };
        for (size_t i{ 0 }, n{ series.size() }; i < n; ++i) {
            U(i) = cdf(stud, (series(i) - mu) / sqrt(garchVec(i)));
            if (U(i) == 1) {
                U(i) = U(i - 1);
            }
        }
    }

    return U;
}

ublas::vector<double> GARCH_vec(ublas::vector<double> const& time_series, ublas::vector<double> const& x) {
    ublas::vector<double> garch_vec(time_series.size() + 1);

    //Calc variance for timeseries
    boost::accumulators::accumulator_set<double, boost::accumulators::stats<boost::accumulators::tag::variance>> acc;
    for_each(time_series.begin(), time_series.end(), boost::bind<void>(boost::ref(acc), _1));

    //Set variance as first element
    garch_vec(0) = boost::accumulators::variance(acc);

    for (size_t i{ 0 }, n{ time_series.size() }; i < n; ++i) {
        garch_vec(i + 1) = x(0) + x(1) * pow(time_series(i), 2) + x(2) * garch_vec(i);
    }

    return garch_vec;
}

ublas::matrix<double> gen_copula_params(size_t const n, int const nRiskFactors, std::string dist) {
    size_t nParams{ 0 };

    if (dist == "normal") {
        nParams = nRiskFactors;
    } else if (dist == "t") {
        nParams = nRiskFactors + 1;
    }

    ublas::vector<double> params(nParams);
    ublas::matrix<double> param_matrix(nParams, n);

    std::random_device rd;
    std::default_random_engine generator{ rd() };
    std::uniform_real_distribution<double> distribution{ -0.6, 0.6 };

    for (size_t i{ 0 }; i < n; ++i) {
        for (size_t j{ 0 }; j < nParams - 1; ++j) {
            params(j) = distribution(generator);
        }

        if (dist == "t") {
            params(nParams - 1) = 5;
        } else {
            params(nParams - 1) = distribution(generator);
        }

        column(param_matrix, i) = params;
    }

    return param_matrix;
}

ublas::vector<double> runBFGS_TCopula(size_t const n, T_Copula* distribution, int const max_iter, double const epsilon) {
    size_t nRiskFactors{ (distribution->time_series).size2() };
    ublas::identity_matrix<double> I(nRiskFactors + 1);
    ublas::matrix<double> H_inv{ I };
    ublas::matrix<double> params(nRiskFactors + 1, n);
    ublas::vector<double> FV(n);

    params = gen_copula_params(n, nRiskFactors, "t");

    // Run optimization problem n times.
    for (size_t i{ 0 }, cols{ params.size2() }; i < cols; ++i) {
        ublas::vector<double> results{ bfgs::minimize(column(params, i), H_inv, max_iter, epsilon, distribution) };

        for (size_t j{ 0 }, rows{ params.size1() }; j < rows; ++j) {
            column(params, i)(j) = results(j);
        }

        FV(i) = results(params.size1());
    }

    //Get smallest loglikelihood value and corresponding parameters
    int index{ std::min_element(FV.begin(), FV.end()) - FV.begin() };
    double smallest{ FV(index) };

    std::cout << "All FV:S  : " << FV << std::endl;
    std::cout << "OPT FV = " << smallest << " at params = " << column(params, index) << std::endl;

    //Save optimal parameters and function value in opt_params
    ublas::vector<double> opt_params(nRiskFactors + 1);

    for (size_t i{ 0 }, len{ opt_params.size() - 1 }; i < len; ++i) {
        opt_params(i) = params(i, index);
    }
    opt_params(opt_params.size() - 1) = smallest;

    return opt_params;
}

ublas::vector<double> runBFGS_normCopula(size_t const n, Gaussian_Copula* distribution, int const max_iter, double const epsilon) {
    size_t nRiskFactors{ (distribution->time_series).size2() };
    ublas::identity_matrix<double> I(nRiskFactors);
    ublas::matrix<double> H_inv{ I };
    ublas::matrix<double> params(nRiskFactors, n);
    ublas::vector<double> FV(n);

    params = gen_copula_params(n, nRiskFactors, "normal");
    // Run optimization problem n times.
    for (size_t i{ 0 }, cols{ params.size2() }; i < cols; ++i) {
        ublas::vector<double> results{ bfgs::minimize(column(params, i), H_inv, max_iter, epsilon, distribution) };

        for (size_t j{ 0 }, rows{ params.size1() }; j < rows; ++j) {
            column(params, i)(j) = results(j);
        }

        FV(i) = results(params.size1());
    }

    //Get smallest loglikelihood value and corresponding parameters
    int index{ std::min_element(FV.begin(), FV.end()) - FV.begin() };
    double smallest{ FV(index) };

    std::cout << "All FV:S  : " << FV << std::endl;
    std::cout << "OPT FV = " << smallest << " at params = " << column(params, index) << std::endl;

    //Save optimal parameters and function value in opt_params
    ublas::vector<double> opt_params(nRiskFactors);

    for (size_t i{ 0 }, len{ opt_params.size() - 1 }; i < len; ++i) {
        opt_params(i) = params(i, index);
    }
    opt_params(opt_params.size() - 1) = smallest;

    return opt_params;
}

ublas::vector<double> runBFGS_Gaussian(size_t const n, ublas::vector<double> const& series, int const max_iter, double const epsilon) {
    //Convert vector as matrix for BFGS minimize function
    ublas::matrix<double> timeseries(series.size(), 1);
    ublas::matrix_column<ublas::matrix<double>> mc(timeseries, 0);
    mc = series;

    //Create Student_t object and set Hessian
    Gaussian dist(timeseries);
    Gaussian* distribution = &dist;
    ublas::identity_matrix<double> I(4);
    ublas::matrix<double> H_inv{ I };
    ublas::matrix<double> params{ gen_start_params(n, "normal") };
    ublas::vector<double> FV(n);

    //params = gen_start_params(n, "normal");
    // Run optimization problem n times.
    for (size_t i{ 0 }, cols{ params.size2() }; i < cols; ++i) {
        ublas::vector<double> results{ bfgs::minimize(column(params, i), H_inv, max_iter, epsilon, distribution) };

        for (size_t j{ 0 }, rows{ params.size1() }; j < rows; ++j) {
            column(params, i)(j) = results(j);
        }
        FV(i) = results(params.size1());
    }

    //Get smallest loglikelihood value and corresponding parameters
    int index{ std::min_element(FV.begin(), FV.end()) - FV.begin() };
    double smallest{ FV(index) };

    std::cout << "All FV:S  : " << FV << std::endl;
    std::cout << "OPT FV = " << smallest << " at params = " << column(params, index) << std::endl;

    //Save optimal parameters and function value in opt_params
    ublas::vector<double> opt_params(5);
    opt_params(0) = params(0, index);
    opt_params(1) = params(1, index);
    opt_params(2) = params(2, index);
    opt_params(3) = params(3, index);
    opt_params(4) = smallest;

    return opt_params;
}

ublas::vector<double> runBFGS_T(size_t const n, ublas::vector<double> const& series, int const max_iter, double const epsilon) {
    //Convert vector as matrix for BFGS minimize function
    ublas::matrix<double> timeseries(series.size(), 1);
    ublas::matrix_column<ublas::matrix<double>> mc{ timeseries, 0 };
    mc = series;

    //Create Student_t object and set Hessian
    Student_t dist_t{ timeseries };
    Student_t* distribution_t{ &dist_t };
    ublas::identity_matrix<double> I(5);
    ublas::matrix<double> H_inv{ I };
    ublas::matrix<double> params{ gen_start_params(n, "t") };
    ublas::vector<double> FV(n);

    // Run optimization problem n times.
    for (size_t i{ 0 }, cols{ params.size2() }; i < cols; ++i) {
        ublas::vector<double> results{ bfgs::minimize(column(params, i), H_inv, max_iter, epsilon, distribution_t) };

        for (size_t j{ 0 }, rows{ params.size1() }; j < rows; ++j) {
            column(params, i)(j) = results(j);
        }
        FV(i) = results(params.size1());
    }

    //Get smallest loglikelihood value and corresponding parameters
    int index{ std::min_element(FV.begin(), FV.end()) - FV.begin() };
    double smallest{ FV(index) };

    std::cout << "All FV:S Student t: " << FV << std::endl;
    std::cout << "OPT FV Student T= " << smallest << " at params = " << column(params, index) << std::endl;

    //Save optimal parameters and function value in opt_params
    ublas::vector<double> opt_params(6);
    opt_params(0) = params(0, index);
    opt_params(1) = params(1, index);
    opt_params(2) = params(2, index);
    opt_params(3) = params(3, index);
    opt_params(4) = params(4, index);
    opt_params(5) = smallest;

    return opt_params;
}

ublas::matrix<double> gen_start_params(size_t const n, std::string dist) { // Todo: Split into gen_start_t and gen_start_normal
    int nParams{ 0 };
    if (dist == "normal") {
        nParams = 4;
    } else if (dist == "t") {
        nParams = 5;
    }

    ublas::vector<double> params(nParams);
    ublas::matrix<double> param_matrix(nParams, n);

    std::default_random_engine generator;
    std::uniform_real_distribution<double> distribution{ 0.7, 0.95 };

    for (size_t i{ 0 }; i < n; ++i) {
        params(0) = 0.0001;
        params(2) = distribution(generator);
        params(1) = 1 - params(2) - 0.001;
        params(3) = 0.0001;

        if (dist == "t") {
            params(4) = 5;
        }

        column(param_matrix, i) = params;
    }

    return param_matrix;
}
