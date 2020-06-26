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
#include <boost/container/vector.hpp>
#include <boost/fusion/algorithm/transformation/push_back.hpp>
#include <boost/fusion/include/push_back.hpp>
#include <boost/range/algorithm.hpp>
#include <boost/range/algorithm_ext.hpp>
#include <boost/assign/std/vector.hpp>
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/variance.hpp>
#include <boost/bind.hpp>

#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <strstream>
#include <sstream>
#include <random>
#include <chrono>

typedef std::chrono::high_resolution_clock Clock;
using namespace boost::numeric;

ublas::vector<double> read_time_series_test(std::string const& file_name);
ublas::matrix<double> gen_start_params(size_t n, std::string dist);
ublas::vector<double> runBFGS_Gaussian(size_t n, ublas::vector<double> series, int max_iter, double epsilon);
ublas::vector<double> runBFGS_T(size_t n, ublas::vector<double> series, int max_iter, double epsilon);
ublas::vector<double> getUniformTimeseries(ublas::vector<double> series, ublas::vector<double> params);
ublas::vector<double> GARCH_vec(ublas::vector<double> time_series, ublas::vector<double> x);
ublas::matrix<double> gen_copula_params(size_t n, int nRiskFactors, std::string dist);
ublas::vector<double> runBFGS_TCopula(size_t n, T_Copula* t, int max_iter, double epsilon);
ublas::vector<double> runBFGS_normCopula(size_t n, Gaussian_Copula* norm, int max_iter, double epsilon);

int main() {
    ublas::vector<double> time_series;
    ublas::matrix<double> hist_rf;
    ublas::matrix<double> E_rf;

    try {
        hist_rf = read_csv_matrix("fHist.csv"); // Senaste kurvan sist
    } catch (std::exception & ex) {
        std::cout << "Error:" << ex.what() << std::endl;
        return 1;
    }

    //L�s in egenmatrisen f�r riskfria kurvan
    try {
        E_rf = read_csv_matrix("EZero.csv"); // Senaste kurvan sist
    } catch (std::exception & ex) {
        std::cout << "Error:" << ex.what() << std::endl;
        return 1;
    }
    auto t1 = Clock::now();

    // Ber�kna riskfaktorer p� riskfria kurvan.
    ublas::matrix<double> delta_f(matrixOperations::diff_matrix(hist_rf));
    ublas::matrix<double> hist_riskfactors = prod(trans(E_rf), trans(delta_f));
    hist_riskfactors = trans(hist_riskfactors);

    //Save uniform variables in U matrix for copula estimation
    ublas::matrix<double> U(hist_riskfactors.size1(), hist_riskfactors.size2());
    size_t nRiskfactors{ U.size2() };

    //Save optimal parameters for all riskfactors
    boost::container::vector<ublas::vector<double>> OptParamsAll(0);

    // Optimize parameters for all riskfactors
    int nSolutions{ 2 };
    int max_iter{ 100 };
    double epsilon{ pow(10, -7) };

    for (size_t i{ 0 }; i < nRiskfactors; ++i) {
        ublas::matrix_column<ublas::matrix<double>> xi(hist_riskfactors, i); // Get risk factor nr i.

        ublas::vector<double> optGaussian_xi = runBFGS_Gaussian(nSolutions, xi, max_iter, epsilon);
        ublas::vector<double> optStudent_xi = runBFGS_T(nSolutions, xi, max_iter, epsilon);

        //Get uniform Timeseries
        if (optGaussian_xi(4) < optStudent_xi(5)) {    // Choose Gaussian marginal dist
            ublas::matrix_column<ublas::matrix<double>> U_column(U, i);
            U_column = getUniformTimeseries(xi, optGaussian_xi);
            OptParamsAll.push_back(optGaussian_xi);
        }
        else {                                        // Choose Student t marginals
            ublas::matrix_column<ublas::matrix<double>> U_column(U, i);
            U_column = getUniformTimeseries(xi, optStudent_xi);
            OptParamsAll.push_back(optStudent_xi);
        }
    }

    //Copula estimation
    T_Copula dist(U);
    T_Copula* TC = &dist;

    Gaussian_Copula distG(U);
    Gaussian_Copula* gaussianC = &distG;

    ublas::vector<double> t_copula_results = runBFGS_TCopula(10, TC, max_iter, epsilon);
    ublas::vector<double> norm_copula_results = runBFGS_normCopula(10, gaussianC, max_iter, epsilon);

    //ublas::matrix<double> P_t = TC->buildP(t_copula_results);

    std::cout << "OptParams first eigenvector = " << OptParamsAll[0] << std::endl;
    std::cout << "OptParams second eigenvector = " << OptParamsAll[1] << std::endl;
    std::cout << "OptParams third eigenvector = " << OptParamsAll[2] << std::endl;

    //std::cout << "P Students t = " << P_t << std::endl;
    std::cout << "FV Students t copula = " << t_copula_results(t_copula_results.size() - 1) << std::endl;

    //ublas::matrix<double> P_norm = gaussianC->buildP(norm_copula_results);
    //std::cout << "P Gaussian = " << P_norm << std::endl;
    std::cout << "FV Gaussian copula = " << norm_copula_results(norm_copula_results.size() - 1) << std::endl;
    
    auto t2 = Clock::now();
    std::cout << "Delta t2-t1: "
        << std::chrono::duration_cast<std::chrono::seconds>(t2 - t1).count()
        << " seconds" << std::endl;
}

ublas::vector<double> getUniformTimeseries(ublas::vector<double> series, ublas::vector<double> params) {
    ublas::vector<double> garchVec(GARCH_vec(series, params));
    double mu{ params(3) };
    ublas::vector<double> U(series.size());

    if (params.size() == 5) {
        for (size_t i{ 0 }; i < series.size(); ++i) {
            boost::math::normal norm = boost::math::normal::normal_distribution(0, 1);
            U(i) = cdf(norm, (series(i) - mu) / sqrt(garchVec(i)));
            if (U(i) == 1) {
                U(i) = U(i - 1);
            }
        }
    }
    else {
        boost::math::students_t stud = boost::math::students_t::students_t_distribution(params(4));
        for (size_t i{ 0 }; i < series.size(); ++i) {
            U(i) = cdf(stud, (series(i) - mu) / sqrt(garchVec(i)));
            if (U(i) == 1) {
                U(i) = U(i - 1);
            }
        }
    }

    return U;
}

ublas::vector<double> GARCH_vec(ublas::vector<double> time_series, ublas::vector<double> x) {
    ublas::vector<double> garch_vec(time_series.size() + 1);

    //Calc variance for timeseries
    boost::accumulators::accumulator_set<double, boost::accumulators::stats<boost::accumulators::tag::variance>> acc;
    for_each(time_series.begin(), time_series.end(), boost::bind<void>(boost::ref(acc), _1));

    //Set variance as first element
    garch_vec(0) = boost::accumulators::variance(acc);

    for (size_t i{ 0 }; i < time_series.size(); i++) {
        garch_vec(i + 1) = x(0) + x(1) * pow(time_series(i), 2) + x(2) * garch_vec(i);
    }

    return garch_vec;
}

ublas::matrix<double> gen_copula_params(size_t n, int nRiskFactors, std::string dist) {
    size_t nParams{ 0 };

    if (dist == "normal") {
        nParams = nRiskFactors;
    } else if (dist == "t") {
        nParams = nRiskFactors + 1;
    }

    ublas::vector<double> params(nParams);
    ublas::matrix<double> param_matrix(nParams, n);

    std::random_device rd;
    std::default_random_engine generator(rd());
    std::uniform_real_distribution<double> distribution(-0.6, 0.6);

    for (size_t i{ 0 }; i < n; ++i) {
        for (size_t i{ 0 }; i < nParams - 1; ++i) {
            params(i) = distribution(generator);
        }

        if (dist == "t") {
            params(nParams - 1) = 5;
        }
        else {
            params(nParams - 1) = distribution(generator);
        }

        column(param_matrix, i) = params;
    }

    return param_matrix;
}

ublas::vector<double> runBFGS_TCopula(size_t n, T_Copula* distribution, int max_iter, double epsilon) {
    size_t nRiskFactors{ (distribution->time_series).size2() };

    //Create Student_t object and set Hessian
    //T_Copula dist(series);
    //T_Copula* distribution = &dist;
    ublas::matrix<double> H_inv(nRiskFactors + 1, nRiskFactors + 1);
    ublas::identity_matrix<double> I(nRiskFactors + 1);
    H_inv = I;
    ublas::matrix<double> params(nRiskFactors + 1, n);
    ublas::vector<double> FV(n);

    params = gen_copula_params(n, nRiskFactors, "t");

    // Run optimization problem n times.
    for (size_t i{ 0 }; i < params.size2(); ++i) {
        ublas::vector<double> results = bfgs::minimize(column(params, i), H_inv, max_iter, epsilon, distribution);

        for (size_t j{ 0 }; j < params.size1(); ++j) {
            column(params, i)(j) = results(j);
        }

        FV(i) = results(params.size1());
    }

    //Get smallest loglikelihood value and corresponding parameters
    double smallest{ FV(0) };
    int index{ 0 };
    for (ublas::vector<double>::iterator it = FV.begin(); it != FV.end(); it++) {
        if (*it < smallest) {
            smallest = *it;
            index = it - FV.begin();
        }
    }

    std::cout << "All FV:S  : " << FV << std::endl;
    std::cout << "OPT FV = " << smallest << " at params = " << ublas::column(params, index) << std::endl;

    //Save optimal parameters and function value in opt_params
    ublas::vector<double> opt_params(nRiskFactors + 1);

    for (size_t i{ 0 }; i < opt_params.size() - 1; ++i) {
        opt_params(i) = params(i, index);
    }
    opt_params(opt_params.size() - 1) = smallest;

    return opt_params;
}

ublas::vector<double> runBFGS_normCopula(size_t n, Gaussian_Copula* distribution, int max_iter, double epsilon) {
    size_t nRiskFactors{ (distribution->time_series).size2() };

    //Create Student_t object and set Hessian
    //Gaussian_Copula dist(series);
    //Gaussian_Copula* distribution = &dist;
    ublas::matrix<double> H_inv(nRiskFactors, nRiskFactors);
    ublas::identity_matrix<double> I(nRiskFactors);
    H_inv = I;
    ublas::matrix<double> params(nRiskFactors, n);
    ublas::vector<double> FV(n);

    params = gen_copula_params(n, nRiskFactors, "normal");
    // Run optimization problem n times.
    for (size_t i{ 0 }; i < params.size2(); ++i) {
        ublas::vector<double> results = bfgs::minimize(column(params, i), H_inv, max_iter, epsilon, distribution);

        for (size_t j{ 0 }; j < params.size1(); ++j) {
            column(params, i)(j) = results(j);
        }
        FV(i) = results(params.size1());
    }

    //Get smallest loglikelihood value and corresponding parameters
    double smallest = FV(0);
    int index{ 0 };
    for (ublas::vector<double>::iterator it = FV.begin(); it != FV.end(); it++) {
        if (*it < smallest) {
            smallest = *it;
            index = it - FV.begin();
        }
    }

    std::cout << "All FV:S  : " << FV << std::endl;
    std::cout << "OPT FV = " << smallest << " at params = " << column(params, index) << std::endl;

    //Save optimal parameters and function value in opt_params
    ublas::vector<double> opt_params(nRiskFactors);

    for (size_t i{ 0 }; i < opt_params.size() - 1; ++i) {
        opt_params(i) = params(i, index);
    }
    opt_params(opt_params.size() - 1) = smallest;

    return opt_params;
}

ublas::vector<double> runBFGS_Gaussian(size_t n, ublas::vector<double> series, int max_iter, double epsilon) {
    //Convert vector as matrix for BFGS minimize function
    ublas::matrix<double> timeseries(series.size(), 1);
    ublas::matrix_column<ublas::matrix<double>> mc(timeseries, 0);
    mc = series;

    //Create Student_t object and set Hessian
    Gaussian dist(timeseries);
    Gaussian* distribution = &dist;
    ublas::matrix<double> H_inv(4, 4);
    ublas::identity_matrix<double> I(4);
    H_inv = I;
    ublas::matrix<double> params(4, n);
    ublas::vector<double> FV(n);

    params = gen_start_params(n, "normal");
    // Run optimization problem n times.
    for (size_t i{ 0 }; i < params.size2(); ++i) {
        ublas::vector<double> results = bfgs::minimize(column(params, i), H_inv, max_iter, epsilon, distribution);

        for (size_t j{ 0 }; j < params.size1(); ++j) {
            column(params, i)(j) = results(j);
        }
        FV(i) = results(params.size1());
    }

    //Get smallest loglikelihood value and corresponding parameters
    double smallest{ FV(0) };
    int index{ 0 };
    for (ublas::vector<double>::iterator it = FV.begin(); it != FV.end(); it++) {
        if (*it < smallest) {
            smallest = *it;
            index = it - FV.begin();
        }
    }

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

ublas::vector<double> runBFGS_T(size_t n, ublas::vector<double> series, int max_iter, double epsilon) {
    //Convert vector as matrix for BFGS minimize function
    ublas::matrix<double> timeseries(series.size(), 1);
    ublas::matrix_column<ublas::matrix<double>> mc(timeseries, 0);
    mc = series;

    //Create Student_t object and set Hessian
    Student_t dist_t(timeseries);
    Student_t* distribution_t = &dist_t;
    ublas::matrix<double> H_inv(5, 5);
    ublas::identity_matrix<double> I(5);
    H_inv = I;
    ublas::matrix<double> params(5, n);
    ublas::vector<double> FV(n);

    params = gen_start_params(n, "t");
    // Run optimization problem n times.
    for (size_t i{ 0 }; i < params.size2(); ++i) {
        ublas::vector<double> results = bfgs::minimize(column(params, i), H_inv, max_iter, epsilon, distribution_t);

        for (size_t j{ 0 }; j < params.size1(); ++j) {
            column(params, i)(j) = results(j);
        }
        FV(i) = results(params.size1());
    }

    //Get smallest loglikelihood value and corresponding parameters
    double smallest{ FV(0) };
    int index{ 0 };

    for (ublas::vector<double>::iterator it = FV.begin(); it != FV.end(); it++) {
        if (*it < smallest) {
            smallest = *it;
            index = it - FV.begin();
        }
    }

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

ublas::matrix<double> gen_start_params(size_t n, std::string dist) {
    int nParams{ 4 };
    if (dist == "normal") {
        nParams = 4;
    }
    else if (dist == "t") {
        nParams = 5;
    }
    ublas::vector<double> params(nParams);
    ublas::matrix<double> param_matrix(nParams, n);

    std::default_random_engine generator;
    std::uniform_real_distribution<double> distribution(0.7, 0.95);

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

ublas::vector<double> read_time_series_test(std::string const& file_name) {
    std::fstream fin;
    std::string line;
    std::string value_string;
    double value{ 0.0 };

    fin.open(file_name, std::ios::in);

    if (!fin) {
        throw std::runtime_error("Could not open file");
    }

    ublas::vector<double> series(7981);
    int k{ 0 };

    while (getline(fin, line)) {
        std::stringstream s(line);

        //Get value
        getline(s, value_string, '\n');
        const char* decimal = value_string.c_str();  // String to const char
        value = std::strtod(decimal, NULL); // const char to float.

        //series(k) = value;
        series(series.size() - 1 - k) = value;
        k = k + 1;
    }
    fin.close();

    ublas::vector<double> log_returns(series.size() - 1);
    for (size_t i{ 0 }; i < series.size() - 1; ++i) {
        log_returns(i) = log(series(i + 1) / series(i));
    }

    return log_returns;
}
