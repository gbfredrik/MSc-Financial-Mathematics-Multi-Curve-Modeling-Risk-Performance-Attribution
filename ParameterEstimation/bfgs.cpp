#include "bfgs.h"

#include "Distribution.h"
#include "Gaussian.h"
#include "Student_t.h"
#include "T_Copula.h"
#include "Gaussian_Copula.h"

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/math/distributions/normal.hpp>
#include <boost/math/distributions/students_t.hpp>
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/variance.hpp>
#include <boost/bind.hpp>

#include <algorithm>
#include <string>
#include <random>

using namespace boost::numeric;

ublas::vector<double> bfgs::minimize(
    ublas::vector<double> x,
    ublas::matrix<double> H_inv,
    int const max_iter,
    double const epsilon,
    Distribution* dist
) {
    size_t n{ x.size() };
    ublas::matrix<double> help_prod(n, n);
    ublas::vector<double> d(n);
    ublas::vector<double> x_new(n);
    ublas::vector<double> y(n);
    ublas::vector<double> s(n);
    ublas::vector<double> gradient_vec(n);
    double alpha{};
    double scale_H{};
    double l{};
    int k{ 0 };
    //Init hessian inverse matrix
    ublas::identity_matrix<double> I(n);

    gradient_vec = dist->calc_gradients(x);

    while (ublas::norm_2(gradient_vec) > epsilon && k < max_iter) {
        //gradient_vec = dist->calc_gradients(x);
        d = -prod(H_inv, gradient_vec);
        alpha = dist->calc_step_size(x, d);
        s = alpha * d;
        x_new = x + s;
        //x_new = x + alpha * d;
        //std::cout << std::setprecision(16) << "Nya s   = " << s << std::endl;
        //std::cout << std::setprecision(16) << "Gamla s = " << x_new - x << std::endl << std::endl;

        ublas::vector<double> gradient_vec_new{ dist->calc_gradients(x_new) };
        y = gradient_vec_new - gradient_vec;
        //s = x_new - x; // Gammal lösning men verkar sämre
        scale_H = inner_prod(y, s);
        l = 1.0 / scale_H;

        if (scale_H == 0) {
            //std::cout << "Breaking since H has invalid values.\n\n";
            break;
        }

        help_prod = prod((I - l * outer_prod(s, y)), H_inv);
        H_inv = prod(help_prod, (I - l * outer_prod(y, s))) + l * outer_prod(s, s);
        x = x_new;
        gradient_vec = gradient_vec_new;
        ++k;
    }

    ublas::vector<double> results(n + 1);

    for (size_t i{ 0 }; i < results.size() - 1; ++i) {
        results(i) = x(i);
    }
    results(n) = dist->function_value(x);

    return results;
}

ublas::vector<double> bfgs::get_uniform_timeseries(ublas::vector<double> const& series, ublas::vector<double> const& params) {
    ublas::vector<double> garchVec{ garch_vec(series, params) };
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

ublas::vector<double> bfgs::run_bfgs_gaussian(size_t const n, ublas::vector<double> const& series, int const max_iter, double const epsilon) {
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
    size_t index{ static_cast<size_t>(std::min_element(FV.begin(), FV.end()) - FV.begin()) };
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

ublas::vector<double> bfgs::run_bfgs_t(size_t const n, ublas::vector<double> const& series, int const max_iter, double const epsilon) {
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
    size_t index{ static_cast<size_t>(std::min_element(FV.begin(), FV.end()) - FV.begin()) };
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

ublas::vector<double> bfgs::run_bfgs_gaussian_copula(size_t const n, Gaussian_Copula* distribution, int const max_iter, double const epsilon) {
    size_t nRiskFactors{ (distribution->time_series).size2() };
    ublas::identity_matrix<double> I(nRiskFactors);
    ublas::matrix<double> H_inv{ I };
    ublas::matrix<double> params(nRiskFactors, n);
    ublas::vector<double> FV(n);

    params = gen_copula_params(n, nRiskFactors, "normal"); // Todo: Move
    // Run optimization problem n times.
    for (size_t i{ 0 }, cols{ params.size2() }; i < cols; ++i) {
        ublas::vector<double> results{ bfgs::minimize(column(params, i), H_inv, max_iter, epsilon, distribution) };

        for (size_t j{ 0 }, rows{ params.size1() }; j < rows; ++j) {
            column(params, i)(j) = results(j);
        }

        FV(i) = results(params.size1());
    }

    //Get smallest loglikelihood value and corresponding parameters
    size_t index{ static_cast<size_t>(std::min_element(FV.begin(), FV.end()) - FV.begin()) };
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

ublas::vector<double> bfgs::run_bfgs_t_copula(size_t const n, T_Copula* distribution, int const max_iter, double const epsilon) {
    size_t nRiskFactors{ (distribution->time_series).size2() };
    ublas::identity_matrix<double> I(nRiskFactors + 1);
    ublas::matrix<double> H_inv{ I };
    ublas::matrix<double> params(nRiskFactors + 1, n);
    ublas::vector<double> FV(n);

    params = gen_copula_params(n, nRiskFactors, "t"); // Todo: Move

    // Run optimization problem n times.
    for (size_t i{ 0 }, cols{ params.size2() }; i < cols; ++i) {
        ublas::vector<double> results{ bfgs::minimize(column(params, i), H_inv, max_iter, epsilon, distribution) };

        for (size_t j{ 0 }, rows{ params.size1() }; j < rows; ++j) {
            column(params, i)(j) = results(j);
        }

        FV(i) = results(params.size1());
    }

    //Get smallest loglikelihood value and corresponding parameters
    size_t index{ static_cast<size_t>(std::min_element(FV.begin(), FV.end()) - FV.begin()) };
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

ublas::vector<double> bfgs::garch_vec(ublas::vector<double> const& time_series, ublas::vector<double> const& x) {
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

ublas::matrix<double> bfgs::gen_start_params(size_t const n, std::string dist) { // Todo: Split into gen_start_t and gen_start_normal
    size_t nParams{ 0 };
    if (dist == "normal") {
        nParams = 4;
    }
    else if (dist == "t") {
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

ublas::matrix<double> bfgs::gen_copula_params(size_t const n, int const nRiskFactors, std::string dist) {
    size_t nParams{ 0 };

    if (dist == "normal") {
        nParams = nRiskFactors;
    }
    else if (dist == "t") {
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
        }
        else {
            params(nParams - 1) = distribution(generator);
        }

        column(param_matrix, i) = params;
    }

    return param_matrix;
}
