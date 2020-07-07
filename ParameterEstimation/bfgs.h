#pragma once

#include "Distribution.h"
#include "T_Copula.h"
#include "Gaussian_Copula.h"

#include <boost/numeric/ublas/matrix.hpp>

//Rosenbrock function: f(x,y) = (1-x)^2 + 100(y-x^2)^2 

class bfgs {
public:
    static boost::numeric::ublas::vector<double> minimize(
        boost::numeric::ublas::vector<double> x,
        boost::numeric::ublas::matrix<double> H_inv,
        int const max_iter,
        double const epsilon,
        Distribution* dist
    );

    static boost::numeric::ublas::vector<double> get_uniform_timeseries(boost::numeric::ublas::vector<double> const& series, boost::numeric::ublas::vector<double> const& params);

    static boost::numeric::ublas::vector<double> run_bfgs_gaussian(size_t const n, boost::numeric::ublas::vector<double> const& series, int const max_iter, double const epsilon);
    static boost::numeric::ublas::vector<double> run_bfgs_t(size_t const n, boost::numeric::ublas::vector<double> const& series, int const max_iter, double const epsilon);
    
    static boost::numeric::ublas::vector<double> run_bfgs_gaussian_copula(size_t const n, Gaussian_Copula* norm, int const max_iter, double const epsilon);
    static boost::numeric::ublas::vector<double> run_bfgs_t_copula(size_t const n, T_Copula* t, int const max_iter, double const epsilon);

private:
    static boost::numeric::ublas::vector<double> garch_vec(boost::numeric::ublas::vector<double> const& time_series, boost::numeric::ublas::vector<double> const& x);
    static boost::numeric::ublas::matrix<double> gen_start_params(size_t const n, std::string dist);
    static boost::numeric::ublas::matrix<double> gen_copula_params(size_t const n, int const nRiskFactors, std::string dist);
};
