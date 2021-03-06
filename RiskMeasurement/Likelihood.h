#pragma once

#include <boost/numeric/ublas/matrix.hpp>

class Likelihood {
public:
    static int likelihood_ratio_test(
        boost::numeric::ublas::vector<double> const& values_functions1,
        boost::numeric::ublas::vector<double> const& values_functions2,
        double const confidence_level
    );
    static int likelihood_ratio_test_residual(
        boost::numeric::ublas::vector<double> const& vector1,
        boost::numeric::ublas::vector<double> const& vector2,
        boost::numeric::ublas::vector<double> const& prices1,
        boost::numeric::ublas::vector<double> const& prices2,
        double const confidence_level
    );

private:
    static boost::numeric::ublas::vector<double> d_i(
        boost::numeric::ublas::vector<double> const& vector1,
        boost::numeric::ublas::vector<double> const& vector2
    );
    static double d(boost::numeric::ublas::vector<double> const& d_i, int const N);
    static double s(double const sigma, int const N);
    static double sigma(
        double const d,
        boost::numeric::ublas::vector<double> const& d_i,
        int const N
    );

    static boost::numeric::ublas::vector<double> dResidual_i(
        boost::numeric::ublas::vector<double> const& values_functions1,
        boost::numeric::ublas::vector<double> const& values_functions2,
        boost::numeric::ublas::vector<double> const& prices1,
        boost::numeric::ublas::vector<double> const& prices2
    );

    static int is_better(double const d, double const s, double const confidence_level);
};
