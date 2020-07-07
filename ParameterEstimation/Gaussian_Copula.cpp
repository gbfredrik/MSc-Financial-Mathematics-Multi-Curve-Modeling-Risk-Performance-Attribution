#include "Gaussian_Copula.h"

#include "../MathLibrary/matrixOperations.h"
#include "../MathLibrary/statisticsOperations.h"
#include "../RiskFactorCalculation/FactorCalculation.h"

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/math/special_functions/digamma.hpp>

#include <cmath>

using namespace boost::numeric;

Gaussian_Copula::Gaussian_Copula(ublas::matrix<double> series) 
    : Distribution(series), time_series{ series }
{
}

double Gaussian_Copula::function_value(ublas::vector<double> const& x) {
    double sum{ 0.0 };
    size_t n{ time_series.size2() }; //Number of riskfaktors
    ublas::matrix<double> P{ build_p(x) };
    double det_P{ matrixOperations::matrix_det(P) };
    ublas::matrix<double> P_inv{ matrixOperations::matrix_inv(P) };

    for (size_t i{ 0 }, rows{ time_series.size1() }; i < rows; ++i) {
        ublas::matrix_row<ublas::matrix<double>> U_row(time_series, i);
        ublas::vector<double> N_inv(n);

        for (size_t j{ 0 }, k{ N_inv.size() }; j < k; ++j) {
            N_inv(j) = statisticsOperations::invCDFNorm(U_row(j));
        }

        //Help products and summations
        ublas::identity_matrix<double> I(n);
        ublas::vector<double> nPI{ prod(N_inv, P_inv - I) };
        double nPIn{ inner_prod(nPI, N_inv) };
        
        sum += - 0.5 * nPIn - 0.5 * log(det_P);
    }

    return -sum;
}

ublas::vector<double> Gaussian_Copula::calc_gradients(ublas::vector<double> const& x) {
    ublas::vector<double> gradients(x.size());
    size_t n{ time_series.size2() }; //Number of riskfaktors
    ublas::zero_vector<double> zeroVec(n * n);
    ublas::vector<double> dFdP{ zeroVec };

    //Get rho gradients as a vector
    for (size_t i{ 0 }, rows{ time_series.size1() }; i < rows; ++i) {
        //Get T_inv(U)
        ublas::matrix_row<ublas::matrix<double>> U_row(time_series, i);
        ublas::vector<double> N_inv(n);

        for (size_t j{ 0 }, k{ N_inv.size() }; j < k; ++j) {
            N_inv(j) = statisticsOperations::invCDFNorm(U_row(j));
        }

        ublas::matrix<double> P{ build_p(x) };
        ublas::matrix<double> P_inv{ matrixOperations::matrix_inv(P) };
        ublas::vector<double> vec_Pinv{ matrixOperations::matrixToVector(P_inv) };

        //Kroneckers product of two vectors
        ublas::vector<double> temp{ prod(N_inv, P_inv) };
        ublas::vector<double> vKron{ matrixOperations::kron_prod_vec(temp, temp) };

        //Get gradients for time_series(i)
        ublas::vector<double> dFdP_temp{ 0.5 * vKron - 0.5 * vec_Pinv };
        dFdP += dFdP_temp;
    }

    //Get rho gradients as matrix
    ublas::matrix<double> dFdP_mat{ matrixOperations::vectorToMatrix(dFdP, time_series.size2()) };

    //Get optimization parameters from rho matrix
    ublas::vector<double> dfdParams{ get_elements(dFdP_mat) };

    return -dfdParams;
}

ublas::matrix<double> Gaussian_Copula::build_p(ublas::vector<double> const& x) {
    size_t n{ time_series.size2() };
    ublas::matrix<double> P(n, n);
    size_t counter{ 0 };    //keep track of fetched elements

    for (size_t i{ 0 }; i < n; ++i) {
        for (size_t j{ 0 }; j < i + 1; ++j) {
            if (i == j) {
                P(i, j) = 1;
            } else {
                P(i, j) = x(counter);
                P(j, i) = x(counter);
                ++counter;
            }
        }
    }

    return P;
}

ublas::vector<double> Gaussian_Copula::get_elements(ublas::matrix<double> const& P) {
    size_t n{ time_series.size2() };
    ublas::vector<double> resVec(static_cast<size_t>((n - 1) * n * 0.5)); // Vector with optimization parameters in P
    size_t counter{ 0 };    //keep track of fetched elements

    for (size_t i{ 0 }; i < n; ++i) {
        for (size_t j{ 0 }; j < i; ++j) {
            resVec(counter) = P(i, j);
            ++counter;
        }
    }

    return resVec;
}

ublas::vector<double> Gaussian_Copula::calc_num_gradients(ublas::vector<double> const& x) {
    double epsilon{ 2.2 * pow(10, -16) };
    ublas::vector<double> increment{ sqrt(epsilon) * x };
    ublas::vector<double> num_gradients(x.size());

    //ublas::vector<double> x_0diff(x);
    //x_0diff(0) += increment(0);
    //ublas::vector<double> x_1diff(x);
    //x_1diff(1) += increment(1);
    //ublas::vector<double> x_2diff(x);
    //x_2diff(2) += increment(2);
    //ublas::vector<double> x_3diff(x);
    //x_3diff(3) += increment(3);
    //ublas::vector<double> x_4diff(x);
    //x_4diff(4) += increment(4);
    //
    //num_gradients(0) = (function_value(x_0diff) - function_value(x)) / (x_0diff(0) - x(0));
    //num_gradients(1) = (function_value(x_1diff) - function_value(x)) / (x_1diff(1) - x(1));
    //num_gradients(2) = (function_value(x_2diff) - function_value(x)) / (x_2diff(2) - x(2));
    //num_gradients(3) = (function_value(x_3diff) - function_value(x)) / (x_3diff(3) - x(3));
    //num_gradients(4) = (function_value(x_4diff) - function_value(x)) / (x_4diff(4) - x(4));

    return num_gradients;
}

double Gaussian_Copula::calc_step_size(
    ublas::vector<double> const& x, 
    ublas::vector<double> const& d
) {
    double a{ 1.0 };
    bool accepted{ false };
    //Kontrollera att rho är positiv definit, dvs minsta egenvärdet är positivt, annars halvera steglängden.

    while (!accepted) {
        accepted = true;
        for (size_t i{ 0 }; i < x.size(); ++i) { // Todo: Check condition x.size()
            if (x(i) + a * d(i) < -1 || x(i) + a * d(i) > 1) {
                accepted = false;
                break;
            }
        }

        ublas::matrix<double> Pnext{ build_p(x + a * d) };
        double minEigenvalue{ FactorCalculation::smallest_eigval(Pnext) };

        if (minEigenvalue <= 0) {
            accepted = false;
        }

        if (!accepted) {
            a *= 0.5;
        }

        if (a == 0.0) {
            break;
        }
    }

    return a;
}
