#include "Gaussian_Copula.h"

#include "../MathLibrary/matrixOperations.h"
#include "../MathLibrary/statisticsOperations.h"
#include "../RiskFactorCalculation/FactorCalculation.h"

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/math/statistics/bivariate_statistics.hpp>
#include <boost/math/constants/constants.hpp>
#include <boost/qvm/mat_operations.hpp>
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/variance.hpp>
#include <boost/bind.hpp>
#include <boost/math/special_functions/digamma.hpp>

#include <Eigen/Dense>

#include <iostream>
#include <numeric>
#include <cmath>

using namespace boost::numeric;

Gaussian_Copula::Gaussian_Copula(ublas::matrix<double> series) : Distribution(series) {
    time_series = series;
}

void Gaussian_Copula::getSeries() {
    //std::cout << time_series << "\n";
}

double Gaussian_Copula::function_value(ublas::vector<double> const& x) {
    double sum{ 0.0 };
    size_t n{ time_series.size2() }; //Number of riskfaktors
    ublas::matrix<double> P{ buildP(x) };
    double det_P{ matrixOperations::ublasToMatrixXd(P).determinant() };
    ublas::matrix<double> P_inv{ matrixOperations::matrixXdToUblas(matrixOperations::ublasToMatrixXd(P).inverse()) };

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

ublas::vector<double> Gaussian_Copula::calcGradients(ublas::vector<double> const& x) {
    ublas::vector<double> gradients(x.size());
    size_t n{ time_series.size2() }; //Number of riskfaktors
    double sum{ 0.0 };
    ublas::zero_vector<double> zeroVec(n * n);
    ublas::vector<double> dFdP{ zeroVec };
    double dFdnu{ 0.0 };

    //Get rho gradients as a vector
    for (size_t i{ 0 }, rows{ time_series.size1() }; i < rows; ++i) {
        //Get T_inv(U)
        ublas::matrix_row<ublas::matrix<double>> U_row(time_series, i);
        ublas::vector<double> N_inv(n);

        for (size_t j{ 0 }, k{ N_inv.size() }; j < k; ++j) {
            N_inv(j) = statisticsOperations::invCDFNorm(U_row(j));
        }

        ublas::matrix<double> P{ buildP(x) };
        ublas::matrix<double> P_inv{ matrixOperations::matrixXdToUblas(matrixOperations::ublasToMatrixXd(P).inverse()) };

        //Get P_inv as a vector of the columns
        ublas::vector<double> vec_Pinv{ matrixToVector(P_inv) };

        //Kroneckers product of two vectors
        ublas::vector<double> vKron{ kronOfVectors(prod(N_inv, P_inv), prod(N_inv, P_inv)) };

        //Get gradients for time_series(i)
        ublas::vector<double> dFdP_temp{ 0.5 * vKron - 0.5 * vec_Pinv };
        dFdP += dFdP_temp;
    }

    //Get rho gradients as matrix
    ublas::matrix<double> dFdP_mat{ vectorToMatrix(dFdP) };

    //Get optimization parameters from rho matrix
    ublas::vector<double> dfdParams{ getElements(dFdP_mat) };

    return -dfdParams;
}

ublas::matrix<double> Gaussian_Copula::vectorToMatrix(ublas::vector<double> const& vec) {
    size_t n{ time_series.size2() };
    ublas::matrix<double> resMatrix(n, n);
    size_t counter{ 0 };

    for (size_t j{ 0 }; j < n; ++j) {
        for (size_t i{ 0 }; i < n; ++i) {
            resMatrix(i, j) = vec(counter);
            ++counter;
        }
    }

    return resMatrix;
}

ublas::matrix<double> Gaussian_Copula::buildP(ublas::vector<double> const& x) {
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


ublas::vector<double> Gaussian_Copula::getElements(ublas::matrix<double> const& P) {
    size_t n{ time_series.size2() };
    ublas::vector<double> resVec((n - 1) * n * 0.5); // Vector with optimization parameters in P
    size_t counter{ 0 };    //keep track of fetched elements

    for (size_t i{ 0 }; i < n; ++i) {
        for (size_t j{ 0 }; j < i; ++j) {
            resVec(counter) = P(i, j);
            ++counter;
        }
    }

    return resVec;
}

ublas::vector<double> Gaussian_Copula::kronOfVectors(ublas::vector<double> const& v1, ublas::vector<double> const& v2) {
    ublas::vector<double> vKron(v1.size() * v2.size());
    size_t counter{ 0 };

    for (size_t i{ 0 }; i < v1.size(); ++i) {
        for (size_t j{ 0 }; j < v2.size(); ++j) {
            vKron(counter) = v1(i) * v2(j);
            ++counter;
        }
    }

    return vKron;
}

ublas::vector<double> Gaussian_Copula::matrixToVector(ublas::matrix<double> const& matrix) {
    size_t n{ matrix.size1() };
    ublas::vector<double> resVec(n * n);
    size_t counter{ 0 };

    for (size_t j{ 0 }; j < n; ++j) {
        for (size_t i{ 0 }; i < n; ++i) {
            resVec(counter) = matrix(i, j);
            ++counter;
        }
    }

    return resVec;
}

ublas::vector<double> Gaussian_Copula::calcNumGradients(ublas::vector<double> const& x) {
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

double Gaussian_Copula::calcStepSize(ublas::vector<double> const& x, ublas::vector<double> const& d) {
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

        ublas::matrix<double> Pnext{ buildP(x + a * d) };
        double minEigenvalue{ FactorCalculation::smallest_eigval(Pnext) };

        if (minEigenvalue <= 0) {
            accepted = false;
        }

        if (!accepted) {
            a *= 0.5;
        }

        if (a == 0) {
            break;
        }
    }

    return a;
}
