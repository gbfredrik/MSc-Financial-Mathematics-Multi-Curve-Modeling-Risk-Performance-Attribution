#include "T_Copula.h"

#include "../RiskFactorCalculation/FactorCalculation.h"
#include "../MathLibrary/matrixOperations.h"
#include <boost/math/distributions/students_t.hpp>

#include <Eigen/Dense>

#include <iostream>
#include <numeric>
#include <cmath>

//Boost packages for numeric represenation
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/io.hpp>

//Boost packages for statistics
#include <boost/math/statistics/bivariate_statistics.hpp>
#include <boost/math/distributions/normal.hpp>

#include <boost/math/constants/constants.hpp>
#include <boost/qvm/mat_operations.hpp>

#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/variance.hpp>
#include <boost/bind.hpp>

#include <boost/math/special_functions/digamma.hpp>

using namespace boost::numeric;

T_Copula::T_Copula(ublas::matrix<double> series) : Distribution(series) {
	time_series = series;
}

void T_Copula::getSeries() {
	//std::cout << time_series << "\n";
}

double T_Copula::function_value(ublas::vector<double> const& x) {
    double nu{ x(x.size() - 1) };
    double sum{ 0 };
    size_t n{ time_series.size2() }; //Number of riskfaktors
    ublas::matrix<double> P{ buildP(x) };
    double det_P{ matrixOperations::ublasToMatrixXd(P).determinant() };
    ublas::matrix<double> P_inv{ matrixOperations::matrixXdToUblas(matrixOperations::ublasToMatrixXd(P).inverse()) };

	for (size_t i = 0; i < time_series.size1(); ++i) {
		ublas::matrix_row<ublas::matrix<double>> U_row(time_series, i);

		//Get N_inv(U)
		ublas::vector<double> T_inv(n);
		for (size_t j = 0; j < T_inv.size(); ++j) {
			boost::math::students_t T = boost::math::students_t::students_t_distribution(nu);
			T_inv(j) = quantile(T, U_row(j));
		}
			
		ublas::vector<double> tP = prod(T_inv, P_inv);
		double tPt = inner_prod(tP, T_inv);
	
		double inner_sum = 0;
		for (size_t k= 0; k < n; ++k) {
			inner_sum = inner_sum + log(1 + pow(T_inv(k),2)/nu);
		}

		sum = sum + (n - 1) * std::lgamma(nu * 0.5) + std::lgamma((nu + n) * 0.5) 
				- (nu + n) * 0.5 * log(1 + tPt / nu) - n * std::lgamma((nu + 1) * 0.5) 
				- 0.5 * log(det_P) + (nu + 1) * 0.5 * inner_sum;
	}

	return -sum;
}


ublas::vector<double> T_Copula::calcGradients(ublas::vector<double> const& x) {
	ublas::vector<double> gradients(x.size());
	double nu = x(x.size() - 1);
	int n = time_series.size2(); //Number of riskfaktors
	double sum = 0;
    ublas::zero_vector<double> zeroVec(n * n);
    ublas::vector<double> dFdP = zeroVec;
	double dFdnu = 0;

	//Get rho gradients as a vector
	for (size_t i = 0; i < time_series.size1(); ++i) {
		//Get T_inv(U)
		ublas::matrix_row<ublas::matrix<double>> U_row(time_series, i);

        ublas::vector<double> T_inv(n);
		for (size_t j = 0; j < T_inv.size(); ++j) {
			boost::math::students_t T = boost::math::students_t::students_t_distribution(nu);
			T_inv(j) = quantile(T, U_row(j));
		}

		ublas::matrix<double> P = buildP(x);
		ublas::matrix<double> P_inv = matrixOperations::matrixXdToUblas(matrixOperations::ublasToMatrixXd(P).inverse());

		//Get P_inv as a vector of the columns
		ublas::vector<double> vec_Pinv = matrixToVector(P_inv);

		//Kroneckers product of two vectors
		ublas::vector<double> vKron = kronOfVectors(prod(T_inv, P_inv), prod(T_inv, P_inv));

		//Help products and summations
		ublas::vector<double> tP = prod(T_inv, P_inv);
		double tPt = inner_prod(tP, T_inv);
		double inner_sum = 0;
		double inner_sum2 = 0;

		for (size_t k = 0; k < n; ++k) {
			inner_sum = inner_sum + log(1 + pow(T_inv(k), 2) / nu);
			inner_sum2 = inner_sum2 + pow(T_inv(k), 2) / pow(nu, 2) / (1 + pow(T_inv(k), 2) / nu);
		}


		//Get gradients for time_series(i)
		ublas::vector<double> dFdP_temp = 0.5 * (nu + n) / (nu + tPt) * vKron - 0.5 * matrixToVector(P_inv);
		dFdP = dFdP + dFdP_temp;

		dFdnu = dFdnu + (n - 1) / (2 * tgamma(nu * 0.5)) * dGamma(nu *0.5) + 1 / (2 * tgamma((nu + n) * 0.5)) * dGamma((nu + n) * 0.5) - 0.5 * log(1 + tPt / nu)
			+ (nu + n) * 0.5 * 1 / (1 + tPt / nu) * tPt / pow(nu,2) - n / (2 * tgamma((nu + 1) * 0.5)) * dGamma((nu + 1) * 0.5)
			+ 0.5 * inner_sum - (nu + 1) * 0.5 * inner_sum2;
	}

	//Get rho gradients as matrix
	ublas::matrix<double> dFdP_mat = vectorToMatrix(dFdP);

	//Get optimization parameters from rho matrix
	ublas::vector<double> dfdParams = getElements(dFdP_mat);

	for (size_t d = 0; d < gradients.size(); ++d) {
		if (d == gradients.size() - 1) {
			gradients(d) = dFdnu;
		} else {
			gradients(d) = dfdParams(d);
		}
	}

    return -gradients;
}


double T_Copula::dGamma(double t) {
	double dG = tgamma(t) * boost::math::digamma(t);

	return dG;
}

ublas::matrix<double> T_Copula::vectorToMatrix(ublas::vector<double> const& vec) {
	int n = time_series.size2();
	ublas::matrix<double> resMatrix(n,n);
	int counter = 0;

	for (size_t j = 0; j < n; ++j) {
		for (size_t i = 0; i < n; ++i) {
            resMatrix(i, j) = vec(counter);
			counter = counter + 1;
		}
	}

	return resMatrix;
}

ublas::matrix<double> T_Copula::buildP(ublas::vector<double> const& x) {
	int n = time_series.size2();
	ublas::matrix<double> P(n, n);
	int counter = 0;	//keep track of fetched elements

	for (size_t i = 0; i < n; ++i) {
		for (size_t j = 0; j < i + 1; ++j) {
			if (i == j) {
				P(i, j) = 1;
			} else {
				P(i, j) = x(counter);
				P(j, i) = x(counter);
				counter = counter + 1;
			}
		}
	}

	return P;
}


ublas::vector<double> T_Copula::getElements(ublas::matrix<double> const& P) {
	int n = time_series.size2();
	ublas::vector<double> resVec((n-1)*n*0.5); // Vector with optimization parameters in P
	int counter = 0;	//keep track of fetched elements

	for (size_t i = 0; i < n; ++i) {
		for (size_t j = 0; j < i; ++j) {
			resVec(counter) = P(i, j);
			counter = counter + 1;
		}
	}

	return resVec;
}

ublas::vector<double> T_Copula::kronOfVectors(ublas::vector<double> const& v1, ublas::vector<double> const& v2) {
	ublas::vector<double> vKron(v1.size() * v2.size());

	int counter = 0;

	for (size_t i = 0; i < v1.size(); ++i) {
		for (size_t j = 0; j < v2.size(); ++j) {
				vKron(counter) = v1(i)*v2(j);
				counter = counter + 1;
		}
	}

	return vKron;
}

ublas::vector<double> T_Copula::matrixToVector(ublas::matrix<double> const& matrix) {
	int n = matrix.size1();
	ublas::vector<double> resVec(n*n);
	int counter = 0;

	for (size_t j = 0; j < n; ++j) {
		for (size_t i = 0; i < n; ++i) {
			resVec(counter) = matrix(i, j);
			counter = counter + 1;
		}
	}

	return resVec;
}

ublas::vector<double> T_Copula::calcNumGradients(ublas::vector<double> const& x) {

	double epsilon = 2.2 * pow(10, -16);
	ublas::vector<double> increment(sqrt(epsilon) * x);
	ublas::vector<double> num_gradients(x.size());

	/*
	ublas::vector<double> x_0diff = x;
	x_0diff(0) += increment(0);
	ublas::vector<double> x_1diff = x;
	x_1diff(1) += increment(1);
	ublas::vector<double> x_2diff = x;
	x_2diff(2) += increment(2);
	ublas::vector<double> x_3diff = x;
	x_3diff(3) += increment(3);
	ublas::vector<double> x_4diff = x;
	x_4diff(4) += increment(4);

	num_gradients(0) = (function_value(x_0diff) - function_value(x)) / (x_0diff(0) - x(0));
	num_gradients(1) = (function_value(x_1diff) - function_value(x)) / (x_1diff(1) - x(1));
	num_gradients(2) = (function_value(x_2diff) - function_value(x)) / (x_2diff(2) - x(2));
	num_gradients(3) = (function_value(x_3diff) - function_value(x)) / (x_3diff(3) - x(3));
	num_gradients(4) = (function_value(x_4diff) - function_value(x)) / (x_4diff(4) - x(4));
	*/

	return num_gradients;
}

double T_Copula::calcStepSize(ublas::vector<double> const& x, ublas::vector<double> const& d) {
	double a = 1;
	bool accepted = false;
	//Kontrollera att rho är positiv definit, dvs minsta egenvärdet är positivt, annars halvera steglängden.

	while (x(x.size() - 1) + a * d(d.size() - 1) <= 2 || x(x.size() - 1) + a * d(d.size() - 1) > 10) {
		a *= 0.5;
	}
	
	while (!accepted) {
		accepted = true;
        for (size_t i = 0; i < x.size() - 1; ++i) {
            if (x(i) + a * d(i) < -1 || x(i) + a * d(i) > 1) {
				accepted = false;
                break;
			}
		}

		ublas::matrix<double> Pnext = buildP(x + a * d);
		double minEigenvalue = FactorCalculation::smallest_eigval(Pnext);
        //std::cout << "Smallest eigval = " << std::to_string(minEigenvalue) << std::endl;
		if (minEigenvalue <= 0) {
			accepted = false;
		}

		if (isnan(function_value(x + a * d))) {
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