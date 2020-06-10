#include "../../CurveLibrary/sample_handler.h"
#include "../../ParameterEstimation/bfgs.h"
#include "../../ParameterEstimation/Distribution.h"
#include "../../ParameterEstimation/Gaussian.h"
#include "../../ParameterEstimation/Student_t.h"
#include "../../ParameterEstimation/T_Copula.h"
#include "../../ParameterEstimation/Gaussian_Copula.h"

#include <boost/math/distributions/normal.hpp>
#include <boost/math/distributions/students_t.hpp>
#include <boost/container/vector.hpp>

#include <boost/fusion/algorithm/transformation/push_back.hpp>
#include <boost/fusion/include/push_back.hpp>
#include <boost/range/algorithm.hpp>
#include <boost/range/algorithm_ext.hpp>

#include <boost/assign/std/vector.hpp>

#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <strstream>
#include <sstream>
#include <random>
//#include <boost\filesystem\operations.hpp>

matrix<double> read_matrix(std::string const& file_name, int rows, int columns);
matrix<double> delta_curves(matrix<double> curves);
vector<double> read_time_series_test(std::string const& file_name);
matrix<double> gen_start_params(int n, std::string dist);
vector<double> runBFGS_Gaussian(int n, vector<double> series, int max_iter, double epsilon);
vector<double> runBFGS_T(int n, vector<double> series, int max_iter, double epsilon);
vector<double> getUniformTimeseries(vector<double> series, vector<double> params);
vector<double> GARCH_vec(vector<double> time_series, vector<double> x);

int main() {
	/*
	vector<double> time_series;
	matrix<double> hist_rf;
	matrix<double> E_rf;
	
	//Läs in riskfria kurvan
	try {
		//Senaste kurvan sist
		hist_rf = read_matrix("fHist.csv", 3451, 730);
	}
	catch (std::exception & ex) {
		std::cout << "Error:" << ex.what() << "\n";
		return 1;
	}

	//Läs in egenmatrisen för riskfria kurvan
	try {
		//Senaste kurvan sist
		E_rf = read_matrix("EZero.csv", 730, 3);
	}
	catch (std::exception & ex) {
		std::cout << "Error:" << ex.what() << "\n";
		return 1;
	}

	// Beräkna riskfaktorer på riskfria kurvan.
	matrix<double> delta_f = delta_curves(hist_rf);
	matrix<double> hist_riskfactors = prod(trans(E_rf), trans(delta_f));
	hist_riskfactors = trans(hist_riskfactors);

	//Save uniform variables in U matrix for copula estimation
	matrix<double> U(hist_riskfactors.size1(), hist_riskfactors.size2());
	int nRiskfactors = U.size2();

	//Save optimal parameters for all riskfactors
	boost::container::vector<vector<double> > OptParamsAll(0);


	// Optimize parameters for all riskfactors
	int nSolutions = 1;
	int max_iter = 100;
	double epsilon = pow(10, -7);

	for (size_t i = 0; i < nRiskfactors; ++i) {
		matrix_column<matrix<double> > xi(hist_riskfactors, i); //Get risk factor nr i.

		vector<double> optGaussian_xi = runBFGS_Gaussian(nSolutions, xi, max_iter, epsilon);
		vector<double> optStudent_xi = runBFGS_T(nSolutions, xi, max_iter, epsilon);
		std::cout << "OptGaussian = " << optGaussian_xi << "\n\n";
		std::cout << "OptStudent = " << optStudent_xi << "\n\n";

		//Get uniform Timeseries
		if (optGaussian_xi(4) < optStudent_xi(5)) {	//Choose Gaussian marginal dist
		
			matrix_column<matrix<double> > U_column(U, i);
			U_column = getUniformTimeseries(xi, optGaussian_xi); 
			OptParamsAll.push_back(optGaussian_xi);
		} else {										//Choose Student t marginals
			matrix_column<matrix<double> > U_column(U, i);
			U_column = getUniformTimeseries(xi, optStudent_xi);
			OptParamsAll.push_back(optStudent_xi);
		}
	}

	//matrix_column<matrix<double> > U_column(U, 0);
	//std::cout << "U_column = " << U_column << "\n\n";

	std::cout << "OptParams first eigenvector = " << OptParamsAll[0] << "\n\n";
	std::cout << "OptParams second eigenvector = " << OptParamsAll[1] << "\n\n";
	std::cout << "OptParams third eigenvector = " << OptParamsAll[2] << "\n\n";
	*/

	//Läs in riskfria kurvan
    matrix<double> U;
	try {
		//Senaste kurvan sist
		U = read_csv_matrix("U.csv");
	}
	catch (std::exception& ex) {
		std::cout << "Error:" << ex.what() << "\n";
		return 1;
	}

	//vector<double> testU(3);
	//testU(0) = 0.8;
	//testU(1) = 0.24;
	//testU(2) = 0.92;

//	matrix<double> Umat(1, 3);
	//matrix_row<matrix<double> > U_row(Umat, 0);
	//U_row = testU;

	//double FV = T->function_value(P_e);

	//vector<double> gradients = gaussianC->calcGradients(P_e);
	//double FVgaussian = gaussianC->function_value(P_e);
	//std::cout << "FV gaussian = " << gradients << "\n\n";
	
	T_Copula dist(U);
	T_Copula* TC = &dist;

	Gaussian_Copula distG(U);
	Gaussian_Copula* gaussianC = &distG;

	vector<double> t_params(4);
    t_params(0) = 0.2;
    t_params(1) = 0.3;
    t_params(2) = 0.4;
    t_params(3) = 5;


	vector<double> norm_params(4);
    norm_params(0) = 0.2;
    norm_params(1) = 0.3;
    norm_params(2) = 0.4;
	
	identity_matrix<double> I(4);
	int max_iter = 100;
	double epsilon = pow(10, -7);

	vector<double> t_copula_results = bfgs::minimize(t_params, I, max_iter, epsilon, TC);
	vector<double> norm_copula_results = bfgs::minimize(norm_params, I, max_iter, epsilon, gaussianC);

	matrix<double> P_t = TC->buildP(t_copula_results);
	std::cout << "P Students t = " << P_t << "\n\n";
	std::cout << "FV Students t copula = " << t_copula_results(4) << "\n\n";

	matrix<double> P_norm = gaussianC->buildP(norm_copula_results);
	std::cout << "P Students t = " << P_norm << "\n\n";
	std::cout << "FV Students t copula = " << norm_copula_results(4) << "\n\n";
}

vector<double> getUniformTimeseries(vector<double> series, vector<double> params) {
	vector<double> garchVec = GARCH_vec(series, params);
	double mu = params(3);
	vector<double> U(series.size());

	if (params.size() == 5) {
		for (size_t i = 0; i < series.size(); ++i) {
			boost::math::normal norm = boost::math::normal::normal_distribution(0, 1);
			U(i) = cdf(norm,((series(i) - mu) / sqrt(garchVec(i))));
		}
	}
	else {
		boost::math::students_t stud = boost::math::students_t::students_t_distribution(params(4));
		for (size_t i = 0; i < series.size(); ++i) {
			U(i) = cdf(stud, ((series(i) - mu) / sqrt(garchVec(i))));
		}
	}

	return U;
}



vector<double> GARCH_vec(vector<double> time_series, vector<double> x) {
	vector<double> garch_vec(time_series.size() + 1);

	//Calc variance for timeseries
	boost::accumulators::accumulator_set<double, boost::accumulators::stats<boost::accumulators::tag::variance> > acc;
	for_each(time_series.begin(), time_series.end(), boost::bind<void>(boost::ref(acc), _1));

	//Set variance as first element
	garch_vec(0) = boost::accumulators::variance(acc);

	for (size_t i = 0; i < time_series.size(); i++) {
		garch_vec(i + 1) = x(0) + x(1) * pow(time_series(i), 2) + x(2) * garch_vec(i);
	}

	return garch_vec;
}

vector<double> runBFGS_Gaussian(int n, vector<double> series, int max_iter, double epsilon) {
	//Convert vector as matrix for BFGS minimize function
	matrix<double> timeseries(series.size(), 1);
	matrix_column<matrix<double> > mc(timeseries, 0);
	mc = series;

	//Create Student_t object and set Hessian
	Gaussian dist(timeseries);
	Gaussian* distribution = &dist;
	matrix<double> H_inv(4, 4);
	identity_matrix<double> I(4);
	H_inv = I;
	matrix<double> params(4, n);
	vector<double> FV(n);

	params = gen_start_params(n, "normal");
	// Run optimization problem n times.
	for (size_t i = 0; i < params.size2(); ++i) {
		vector<double> results = bfgs::minimize(column(params, i), H_inv, max_iter, epsilon, distribution);

		for (size_t j = 0; j < params.size1(); ++j) {
			column(params, i)(j) = results(j);
		}
		FV(i) = results(params.size1());
	}

	//Get smallest loglikelihood value and corresponding parameters
	double smallest = FV(0);
	int index = 0;
	for (vector<double>::iterator it = FV.begin(); it != FV.end(); it++) {

		if (*it < smallest) {
			smallest = *it;
			index = it - FV.begin();
		}
	}

	std::cout << "All FV:S  : " << FV << "\n\n";
	std::cout << "OPT FV = " << smallest << " at params = " << column(params, index) << "\n\n";

	//Save optimal parameters and function value in opt_params
	vector<double> opt_params(5);
	opt_params(0) = params(0, index);
	opt_params(1) = params(1, index);
	opt_params(2) = params(2, index);
	opt_params(3) = params(3, index);
	opt_params(4) = smallest;

	return opt_params;
}

vector<double> runBFGS_T(int n, vector<double> series, int max_iter, double epsilon) {
	//Convert vector as matrix for BFGS minimize function
	matrix<double> timeseries(series.size(), 1);
	matrix_column<matrix<double> > mc(timeseries, 0);
	mc = series;

	//Create Student_t object and set Hessian
	Student_t dist_t(timeseries);
	Student_t* distribution_t = &dist_t;
	matrix<double> H_inv(5, 5);
	identity_matrix<double> I(5);
	H_inv = I;
	matrix<double> params(5, n);
	vector<double> FV(n);

	params = gen_start_params(n, "t");
	// Run optimization problem n times.
	for (size_t i = 0; i < params.size2(); ++i) {
		vector<double> results = bfgs::minimize(column(params, i), H_inv, max_iter, epsilon, distribution_t);

		for (size_t j = 0; j < params.size1(); ++j) {
			column(params, i)(j) = results(j);
		}
		FV(i) = results(params.size1());
	}

	//Get smallest loglikelihood value and corresponding parameters
	double smallest = FV(0);
	int index = 0;
	for (vector<double>::iterator it = FV.begin(); it != FV.end(); it++) {

		if (*it < smallest) {
			smallest = *it;
			index = it - FV.begin();
		}
	}

	std::cout << "All FV:S Student t: " << FV << "\n\n";
	std::cout << "OPT FV Student T= " << smallest << " at params = " << column(params, index) << "\n\n";
	//Save optimal parameters and function value in opt_params
	vector<double> opt_params(6);
	opt_params(0) = params(0, index);
	opt_params(1) = params(1, index);
	opt_params(2) = params(2, index);
	opt_params(3) = params(3, index);
	opt_params(4) = params(4, index);
	opt_params(5) = smallest;
	return opt_params;
}


matrix<double> gen_start_params(int n, std::string dist){
	int nParams = 4;
	if (dist == "normal") {
		nParams = 4;
	}
	else if (dist == "t") {
		nParams = 5;
	}
	vector<double> params(nParams);
	matrix<double> param_matrix(nParams, n);

	std::default_random_engine generator;
	std::uniform_real_distribution<double> distribution(0.7,0.95);

	for (size_t i = 0; i < n; ++i) {
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


vector<double> read_time_series_test(std::string const& file_name) {
	std::fstream fin;
	std::string line;
	std::string value_string;
	double value = 0;

	fin.open(file_name, std::ios::in);

	if (!fin) {
		throw std::runtime_error("Could not open file");
	}

	vector<double> series(7981);
	int k = 0;

	while (getline(fin, line)) {
		std::stringstream s(line);

		//Get value
		getline(s, value_string, '\n');
		const char* decimal = value_string.c_str();  // String to const char
		value = std::strtod(decimal,NULL); // const char to float.

		//series(k) = value;
		series(series.size() - 1 - k) = value;
		k = k + 1;
	}
	fin.close();

	vector<double> log_returns(series.size() - 1);
	for (size_t i = 0; i < series.size() - 1; ++i) {
		log_returns(i) = log(series(i + 1) / series(i));
	}

	return log_returns;
}


matrix<double> read_matrix(std::string const& file_name, int rows, int columns) {

	std::fstream fin;
	std::string line;
	std::string value_string;
	double value = 0;
	matrix<double> matrix(rows, columns);

	fin.open(file_name, std::ios::in);

	//Read curves from file
	//One row is one curve. Number of columns is number of discretization points
	if (!fin) {
		throw std::runtime_error("Could not open file");
	}

	for (size_t i = 0; i < matrix.size1(); ++i) {
	getline(fin, line);
		std::stringstream s(line);

		for (size_t j = 0; j < matrix.size2(); ++j) {
			(getline(s, value_string, ','));
			const char* decimal = value_string.c_str();  // String to const char
			value = std::atof(decimal); // const char to float.
			matrix(i,j) = value;
		}
	}
	fin.close();

	return matrix;
}

matrix<double> delta_curves(matrix<double> curves) {
	//Calculate delta of curves from the curves matrix
	matrix<double> delta_f(curves.size1() - 1, curves.size2());

	for (size_t i = 0; i < delta_f.size1(); ++i) {

		for (size_t j = 0; j < delta_f.size2(); ++j) {
			delta_f(i, j) = curves(i + 1, j) - curves(i, j);
		}
	}

	return delta_f;
}


