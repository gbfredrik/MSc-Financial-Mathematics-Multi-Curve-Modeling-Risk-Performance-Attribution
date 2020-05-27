#include "../../ParameterEstimation/bfgs.h"
#include "../../ParameterEstimation/Distribution.h"
#include "../../ParameterEstimation/Gaussian.h"
#include "../../ParameterEstimation/Student_t.h"

#include <boost/fusion/algorithm/transformation/push_back.hpp>
#include <boost/fusion/include/push_back.hpp>

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

int main() {
	 
	

	//FÖR XI1
	//omega
	//start(0) = 0.001;
	//alpha
	//start(1) = 0.05;
	//beta
	//start(2) = 0.9;
	//väntevärde
	//start(3) = 0;

	////omega xi2
	//start(0) = 0.00000001;
	////alpha
	//start(1) = 0.0496;
	////beta
	//start(2) = 0.950399;
	////väntevärde
	//start(3) = -0.00003199;

	////omega
	//start(0) = 0.001;
	////alpha
	//start(1) = 0.05;
	////beta
	//start(2) = 0.9;
	////väntevärde
	//start(3) = 0;


	/*
	vector<double> start_r(2);
	start_r(0) = 5;
	start_r(1) = -5;
	vector<double> vec(2);
	vec(0) = 0.5;
	vec(1) = 0.5;
	Distribution dist(vec);
	Distribution* rosenbrock = &dist;
	*/

	double dt = 1.00;
	vector<double> time_series;
	matrix<double> hist_rf;
	matrix<double> E_rf;

	
	

	std::cout << "dt = " << dt << "\n";
	
	try {
		//Senaste kurvan sist
		time_series = read_time_series_test("aapl.csv");
	}

	catch (std::exception & ex) {
		std::cout << "Error:" << ex.what() << "\n";
		return 1;
	}
	
	//Läs in riskfria kurvan
	try {
		//Senaste kurvan sist
		hist_rf = read_matrix("fHist.csv", 3451, 730);
	}

	catch (std::exception & ex) {
		std::cout << "Error:" << ex.what() << "\n";
		return 1;
	}
	// Beräkna delta f på riskfria kurvan.
	matrix<double> delta_f = delta_curves(hist_rf);

	//Läs in egenmatrisen för riskfria kurvan
	try {
		//Senaste kurvan sist
		E_rf = read_matrix("EZero.csv", 730, 3);
	}

	catch (std::exception & ex) {
		std::cout << "Error:" << ex.what() << "\n";
		return 1;
	}


	matrix<double> hist_risk_faktors = prod(trans(E_rf), trans(delta_f));


	matrix_row<matrix<double> > xi1(hist_risk_faktors, 0);
	matrix_row<matrix<double> > xi2(hist_risk_faktors, 1);
	matrix_row<matrix<double> > xi3(hist_risk_faktors, 2);


	
	//vector<double> start(4);
	
	
	
	
	
	// Create start parameters 


	std::string dist_choice = "t";
	int nSolutions = 20;
	vector<double> series = xi1;
	int max_iter = 100;
	double epsilon = pow(10, -7);

	//Initiate gaussian
	int nParams = 4;
	Gaussian dist(series);
	Gaussian* distribution = &dist;
	
	Student_t dist_t(series);
	Student_t* distribution_t = &dist_t;


	if (dist_choice == "t") {
		nParams = 5;
		
	}

	matrix<double> H_inv(nParams, nParams);
	identity_matrix<double> I(nParams);
	H_inv = I;
	
	matrix<double> params(nParams, nSolutions);
	boost::numeric::ublas::vector<double> FV(nSolutions);
	matrix<double> opt_params(nParams, nSolutions);
	
	//params = gen_start_params(nSolutions, dist_choice);

	vector<double> test_x(5);
	test_x(0) = 0.00000005379722106733691;
	test_x(1) = 0.043558331112312;
	test_x(2) = 0.935302736808228;
	test_x(3) = -0.0001904288328934271;
	test_x(4) = 5.269578291031785;

	std::cout << "Test LL = " << distribution_t->function_value(test_x);
	/*
	for (size_t i = 0; i < params.size2(); ++i) {

		vector<double> results = bfgs::minimize(column(params,i), H_inv, max_iter, epsilon, distribution_t);
		
		for (size_t j = 0; j < opt_params.size1(); ++j) {
			column(opt_params,i)(j) = results(j);
		}
		FV(i) = results(params.size1());
	}

	double smallest = FV(0);
	int index = 0;

	for (vector<double>::iterator it = FV.begin(); it != FV.end() ; it++) {

		if (*it < smallest) {
			smallest = *it;
			index = it - FV.begin();
		}
	}


	
	std::cout << "OPT FV = " << smallest << " at params = " << column(opt_params, index) <<  "\n\n";
	*/
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

	//matrix_row<matrix<float> > mr(delta_f, 0);

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


