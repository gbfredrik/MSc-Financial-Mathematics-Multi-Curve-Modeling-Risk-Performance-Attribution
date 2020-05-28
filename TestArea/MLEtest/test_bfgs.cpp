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
	
	vector<double> time_series;
	matrix<double> hist_rf;
	matrix<double> E_rf;
	
	
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
	matrix<double> hist_risk_faktors = prod(trans(E_rf), trans(delta_f));

	matrix_row<matrix<double> > xi1(hist_risk_faktors, 0);
	matrix_row<matrix<double> > xi2(hist_risk_faktors, 1);
	matrix_row<matrix<double> > xi3(hist_risk_faktors, 2);

	// OPtimization problem options

	std::string dist_choice = "t";
	matrix<double> series(xi1.size(), 1);
	matrix_column<matrix<double> > mc(series, 0);
	mc = xi1;

	//vector<double> series = xi1;
	int nSolutions = 10;
	int max_iter = 100;
	double epsilon = pow(10, -7);

	//Initiate distribution and 
	int nParams = 4;
	//Gaussian dist(series);
	//Gaussian* distribution = &dist;
	
	Student_t dist_t(series);
	Student_t* distribution_t = &dist_t;


	if (dist_choice == "t") {
		nParams = 5;
	}

	matrix<double> H_inv(nParams, nParams);
	identity_matrix<double> I(nParams);
	H_inv = I;
	
	matrix<double> params(nParams, nSolutions);
	vector<double> FV(nSolutions);
	//matrix<double> opt_params(nParams, nSolutions);
	
	params = gen_start_params(nSolutions, dist_choice);

	for (size_t i = 0; i < params.size2(); ++i) {

		vector<double> results = bfgs::minimize(column(params,i), H_inv, max_iter, epsilon, distribution_t);
		
		for (size_t j = 0; j < params.size1(); ++j) {
			column(params,i)(j) = results(j);
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

	std::cout << "All FV:S : " << FV << "\n\n";
	std::cout << "OPT FV = " << smallest << " at params = " << column(params, index) <<  "\n\n";

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


