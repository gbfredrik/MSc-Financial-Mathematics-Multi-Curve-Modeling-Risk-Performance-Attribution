#include "../../ParameterEstimation/bfgs.h"
#include "../../ParameterEstimation/Distribution.h"
#include "../../ParameterEstimation/Gaussian.h"

#include <boost/fusion/algorithm/transformation/push_back.hpp>
#include <boost/fusion/include/push_back.hpp>

#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <strstream>
#include <sstream>

matrix<double> read_matrix(std::string const& file_name, int rows, int columns);
matrix<double> delta_curves(matrix<double> curves);
vector<double> read_time_series_test(std::string const& file_name);


int main() {
	 
	vector<double> start(4);
	
	/*omega
	start(0) = 0.001;
	//alpha
	start(1) = 0.05;
	//beta
	start(2) = 0.9;
	//väntevärde
	start(3) = 0.05;
	*/
	/* FÖR XI1
	//omega
	start(0) = 0.001;
	//alpha
	start(1) = 0.05;
	//beta
	start(2) = 0.9;
	//väntevärde
	start(3) = 0;*/

	//omega xi2
	//start(0) = 0.00000001;
	//alpha
	//start(1) = 0.0496;
	//beta
	//start(2) = 0.950399;
	//väntevärde
	//start(3) = -0.00003199;

	//omega
	start(0) = 0.00000001;
	//alpha
	start(1) = 0.0496;
	//beta
	start(2) = 0.950399;
	//väntevärde
	start(3) = -0.00003199;

	vector<double> start_r(2);
	start_r(0) = 5;
	start_r(1) = -5;

	int max_iter = 200;
	double epsilon = pow(10,-7);
	double dt = 1.00;
	matrix<double> H_inv(4, 4);
	identity_matrix<double> I(4);
	H_inv = I;
	vector<double> time_series;
	matrix<double> hist_rf;
	matrix<double> E_rf;
	
	vector<double> vec(2);
	vec(0) = 0.5;
	vec(1) = 0.5;
	Distribution dist(vec);
	Distribution* rosenbrock = &dist;
	
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
	//std::cout << "delta_f = " << trans(delta_f) << "\n \n";

	matrix<double> hist_risk_faktors = prod(trans(E_rf), trans(delta_f));


	matrix_row<matrix<double> > xi1(hist_risk_faktors, 0);
	matrix_row<matrix<double> > xi2(hist_risk_faktors, 1);
	matrix_row<matrix<double> > xi3(hist_risk_faktors, 2);

	//std::cout << "Riskfaktor 1 = " << xi1 << "\n\n";
	//std::cout << "Time_series = " << time_series << "\n\n";

	

	Gaussian dist2(xi3);
	Gaussian* gaussian = &dist2;



	//std::cout << "Funktionsvärde startparametrar : " << gaussian->function_value(start, dt) << "\n";


	vector<double> opt_parameters = bfgs::minimize(start, H_inv, max_iter, epsilon, gaussian);

	//std::cout << opt_parameters;



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


