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


vector<float> read_time_series(std::string const& file_name);
vector<double> get_GARCH_process();
vector<double> simulate_GARCH_process(int n);

int main() {
	 
	vector<double> start(4);
	start(0) = 0.001;
	start(1) = 0.05;
	start(2) = 0.9;
	start(3) = 0.05;
	//start(0) = 5;
	//start(1) = 5;
	//start(2) = 0.95;

	int max_iter = 100;
	float epsilon = pow(10,-5);
	matrix<double> H_inv(4, 4);
	identity_matrix<double> I(4);
	H_inv = I;
	
	
	vector<double> vec(2);
	vec(0) = 0.5;
	vec(1) = 0.5;
	Distribution dist(vec);
	Distribution* rosenbrock = &dist;
	
	vector<float> time_series;

	try {
		//time_series = get_GARCH_process2();
		time_series = read_time_series("aapl.csv");
	}

	catch (std::exception & ex) {
		std::cout << "Error:" << ex.what() << "\n";
		return 1;
	}

	std::cout << "timeseries = " << time_series << "\n";

	vector<float> time_series2(1000);
	for (size_t i = 0; i < time_series2.size(); ++i) {
		time_series2(i) = time_series(i);
	}



	Gaussian dist2(time_series);

	Gaussian* gaussian = &dist2;

	std::cout << "\n startit = " << start << "\n";


	std::cout << "Funktionsvärde startparametrar : " << gaussian->function_value(start) << "\n";

	vector<double> theta(3);
	theta(0) = 0.02;
	theta(1) = 0.04;
	theta(2) = 0.95;
	//std::cout << "Optimalt funktionsvärde : " << gaussian->function_value(theta) << "\n\n";

	vector<double> opt_parameters = bfgs::minimize(start, H_inv, max_iter, epsilon, gaussian);

	//std::cout << opt_parameters;



}

vector<float> read_time_series(std::string const& file_name) {

	std::fstream fin;
	std::string line;
	std::string value_string;
	float value = 0;

	fin.open(file_name, std::ios::in);

	if (!fin) {
		throw std::runtime_error("Could not open file");
	}

	vector<float> series(7981);
	int date;
	int k = 0;


	while (getline(fin, line)) {

		std::stringstream s(line);

		//Get value
		getline(s, value_string, '\n');
		const char* decimal = value_string.c_str();  // String to const char
		value = std::atof(decimal); // const char to float.
		//entry.second = value;

		//series(k) = value;
		series(series.size() - 1 - k) = value;
		k = k + 1;
	}

	fin.close();

	vector<double> log_returns(series.size() - 1);
	for (int i = 0; i < series.size() - 1; ++i) {
		log_returns(i) = log(series(i + 1) / series(i));
	}


	return log_returns;
}




vector<double > get_GARCH_process() {

	std::fstream fin;
	std::string line;
	std::string value_string;
	float value = 0;

	fin.open("MLE_Data.csv", std::ios::in);

	if (!fin) {
		throw std::runtime_error("Could not open file");
	}
	//7981
	vector<double> values(10);




	//while (getline(fin, line)) {
	for (int i = 0; i< 10;++i) {
		(getline(fin, line));
		std::cout << "string read : " << line << "\n";

		std::stringstream s(line);

		getline(s, value_string, '\n');
		std::cout << "string read : " << value_string << "\n";
		const char* decimal = value_string.c_str();  // String to const char
		value = std::atof(decimal); // const char to float.
		//values.push_back(value);
		values(values.size()-1-i) = value;
	}

	fin.close(); 


	
	vector<double> log_returns(values.size() - 1);
		for (int i = 0; i < values.size()-1; ++i) {
			log_returns(i) = log(values(i+1)/values(i));
		}

	

	vector<double> temp(2);
	temp(0) = 1;
	temp(1) = 2;
	return log_returns;
}





vector<double > simulate_GARCH_process(int n) {
	
	vector<double> normals(n);
	/*
	//boost::mt19937 rng; // I don't seed it on purpouse (it's not relevant)
	typedef boost::mt19937 base_generator_type;

	base_generator_type generator(42u);
	


	boost::normal_distribution<> nd(0.0, 1.0);

	boost::variate_generator<boost::mt19937&,
	boost::normal_distribution<> > var_nor(generator, nd);

	for (	int i = 0; i < n; ++i)
	{
		normals(i) = var_nor();
		//generator.seed(static_cast<unsigned int>(std::time(0)));
	}

	*/

	normals(0) = 1;
	normals(1) = 1;
	normals(2) = 1;
	normals(3) = 1;
	normals(4) = 1;
	normals(5) = 1;
	normals(6) = 1;
	normals(7) = 1;
	normals(8) = 1;
	normals(9) = 1;
	//normals(0) = -0.532011376808821;
	//normals(1) = 1.68210359466318;
	//normals(2) = -0.875729346160017;
	//normals(3) = -0.483815050110121;
	//normals(4) = -0.712004549027423;
	//normals(5) = -1.17421233145682;
	//normals(6) = -0.192239517539275;
	//normals(7) = -0.274070229932602;
	//normals(8) = 1.53007251442410;
	//normals(9) = -0.249024742513714;

	double w = 0.02;
	double alpha = 0.04;
	double beta = 0.95;

	vector<double> nu(n);

	nu(0) = pow(normals(0), 2);

	for (int i = 1; i < n; ++i) {
		nu(i) = w + alpha * nu(i - 1) * pow(normals(i), 2) + beta * nu(i - 1);
	}

	std::cout << "In garch normals: " << normals << "\n";
	std::cout << "In garch nu: " << nu << "\n";

	
	vector<double> process(n - 1);
	for (int i = 0; i < n - 1; ++i) {

		process(i) = sqrt(nu(i)) * normals(i + 1);
	}

	std::cout << "In garch process: " << process << "\n";
	
	/**
	for (int i = 0; i < n; i++) {

		garch(i) = sqrt(nu(i));
	}
	*/
	
	//vector<double> process = element_prod(nu, normals);

	std::cout << process;

	return process;
}