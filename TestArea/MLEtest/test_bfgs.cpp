#include "../../ParameterEstimation/bfgs.h"
#include "../../ParameterEstimation/Distribution.h"
#include "../../ParameterEstimation/Gaussian.h"

vector<double > simulate_GARCH_process(int n);

int main() {

	vector<double> start(2);
	//start(0) = 0.02;
	//start(1) = 0.04;
	//start(2) = 0.95;
	start(0) = 5;
	start(1) = 5;
	//start(2) = 0.95;

	int max_iter = 100;
	float epsilon = pow(10,-5);
	matrix<double> H_inv(2, 2);
	identity_matrix<double> I(2);
	H_inv = I;
	
	
	vector<double> vec(2);
	vec(0) = 0.5;
	vec(1) = 0.5;
	Distribution dist(vec);
	Distribution* d = &dist;

	vector<double> time_series = simulate_GARCH_process(10);

	Gaussian dist2(time_series);

	Gaussian* gaussian = &dist2;

	std::cout << "\nstartit = " << start << "\n";


	std::cout << "Funktionsvärde startparametrar : " << d->function_value(start) << "\n";

	vector<double> opt_parameters = bfgs::minimize(start, H_inv, max_iter, epsilon, d);

	//std::cout << opt_parameters;



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
	
	/**
	for (int i = 0; i < n; i++) {

		garch(i) = sqrt(nu(i));
	}
	*/
	
	//vector<double> process = element_prod(nu, normals);

	std::cout << process;

	return process;
}