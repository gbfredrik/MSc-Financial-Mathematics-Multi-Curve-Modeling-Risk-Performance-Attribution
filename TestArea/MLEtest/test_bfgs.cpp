#include "../../ParameterEstimation/bfgs.h"
#include "../../ParameterEstimation/Distribution.h"
#include "../../ParameterEstimation/Gaussian.h"

vector<double > simulate_GARCH_process(int n);

int main() {

	vector<double> start(3);
	for (unsigned i = 0; i < start.size(); ++i)
		start(i) = 0.4;

	int max_iter = 4;
	float epsilon = 10 ^ (-4);
	matrix<double> H_inv(3, 3);
	identity_matrix<double> I(3);
	H_inv = -I;
	
	
	vector<double> vec(2);
	vec(0) = 0.5;
	vec(1) = 0.5;
	Distribution dist(vec);
	

	vector<double> time_series = simulate_GARCH_process(10);

	Gaussian dist2(time_series);

	Gaussian* gaussian = &dist2;

	std::cout << "startit = " << start << "\n";

	std::cout << "Funktionsvärde startparametrar : " << gaussian->function_value(start) << "\n";

	vector<double> opt_parameters = bfgs::minimize(start, H_inv, max_iter, epsilon, &dist2);

	std::cout << opt_parameters;



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

	normals(0) = -0.532011376808821;
	normals(1) = 1.68210359466318;
	normals(2) = -0.875729346160017;
	normals(3) = -0.483815050110121;
	normals(4) = -0.712004549027423;
	normals(5) = -1.17421233145682;
	normals(6) = -0.192239517539275;
	normals(7) = -0.274070229932602;
	normals(8) = 1.53007251442410;
	normals(9) = -0.249024742513714;


	double w = 0.0099;
	double alpha = 0.04;
	double beta = 0.95;

	vector<double> garch(n);

	garch(0) = pow(normals(0), 2);

	for (int i=1; i < n; i++) {
		garch(i) = w + alpha * pow(normals(i),2) + beta * garch(i - 1);
	}

	std::cout << "In garch normals: " << normals << "\n";
	std::cout << "In garch garch: " << garch << "\n";


	for (int i = 0; i < n; i++) {
		garch(i) = sqrt(garch(i));
	}

	vector<double> process = element_prod(garch, normals);

	std::cout << process;

	return process;
}