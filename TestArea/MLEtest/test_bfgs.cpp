#include "../../ParameterEstimation/bfgs.h"
#include "../../ParameterEstimation/Distribution.h"

int main() {

	vector<double> start(2);
	for (unsigned i = 0; i < start.size(); ++i)
		start(i) = 0.8;

	int max_iter = 500;
	float epsilon = 10 ^ (-4);
	matrix<double> H_inv(2, 2);
	identity_matrix<double> I(2);
	H_inv = I;
	
	vector<double> vec(2);
	vec(0) = 1;
	vec(1) = 1;
	Distribution dist(vec);
	//Distribution dist;

	vector<double> opt_parameters = bfgs::minimize(start, H_inv, max_iter, epsilon, dist);

	std::cout << opt_parameters;



}