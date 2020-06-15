#include "bfgs.h"


vector<double> bfgs::minimize(boost::numeric::ublas::vector<double> x, matrix<double> H_inv, int max_iter, double epsilon, Distribution* dist) {

	//std::cout << "in minimize: ";
	//dist->getSeries();

	int n = x.size();
	matrix<double> help_prod(n, n);
	vector<double> d(n);
	double alpha;
	int k = 0;
	vector<double> x_new(n);
	vector<double> y(n);
	vector<double> s(n);
	double scale_H;
	double l;
	vector<double> gradient_vec(n);
	//Init hessian inverse matrix
	identity_matrix<double> I(n);
	
	while (norm_2(dist->calcGradients(x)) > epsilon && k < max_iter) {
		std::cout << "New iteration \n\n";
		gradient_vec = dist->calcGradients(x);
		std::cout << "    H_inv: " << H_inv << "\n\n";
		//std::cout << "Num H_inv: " << dist->calcNumHessian(x) << "\n\n";
		std::cout << "gradient_vec: " << gradient_vec << "\n \n";
		//std::cout << "Numerical gradient_vec: " << dist->calcNumGradients(x) << "\n \n";
		d = -prod(H_inv, gradient_vec);

		//std::cout << "d : " << d << "\n\n";
		alpha = dist->calcStepSize(x, d);
		x_new = x + alpha * d;

		y = dist->calcGradients(x_new) - dist->calcGradients(x);
		s = x_new - x;
		scale_H = inner_prod(y, s);
		l = 1 /scale_H;

		//std::cout << "y: " << y << "\n\n";
		//std::cout << "s: " << s << "\n\n";
		//std::cout << "l: " << l << "\n\n";

		if (scale_H == 0) {
			std::cout << "Break since H has invalid values. \n\n";
			break;
		}

		help_prod = prod((I - l * outer_prod(s, y)), H_inv);
		H_inv = prod(help_prod, (I - l * outer_prod(y, s))) + l * outer_prod(s, s);
		x = x_new;
		k = k + 1;
		
		std::cout << "k = " << k << " , Function value = " << dist->function_value(x) << " for parameters : " << x  
					<< " and norm of gradients = " << norm_2(dist->calcGradients(x)) << "\n \n \n\n";
	}

	std::cout << "Final k = " << k << " , Function value = " << dist->function_value(x) << " for parameters : " << x
		<< " and norm of gradients = " << norm_2(dist->calcGradients(x)) << "\n \n \n\n";
	
	vector<double> results(n+1);

	
	for (size_t i = 0; i < results.size() - 1; i++) {
		results(i) = x(i);
	}

	results(n) = dist -> function_value(x);

	return results;
}

