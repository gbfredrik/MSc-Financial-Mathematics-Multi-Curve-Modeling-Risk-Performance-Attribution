#include "bfgs.h"


vector<double> bfgs::minimize(vector<double> x, matrix<double> H_inv, int max_iter, float epsilon, Distribution* dist, double dt) {

	std::cout << "in minimize: ";
	dist->getSeries();


	int n = x.size();
	matrix<double> help_prod(n, n);
	vector<double> d(n);
	double alpha;
	int k = 0;
	vector<double> x_new(n);
	vector<double> y(n);
	vector<double> s(n);
	double l;
	vector<double> gradient_vec(n);
	//Init hessian inverse matrix
	identity_matrix<double> I(n);

	
	
	
	while (norm_2(dist->calcGradients(x, dt)) > epsilon && k < max_iter) {

		gradient_vec = dist->calcGradients(x, dt);
		std::cout << "H_inv: " << H_inv << "\n";
		std::cout << "gradient_vec: " << gradient_vec << "\n \n";
		d = -prod(H_inv, gradient_vec);
		
		alpha = dist->calcStepSize(x, d, dt);
		x_new = x + alpha * d;
		
		y = dist->calcGradients(x_new, dt) - dist->calcGradients(x, dt);
		s = x_new - x;
		l = 1 / inner_prod(y, s);
		
		help_prod = prod((I - l * outer_prod(s, y)), H_inv);
		H_inv = prod(help_prod, (I - l * outer_prod(y, s))) + l * outer_prod(s, s);
		x = x_new;
		k = k + 1;
		
		std::cout << "k = " << k << " ,Function value = " << dist->function_value(x, dt) << " for parameters : " << x  
					<< " and norm of gradients = " << norm_2(dist->calcGradients(x, dt)) << "\n";
	}

	
	return x;
}

