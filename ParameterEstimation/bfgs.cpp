#include "bfgs.h"

int bfgs::function_type;

vector<double> vec(1);
Distribution bfgs::dist(vec);


vector<double> bfgs::minimize(vector<double> x, matrix<double> H_inv, int max_iter, float epsilon, Distribution* dist) {

	std::cout << "in minimize: ";
	dist->getSeries();


	//bfgs::dist = dist;
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

	
	
	
	while (norm_2(dist->calcGradients(x)) > epsilon && k < max_iter) {
	//while ( k < max_iter) {
		gradient_vec = dist->calcGradients(x);
		std::cout << "H_inv: " << H_inv << "\n";
		std::cout << "gradient_vec: " << gradient_vec << "\n \n";
		d = prod(H_inv, gradient_vec);
		
		alpha = calcStepSize(x, d, dist);
		x_new = x + alpha * d;
		
		y = dist->calcGradients(x_new) - dist->calcGradients(x);
		s = x_new - x;
		l = 1 / inner_prod(y, s);
		
		help_prod = prod((I - l * outer_prod(y, s)), H_inv);
		H_inv = prod(help_prod, (I - l * outer_prod(s, y))) + l * outer_prod(y, y);
		x = x_new;
		k = k + 1;
		
		std::cout << "k = " << k << " ,Function value = " << dist->function_value(x) << " for parameters : " << x << "\n \n";
	}

	
	return x;
}



double bfgs::calcStepSize(vector<double> x, vector<double> d, Distribution* dist) {

	double a = 1;
	double c1 = pow(10, -4);
	double c2 = 0.9;

	while (x(0) + a * d(0) < 0 || x(1) + a * d(1) < 0 || x(2) + a * d(2) < 0 || x(1)+a*d(1) + x(2)+a*d(2) >= 1) {
	
		std::cout << "old parameters = " << x << "\n";
		std::cout << "new parameters = " << x(0) + a * d(0) << ", " << x(1) + a * d(1) << ", " << x(2) + a * d(2) << "\n";
		a = a * 0.5;
		std::cout << "steglängd = " << a << "\n";
	}

	std::cout << " old parameters = " << x << "\n";
	std::cout << "final new parameters = " << x(0) + a * d(0) << ", " << x(1) + a * d(1) << ", " << x(2) + a * d(2) << "\n";
	std::cout << "steglängd = " << a << "\n";
	


	//while (dist->function_value(x + a * d) > dist->function_value(x) + c1 * a * inner_prod(dist->calcGradients(x), d)
	while (dist->function_value(x + a * d) > dist->function_value(x))
	{
		std::cout << "new f : " << dist->function_value(x + a * d) << "\n";
		std::cout << "old f : " << dist->function_value(x) + c1 * a * inner_prod(dist->calcGradients(x), d) << "\n";
		std::cout << "steglängd = " << a << "\n";
		a = a * 0.5;
	}

	std::cout << "new f : " << dist->function_value(x + a * d) << "\n";
	std::cout << "old f : " << dist->function_value(x) + c1 * a * inner_prod(dist->calcGradients(x), d) << "\n";
	std::cout << "steglängd = " << a << "\n \n";

	return a;
}



double bfgs::f(vector<double> params) {
	double fun_value;

	fun_value = bfgs::dist.function_value(params);

	return fun_value;
}

vector<double> bfgs::calcGradients(vector<double> x) {

	vector<double> gradients = bfgs::dist.calcGradients(x);
	return gradients;

}

/*
double bfgs::f(vector<double> params) {
	double fun_value;

	switch (function_type) {
	case 1: fun_value = bfgs::rosenbrock(params);
	}

	return fun_value;
}

double bfgs::rosenbrock(vector<double> x) {

	double rosenbrock_value = pow(1 - x(0), 2) + 100 * pow(x(1) - x(0) * x(0), 2);
	return rosenbrock_value;
}

vector<double> bfgs::calcGradients(vector<double> x) {

	double dx = 2 * (200 * pow(x(0), 3) - 200 * x(0) * x(1) + x(0) - 1);
	double dy = 200 * (x(1) - pow(x(0), 2));

	vector<double> gradients(2);
	gradients(0) = dx;
	gradients(1) = dy;
	return gradients;

}

*/