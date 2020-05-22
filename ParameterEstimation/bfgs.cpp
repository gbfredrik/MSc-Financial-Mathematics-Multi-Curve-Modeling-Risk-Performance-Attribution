#include "bfgs.h"

int bfgs::function_type;

vector<double> vec(1);
Distribution bfgs::dist(vec);


vector<double> bfgs::minimize(vector<double> x, matrix<double> H_inv, int max_iter, float epsilon, Distribution dist) {

	bfgs::dist = dist;
	int n = x.size();
	matrix<double> help_prod(n, n);
	vector<double> d;
	double alpha;
	int k = 0;
	vector<double> x_new(n);
	vector<double> y(n);
	vector<double> s(n);
	double l;
	vector<double> gradient_vec(n);

	//Init hessian inverse matrix
	identity_matrix<double> I(n);


	
	while (norm_2(calcGradients(x)) > epsilon && k < max_iter) {

		gradient_vec = calcGradients(x);
		d = -prod(H_inv, gradient_vec);

		alpha = calcStepSize(x, d);
		x_new = x + alpha * d;

		y = calcGradients(x_new) - calcGradients(x);
		s = x_new - x;
		l = 1 / inner_prod(y, s);

		help_prod = prod((I - l * outer_prod(y, s)), H_inv);
		H_inv = prod(help_prod, (I - l * outer_prod(s, y))) + l * outer_prod(y, y);
		x = x_new;
		k = k + 1;
	}

	std::cout << "Function value = " << f(x) << "for parameters : " << x << "\n";
	
	return x;
}



double bfgs::calcStepSize(vector<double> x, vector<double> d) {

	double a = 1;
	float c1 = pow(10, -4);
	float c2 = 0.9;

	while (f(x + a * d) > f(x) + c1 * a * inner_prod(calcGradients(x), d) || inner_prod(calcGradients(x + a * d), d) < c2 * inner_prod(calcGradients(x), d))
	{
		a = a * 0.5;
	}
std::cout << a << "\n";
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