#include "T_Copula.h"


T_Copula::T_Copula(matrix<double> series) : Distribution(series) {
	time_series = series;
}


void T_Copula::getSeries() {
	std::cout << time_series << "\n";
}


double T_Copula::function_value(matrix<double> const& x) {
	std::cout << "in FN copula" << "\n\n";
	double sum = 0;

	return -sum;
}


vector<double> T_Copula::calcGradients(vector<double> const& x) {
	vector<double> gradients(x.size());

	return gradients;
}

vector<double> T_Copula::calcNumGradients(vector<double> const& x) {

	double epsilon = 2.2 * pow(10, -16);
	vector<double> increment(sqrt(epsilon) * x);
	vector<double> num_gradients(x.size());

	/*
	vector<double> x_0diff = x;
	x_0diff(0) += increment(0);
	vector<double> x_1diff = x;
	x_1diff(1) += increment(1);
	vector<double> x_2diff = x;
	x_2diff(2) += increment(2);
	vector<double> x_3diff = x;
	x_3diff(3) += increment(3);
	vector<double> x_4diff = x;
	x_4diff(4) += increment(4);

	num_gradients(0) = (function_value(x_0diff) - function_value(x)) / (x_0diff(0) - x(0));
	num_gradients(1) = (function_value(x_1diff) - function_value(x)) / (x_1diff(1) - x(1));
	num_gradients(2) = (function_value(x_2diff) - function_value(x)) / (x_2diff(2) - x(2));
	num_gradients(3) = (function_value(x_3diff) - function_value(x)) / (x_3diff(3) - x(3));
	num_gradients(4) = (function_value(x_4diff) - function_value(x)) / (x_4diff(4) - x(4));

	*/

	return num_gradients;
}

double T_Copula::calcStepSize(vector<double> const& x, vector<double> const& d) {

	double a = 1;
	double c1 = pow(10, -4);
	double c2 = 0.9;
	/*
	while (x(0) + a * d(0) < 0 || x(1) + a * d(1) < 0 || x(2) + a * d(2) < 0 || x(1) + a * d(1) + x(2) + a * d(2) >= 1 || x(4) + a * d(4) <= 2) {
		//std::cout <<"I bilvillkor, x_nytt = " <<  x + a * d << "\n\n";
		a = a * 0.5;

		if (a == 0) {
			break;
		}
	}

	std::cout << "steglängd efter bivillkor = " << a << "\n";



	while (function_value(x + a * d) > function_value(x) + c1 * a * inner_prod(calcGradients(x), d))
		//while (function_value(x + a * d) > function_value(x))
	{
		a = a * 0.5;
		if (a == 0) {
			break;
		}
	}

	std::cout << "steglängd efter funktionsvärdeskoll = " << a << "\n \n";

	*/
	return a;
}