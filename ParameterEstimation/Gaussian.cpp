#include "Gaussian.h"


Gaussian::Gaussian(vector<double> series) : Distribution(series) {
	time_series = series;
}


void Gaussian::getSeries() {
	std::cout << "In gaussian: " << time_series << "\n";
}

vector<double> Gaussian::create_GARCH_vec(vector<double> x) {
	vector<double> garch_vec(time_series.size());

	garch_vec(0) = pow(time_series(0), 2); // datum växer med index


	for (size_t i = 1; i < garch_vec.size(); i++) {
		garch_vec(i) = x(0) + x(1) * pow(time_series(i), 2) + x(2) * garch_vec(i - 1);
	}

	return garch_vec;
}

double Gaussian::function_value(vector<double> x) {

	GARCH_vec = create_GARCH_vec(x);
	//std::cout << "\n Garchvec = " << GARCH_vec << std::endl;
	
	double sum = 0;

	for (size_t i = 1; i < GARCH_vec.size() - 1; ++i) {
		sum += log(2 * boost::math::constants::pi<double>() * GARCH_vec(i)) + pow(time_series(i + 1),2) / (GARCH_vec(i));

	}
	//sum *= 0.5;
	sum = sum * 0.5;

	return sum;

}


vector<double> Gaussian::calcGradients(vector<double> x) {

	GARCH_vec = create_GARCH_vec(x);
	/*
	double dw = 0;
	double da = 0;
	double db = 0;

	double temp{};
	for (size_t i = 1; i < GARCH_vec.size()-1; i++) {

		temp = 0.5*(GARCH_vec(i) - pow(time_series(i+1), 2)) / (pow(GARCH_vec(i), 2));

		dw = dw + temp;
		da = da + temp * pow(time_series(i), 2);
		db = db + temp * GARCH_vec(i - 1);
	}
	*/
	double increment{ 0.00001 };
	vector<double> gradients(3);
	vector<double> x_0diff(x);
	x_0diff(0) += increment;
	vector<double> x_1diff(x);
	x_1diff(1) += increment;
	vector<double> x_2diff(x);
	x_2diff(2) += increment;


	gradients(0) = (function_value(x_0diff) - function_value(x)) / increment;
	gradients(1) = (function_value(x_1diff) - function_value(x)) / increment;
	gradients(2) = (function_value(x_2diff) - function_value(x)) / increment;
	return gradients;
}

double Gaussian::calcStepSize(vector<double> x, vector<double> d) {

	double a = 1;
	double c1 = pow(10, -4);
	double c2 = 0.9;

	while (x(0) + a * d(0) < 0 || x(1) + a * d(1) < 0 || x(2) + a * d(2) < 0 || x(1) + a * d(1) + x(2) + a * d(2) >= 1) {

		std::cout << "old parameters = " << x << "\n";
		std::cout << "new parameters = " << x(0) + a * d(0) << ", " << x(1) + a * d(1) << ", " << x(2) + a * d(2) << "\n";
		a = a * 0.5;
		std::cout << "steglängd = " << a << "\n";
	}

	std::cout << " old parameters = " << x << "\n";
	std::cout << "final new parameters = " << x(0) + a * d(0) << ", " << x(1) + a * d(1) << ", " << x(2) + a * d(2) << "\n";
	std::cout << "steglängd = " << a << "\n";



	//while (dist->function_value(x + a * d) > dist->function_value(x) + c1 * a * inner_prod(dist->calcGradients(x), d)
	while (function_value(x + a * d) > function_value(x))
	{
		std::cout << "new f : " << function_value(x + a * d) << "\n";
		std::cout << "old f : " << function_value(x) << "\n";
		std::cout << "steglängd = " << a << "\n";
		std::cout << "new parameters2 = " << x(0) + a * d(0) << ", " << x(1) + a * d(1) << ", " << x(2) + a * d(2) << "\n";
		a = a * 0.5;
	}

	std::cout << "new f : " << function_value(x + a * d) << "\n";
	std::cout << "old f : " << function_value(x) << "\n";
	std::cout << "steglängd = " << a << "\n \n";

	return a;
}