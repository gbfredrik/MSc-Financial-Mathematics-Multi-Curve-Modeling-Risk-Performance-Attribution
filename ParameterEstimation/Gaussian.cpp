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


	for (unsigned i = 1; i < garch_vec.size(); i++) {
		garch_vec(i) = x(0) + x(1) * time_series(i - 1) + x(2) * garch_vec(i - 1);
	}



	return garch_vec;
}

double Gaussian::function_value(vector<double> x) {

	GARCH_vec = create_GARCH_vec(x);

	
	double sum = 1;

	for (size_t i = 1; i < GARCH_vec.size() -1; ++i) {
		sum += log(2 * boost::math::constants::pi<double>() * GARCH_vec(i)) + pow(time_series(i + 1),2) / (2 * GARCH_vec(i));

	}
	sum *= 0.5;



	return sum;

}


vector<double> Gaussian::calcGradients(vector<double> x) {
	
	GARCH_vec = create_GARCH_vec(x);

	

	double dw = 0;
	double da = 0;
	double db = 0;

	
	for (unsigned i = 1; i < GARCH_vec.size(); i++) {
	
		double temp = (pow(time_series(i), 2) - GARCH_vec(i)) / (pow(GARCH_vec(i), 2));

		dw = dw + temp;
		da = da + temp * time_series(i - 1);
		db = db + temp * GARCH_vec(i - 1);

	}

	vector<double> gradients(3);
	gradients(0) = dw;
	gradients(1) = da;
	gradients(2) = db;
	return gradients;
}