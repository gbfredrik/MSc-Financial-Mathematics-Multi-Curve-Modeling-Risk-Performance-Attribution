#include "Gaussian.h"


Gaussian::Gaussian(vector<double> series) : Distribution(series) {
	time_series = series;
}


vector<double> Gaussian::create_GARCH_vec(vector<double> x) {

	vector<double> garch_vec(time_series.size() - 1);

	garch_vec(0) = pow(time_series(0), 2); // datum växer med index

	for (unsigned i = 1; i < garch_vec.size(); i++) {
		garch_vec(i) = x(0) + (1) * garch_vec(i - 1) + x(3) * time_series(i - 1);
	}

	return garch_vec;
}

double Gaussian::function_value(vector<double> x) {

	GARCH_vec = create_GARCH_vec(x);
	
	double sum = 0;

	for (unsigned i = 1; i < GARCH_vec.size() -1; i++) {
		sum = sum + log(2 * 3.14 * GARCH_vec(i)) + time_series(i + 1) / (2 * GARCH_vec(i));
	}

	sum = (1 / 2) * sum;

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
	
	return time_series;
}