#include "Student_t.h"


Student_t::Student_t(vector<double> series) : Distribution(series) {
	time_series = series;
	vector<double> garch_vec(time_series.size() + 1);
	m_GARCH_vec = garch_vec;

	//Calc variance for timeseries
	boost::accumulators::accumulator_set<double, boost::accumulators::stats<boost::accumulators::tag::variance> > acc;
	for_each(time_series.begin(), time_series.end(), boost::bind<void>(boost::ref(acc), _1));

	//Set variance as first element
	m_GARCH_vec(0) = boost::accumulators::variance(acc);

	std::cout << "GARCH0 = " << m_GARCH_vec(0) << "\n\n";
}


void Student_t::getSeries() {
	std::cout << "In gaussian: " << time_series << "\n";
}

void Student_t::update_GARCH_vec(vector<double> x) {
	vector<double> garch_vec(time_series.size() + 1);

	for (size_t i = 0; i < time_series.size(); i++) {
		m_GARCH_vec(i + 1) = x(0) + x(1) * pow(time_series(i), 2) + x(2) * m_GARCH_vec(i);
	}
}

double Student_t::function_value(vector<double> x) {

	update_GARCH_vec(x);

	double sum = 0;

	for (size_t i = 1; i < GARCH_vec.size() - 1; ++i) {
		sum += log(GARCH_vec(i))
			- (x(4) + 1) / 2 * log(1 + pow(time_series(i + 1), 2) / (GARCH_vec(i) * (x(4) - 2)));
	} 

	sum = sum - log(std::tgamma((x(4) + 1) / 2) / (std::sqrt(x(4) - 2) * std::tgamma(x(4) / 2)));

	sum = sum * 0.5;

	return sum;

}


vector<double> Student_t::derivative_w(vector<double> x, vector<double> GARCH_vec) {
	vector<double> inner_dw(time_series.size());

	inner_dw(0) = 0;

	for (size_t i = 1; i < inner_dw.size(); ++i) {
		inner_dw(i) = 1 + x(2) * inner_dw(i - 1);
	}

	return inner_dw;
}

vector<double> Student_t::derivative_a(vector<double> x, vector<double> GARCH_vec) {
	vector<double> inner_da(time_series.size());

	inner_da(0) = 0;

	for (size_t i = 1; i < inner_da.size(); ++i) {
		inner_da(i) = pow(time_series(i - 1), 2) * 252 + x(2) * inner_da(i - 1);
	}

	return inner_da;
}

vector<double> Student_t::derivative_b(vector<double> x, vector<double> GARCH_vec) {
	vector<double> inner_db(time_series.size());

	inner_db(0) = 0;

	for (size_t i = 1; i < inner_db.size(); ++i) {
		inner_db(i) = GARCH_vec(i - 1) + x(2) * inner_db(i - 1);
	}

	return inner_db;
}


vector<double> Student_t::calcGradients(vector<double> x) {

	GARCH_vec = create_GARCH_vec(x);

	double dw = 0;
	double da = 0;
	double db = 0;
	double dmu = 0;

	vector<double> inner_w = derivative_w(x, GARCH_vec);
	vector<double> inner_a = derivative_a(x, GARCH_vec);
	vector<double> inner_b = derivative_b(x, GARCH_vec);

	double temp{};
	for (size_t i = 0; i < GARCH_vec.size() - 1; i++) {

		temp = 0.5 * (GARCH_vec(i) - pow(time_series(i) - x(3) / 252, 2) * 252) / (pow(GARCH_vec(i), 2));

		dmu = dmu - 2 * (time_series(i) - x(3) / 252) / GARCH_vec(i);

		dw = dw + temp * inner_w(i);
		da = da + temp * inner_a(i);
		db = db + temp * inner_b(i);
	}

	vector<double> gradients(4);
	gradients(0) = dw;
	gradients(1) = da;
	gradients(2) = db;
	gradients(3) = dmu;


	/*
	double increment{ 0.00001 };
	vector<double> num_gradients(4);
	vector<double> x_0diff(x);
	x_0diff(0) += increment;
	vector<double> x_1diff(x);
	x_1diff(1) += increment;
	vector<double> x_2diff(x);
	x_2diff(2) += increment;
	vector<double> x_3diff(x);
	x_3diff(3) += increment;


	num_gradients(0) = (function_value(x_0diff) - function_value(x)) / increment;
	num_gradients(1) = (function_value(x_1diff) - function_value(x)) / increment;
	num_gradients(2) = (function_value(x_2diff) - function_value(x)) / increment;
	num_gradients(3) = (function_value(x_3diff) - function_value(x)) / increment;

	*/


	return gradients;
}

double Student_t::calcStepSize(vector<double> x, vector<double> d) {

	double a = 1;
	double c1 = pow(10, -4);
	double c2 = 0.9;

	while (x(0) + a * d(0) < 0 || x(1) + a * d(1) < 0 || x(2) + a * d(2) < 0 || x(1) + a * d(1) + x(2) + a * d(2) >= 1) {

		std::cout << "old parameters = " << x << "\n";
		std::cout << "new parameters = " << x(0) + a * d(0) << ", " << x(1) + a * d(1) << ", " << x(2) + a * d(2) << x(3) + a * d(3) << "\n";
		a = a * 0.5;
		std::cout << "steglängd = " << a << "\n";
	}

	std::cout << " old parameters = " << x << "\n";
	std::cout << "final new parameters = " << x(0) + a * d(0) << ", " << x(1) + a * d(1) << ", " << x(2) + a * d(2) << x(3) + a * d(3) << "\n";
	std::cout << "steglängd = " << a << "\n";



	while (function_value(x + a * d) > function_value(x) + c1 * a * inner_prod(calcGradients(x), d))
		//while (function_value(x + a * d) > function_value(x))
	{
		std::cout << "new f : " << function_value(x + a * d) << "\n";
		std::cout << "old f : " << function_value(x) << "\n";
		std::cout << "steglängd = " << a << "\n";
		std::cout << "new parameters2 = " << x(0) + a * d(0) << ", " << x(1) + a * d(1) << ", " << x(2) + a * d(2) << x(3) + a * d(3) << "\n";
		a = a * 0.5;
	}

	std::cout << "new f : " << function_value(x + a * d) << "\n";
	std::cout << "old f : " << function_value(x) << "\n";
	std::cout << "steglängd = " << a << "\n \n";

	return a;
}