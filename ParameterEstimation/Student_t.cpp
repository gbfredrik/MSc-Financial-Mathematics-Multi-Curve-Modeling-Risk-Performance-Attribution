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

}


void Student_t::getSeries() {
	std::cout << time_series << "\n";
}

void Student_t::update_GARCH_vec(vector<double> const& x) {
	vector<double> garch_vec(time_series.size() + 1);

	for (size_t i = 0; i < time_series.size(); i++) {
		m_GARCH_vec(i + 1) = x(0) + x(1) * pow(time_series(i), 2) + x(2) * m_GARCH_vec(i);
	}
}

double Student_t::function_value(vector<double> const& x) {

	update_GARCH_vec(x);

	double sum = 0;

	for (size_t i = 0; i < m_GARCH_vec.size() - 1; ++i) {
		//sum += log(m_GARCH_vec(i))
		//	+ (x(4) + 1) / 2 * log(1 + pow(time_series(i)-x(3), 2) / (m_GARCH_vec(i) * (x(4) - 2)));
	
		sum += log(m_GARCH_vec(i)) + (x(4) + 1) * log(1 + pow(time_series(i) - x(3), 2) / (x(4) * m_GARCH_vec(i)));
	} 

	//sum = 0.5*sum - log(std::tgamma((x(4) + 1) / 2) / (std::sqrt(x(4) - 2) * std::tgamma(x(4) / 2)));
	sum = 0.5 * sum - std::lgamma((x(4) + 1) / 2) + lgamma(x(4) / 2) + 1 / 2 * log(boost::math::constants::pi<double>())
		+ 1 / 2 * log(x(4));

	return sum;

}


vector<double> Student_t::derivative_w(vector<double> const& x) {
	vector<double> inner_dw(time_series.size());

	inner_dw(0) = 0;

	for (size_t i = 1; i < inner_dw.size(); ++i) {
		inner_dw(i) = 1 + x(2) * inner_dw(i - 1);
	}

	return inner_dw;
}

vector<double> Student_t::derivative_a(vector<double> const& x) {
	vector<double> inner_da(time_series.size());

	inner_da(0) = 0;

	for (size_t i = 1; i < inner_da.size(); ++i) {
		inner_da(i) = pow(time_series(i - 1), 2) + x(2) * inner_da(i - 1);
	}

	return inner_da;
}

vector<double> Student_t::derivative_b(vector<double> const& x) {
	vector<double> inner_db(time_series.size());

	inner_db(0) = 0;

	for (size_t i = 1; i < inner_db.size(); ++i) {
		inner_db(i) = m_GARCH_vec(i - 1) + x(2) * inner_db(i - 1);
	}

	return inner_db;
}


vector<double> Student_t::calcGradients(vector<double> const& x) {
	
	update_GARCH_vec(x);

	double dw = 0;
	double da = 0;
	double db = 0;
	double dmu = 0;
	double ddf = 0;

	vector<double> inner_w = derivative_w(x);
	vector<double> inner_a = derivative_a(x);
	vector<double> inner_b = derivative_b(x);

	double temp{};
	for (size_t i = 0; i < m_GARCH_vec.size() - 1; i++) {

		//temp = 0.5 * (1/m_GARCH_vec(i) + ((x(4)+1)/2)*((m_GARCH_vec(i)*(x(4)-2))
		//								/(m_GARCH_vec(i)*(x(4)-2) + pow(time_series(i)-x(3),2)))
		//								*(-pow(time_series(i)-x(3),2))/(pow(m_GARCH_vec(i),2)*(x(4)-2)));

		temp = 0.5 * (1 / m_GARCH_vec(i) + (x(4) + 1) * ((m_GARCH_vec(i) * x(4))
			/ (m_GARCH_vec(i) * x(4) + pow(time_series(i) - x(3), 2)))
			* (-pow(time_series(i) - x(3), 2)) / (pow(m_GARCH_vec(i), 2) * x(4)));

		//dmu = dmu - 0.5 * (x(4) + 1) * (time_series(i) - x(3)) / (m_GARCH_vec(i) * (x(4) - 2) + pow(time_series(i) - x(3), 2));
		dmu = dmu - (x(4) + 1) * (time_series(i) - x(3)) / (m_GARCH_vec(i) * x(4) + pow(time_series(i) - x(3), 2));

		dw = dw + temp * inner_w(i);
		da = da + temp * inner_a(i);
		db = db + temp * inner_b(i);
	}



	vector<double> gradients(x.size());
	gradients(0) = dw;
	gradients(1) = da;
	gradients(2) = db;
	gradients(3) = dmu;
	gradients(4) = calcNumGradients(x)(4);

	return gradients;
}

vector<double> Student_t::calcNumGradients(vector<double> const& x) {

	double epsilon = 2.2 * pow(10, -16);
	vector<double> increment(sqrt(epsilon) * x);
	vector<double> num_gradients(x.size());


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


	return num_gradients;
}

double Student_t::calcStepSize(vector<double> const& x, vector<double> const& d) {

	double a = 1;
	double c1 = pow(10, -4);
	double c2 = 0.9;

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

	return a;
}