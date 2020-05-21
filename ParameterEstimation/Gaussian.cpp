#include "Gaussian.h"
#include "../MathLibrary/matrixOperations.h"

#include <Eigen/Dense>

using namespace boost::numeric::ublas;

Gaussian::Gaussian(vector<double> series) : Distribution(series) {
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


void Gaussian::getSeries() {
	std::cout << "In gaussian: " << time_series << "\n";
}

//Update garch vector with new parameters
void Gaussian::update_GARCH_vec(vector<double> x) {  // datum växer med index

	for (size_t i = 0; i < time_series.size(); i++) {
		m_GARCH_vec(i + 1) = x(0) + x(1) * pow(time_series(i), 2) + x(2) * m_GARCH_vec(i);
	}
}

//Calculate function value given parameters x
double Gaussian::function_value(vector<double> x) {

	update_GARCH_vec(x);

	double sum = 0;

	for (size_t i = 0; i < m_GARCH_vec.size() - 1; ++i) {
		sum += log(2 * boost::math::constants::pi<double>() * m_GARCH_vec(i)) + pow(time_series(i) - x(3), 2) / (m_GARCH_vec(i));
	}

	sum = sum * 0.5;

	return sum;
}

vector<double> Gaussian::calcNumGradients(vector<double> x) {

	double epsilon = 2.2*pow(10, -16);
	vector<double> increment(sqrt(epsilon) * x);
	vector<double> num_gradients(4);

	
	vector<double> x_0diff = x;
	x_0diff(0) += increment(0);
	vector<double> x_1diff = x;
	x_1diff(1) += increment(1);
	vector<double> x_2diff = x;
	x_2diff(2) += increment(2);
	vector<double> x_3diff = x;
	x_3diff(3) += increment(3);

	num_gradients(0) = (function_value(x_0diff) - function_value(x)) / (x_0diff(0) - x(0));
	num_gradients(1) = (function_value(x_1diff) - function_value(x)) / (x_1diff(1) - x(1));
	num_gradients(2) = (function_value(x_2diff) - function_value(x)) / (x_2diff(2) - x(2));
	num_gradients(3) = (function_value(x_3diff) - function_value(x)) / (x_3diff(3) - x(3));


	return num_gradients;
}


// Hessian: https://v8doc.sas.com/sashtml/ormp/chap5/sect28.htm

matrix<double> Gaussian::calcNumHessian(vector<double> x) {
	
	//double epsilon = 2.2 * pow(10, -16);
	//vector<double> increment(0.00001);
	double increment = 0.0000001;
	vector<double> num_gradients(4);

	vector<double> h(4);
	for (size_t i = 0; i < h.size(); ++i) {
		h(i) = 0;
	}

	double epsilon = 4.8 * pow(10, -6);
	
	vector<double> inc(4);
	inc(0) = epsilon * (1 + abs(x(0)));
	inc(1) = epsilon * (1 + abs(x(1)));
	inc(2) = epsilon * (1 + abs(x(2)));
	inc(3) = epsilon * (1 + abs(x(3)));
	
	vector<double> hw = h;
	hw(0) += inc(0);
	vector<double> ha = h;
	ha(1) += inc(1);
	vector<double> hb = h;
	hb(2) += inc(2);
	vector<double> hmu = h;
	hmu(3) += inc(3);

	
	double dw2 = (calcGradients(x + hw)(0) - calcGradients(x - hw)(0)) / (4.0 * inc(0)) + (calcGradients(x + hw)(0) - calcGradients(x - hw)(0)) / (4.0 * inc(0));
	double dwa = (calcGradients(x + ha)(0) - calcGradients(x - ha)(0)) / (4.0 * inc(1)) + (calcGradients(x + hw)(1) - calcGradients(x - hw)(1)) / (4.0 * inc(0));
	double dwb = (calcGradients(x + hb)(0) - calcGradients(x - hb)(0)) / (4.0 * inc(2)) + (calcGradients(x + hw)(2) - calcGradients(x - hw)(2)) / (4.0 * inc(0));
	double dwmu = (calcGradients(x + hmu)(0) - calcGradients(x - hmu)(0)) / (4.0 * inc(3)) + (calcGradients(x + hw)(3) - calcGradients(x - hw)(3)) / (4.0 * inc(0));
	double dab = (calcGradients(x + hb)(1) - calcGradients(x - hb)(1)) / (4.0 * inc(2)) + (calcGradients(x + ha)(2) - calcGradients(x - ha)(2)) / (4.0 * inc(1));
	double damu = (calcGradients(x + hmu)(1) - calcGradients(x - hmu)(1)) / (4.0 * inc(3)) + (calcGradients(x + ha)(3) - calcGradients(x - ha)(3)) / (4.0 * inc(1));
	double dbmu = (calcGradients(x + hmu)(2) - calcGradients(x - hmu)(2)) / (4.0 * inc(3)) + (calcGradients(x + hb)(3) - calcGradients(x - hb)(3)) / (4.0 * inc(2));
	double da2 = (calcGradients(x + ha)(1) - calcGradients(x - ha)(1)) / (4.0 * inc(1)) + (calcGradients(x + ha)(1) - calcGradients(x - ha)(1)) / (4.0 * inc(1));
	double db2 = (calcGradients(x + hb)(2) - calcGradients(x - hb)(2)) / (4.0 * inc(2)) + (calcGradients(x + hb)(2) - calcGradients(x - hb)(2)) / (4.0 * inc(2));
	double dmu2 = (calcGradients(x + hmu)(3) - calcGradients(x - hmu)(3)) / (4.0 * inc(3)) + (calcGradients(x + hmu)(3) - calcGradients(x - hmu)(3)) / (4.0 * inc(3));

	matrix<double> Hessian(4, 4);
	Hessian(0, 0) = dw2;
	Hessian(0, 1) = dwa;
	Hessian(0, 2) = dwb;
	Hessian(0, 3) = dwmu;
	Hessian(1, 0) = dwa;
	Hessian(1, 1) = da2;
	Hessian(1, 2) = dab;
	Hessian(1, 3) = damu;
	Hessian(2, 0) = dwb;
	Hessian(2, 1) = dab;
	Hessian(2, 2) = db2;
	Hessian(2, 3) = dbmu;
	Hessian(3, 0) = dwmu;
	Hessian(3, 1) = damu;
	Hessian(3, 2) = dbmu;
	Hessian(3, 3) = dmu2;



	matrix<double> Hessian_inv = matrixOperations::matrixXdToUblas( matrixOperations::ublasToMatrixXd(Hessian).inverse() );



	return Hessian_inv;
}


//Calculate gradients, x = [omega alpha beta mu]
vector<double> Gaussian::calcGradients(vector<double> x) {

	update_GARCH_vec(x);
	
	double dw = 0;
	double da = 0;
	double db = 0;
	double dmu = 0;

	vector<double> inner_w  = derivative_w(x);
	vector<double> inner_a = derivative_a(x);
	vector<double> inner_b = derivative_b(x);

	double temp{};
	for (size_t i = 0; i < m_GARCH_vec.size()-1; i++) {
	
		temp = 0.5 * (m_GARCH_vec(i) - pow(time_series(i) - x(3), 2)) / (pow(m_GARCH_vec(i), 2));

		dw = dw + temp * inner_w(i);
		da = da + temp * inner_a(i);
		db = db + temp * inner_b(i);
		dmu = dmu - (time_series(i) - x(3)) / m_GARCH_vec(i);
	}

	vector<double> gradients(4);
	gradients(0) = dw;
	gradients(1) = da;
	gradients(2) = db;
	gradients(3) = dmu;

	return gradients;
}

//Calculate vector with derivatives of garch with respect to omega. dv/dw.
vector<double> Gaussian::derivative_w(vector<double> x) {
	vector<double> inner_dw(time_series.size());

	inner_dw(0) = 0;

	for (size_t i = 1; i < inner_dw.size(); ++i) {
		inner_dw(i) = 1 + x(2) * inner_dw(i-1);
	}

	return inner_dw;
}

//Calculate vector with derivatives of garch with respect to alpha. dv/da.
vector<double> Gaussian::derivative_a(vector<double> x) {
	vector<double> inner_da(time_series.size());

	inner_da(0) = 0;

	for (size_t i = 1; i < inner_da.size(); ++i) {
		inner_da(i) = pow(time_series(i - 1), 2) + x(2) * inner_da(i - 1);
	}

	return inner_da;
}
//Calculate vector with derivatives of garch with respect to beta. dv/db
vector<double> Gaussian::derivative_b(vector<double> x) {
	vector<double> inner_db(time_series.size());

	inner_db(0) = 0;

	for (size_t i = 1; i < inner_db.size(); ++i) {
		inner_db(i) = m_GARCH_vec(i - 1) + x(2) * inner_db(i - 1);
	}

	return inner_db;
}


double Gaussian::calcStepSize(vector<double> x, vector<double> d) {

	double a = 1;
	double c1 = pow(10, -3);
	double c2 = 0.9;

	while (x(0) + a * d(0) < 0 || x(1) + a * d(1) < 0 || x(2) + a * d(2) < 0 || x(1) + a * d(1) + x(2) + a * d(2) >= 1) {
		a = a * 0.5;

		if (a == 0) {
			break;
		}
	}
	std::cout << "alpha efter bivillkorsuppfyllning = " << a << "\n\n";



	//while (function_value(x + a * d) > function_value(x) + c1 * a * inner_prod(calcGradients(x), d))
	while (function_value(x + a * d) > function_value(x)+0.000001)
	{
		a = a * 0.5;

	}

	std::cout << "alpha efter funktionsvärdestest = " << a << "\n\n";

	return a;
}