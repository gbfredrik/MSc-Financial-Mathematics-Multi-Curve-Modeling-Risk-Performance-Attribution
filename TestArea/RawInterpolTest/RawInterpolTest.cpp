// RawInterpolTest.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/range/counting_range.hpp>

using namespace boost::numeric::ublas;
//using namespace boost::range;

void price_rawInterpol(double const spot_price, double const maturity_instrument, size_t const length, vector<double> maturity_input, vector<double> mid_price_input);
void curve_rawInterpol(size_t length, vector<double> maturity_input, vector<double> mid_price_input);
vector<double> basis_calc(size_t const length, vector<double> maturity_input, vector<double> mid_price_input);
vector<double> delta_T_calc(size_t const length, vector<double> maturity_input);

int main()
{	
	double spot_price = 9.8725;
	size_t length = 10; //Calculate
	double maturity_instrument = 10; //Calculate
	vector<double> maturity_input(length + 1);
	vector<double> mid_price_input(length + 1);

	maturity_input(0) = 0;
	maturity_input(1) = 1;
	maturity_input(2) = 5;
	maturity_input(3) = 6;
	maturity_input(4) = 12;
	maturity_input(5) = 19;
	maturity_input(6) = 27;
	maturity_input(7) = 36;
	maturity_input(8) = 68;
	maturity_input(9) = 97;
	maturity_input(10) = 128;

	mid_price_input(0) = 0; 
	mid_price_input(1) = -0.59;
	mid_price_input(2) = -3.05;
	mid_price_input(3) = -0.69;
	mid_price_input(4) = -5.38;
	mid_price_input(5) = -13.24;
	mid_price_input(6) = -21.665;
	mid_price_input(7) = -31.76;
	mid_price_input(8) = -59.95;
	mid_price_input(9) = -91.98;
	mid_price_input(10) = -116.765;
	mid_price_input = mid_price_input / 1000;

	price_rawInterpol(spot_price, maturity_instrument, length, maturity_input, mid_price_input);
	curve_rawInterpol(length, maturity_input, mid_price_input);
}


//Calculates forward price
void price_rawInterpol(double const spot_price, double const maturity_instrument, size_t const length, vector<double> maturity_input, vector<double> mid_price_input) {
	double basis_m = 0;
	double forward_price = 0;
	size_t maturity_pos = 0;
	vector<double> basis(length);
	vector<double> delta_T(length);

	for (size_t i = 0; i < length; i++) {
		if (maturity_instrument < maturity_input(i)){
			double T = maturity_input(i) - maturity_input(i - 1) - (maturity_input(i) - maturity_instrument);
			basis_m = (mid_price_input(i) - mid_price_input(i - 1)) / (maturity_input(i) - maturity_input(i - 1))*T;
			maturity_pos--;
			break;
		}
		else if (maturity_instrument <= maturity_input(i)) {
			basis_m = 0;
			break;
		}
		maturity_pos++;
	}

	basis = basis_calc(maturity_pos, maturity_input, mid_price_input);
	delta_T = delta_T_calc(maturity_pos, maturity_input);

	forward_price = spot_price+inner_prod(basis, delta_T)+basis_m;

	std::cout << "-----" <<std::endl;
	std::cout << "forward_price: " << forward_price << std::endl;
	std::cout << "-----" <<std::endl;

}

//Calculates raw interpolation curve. 
void curve_rawInterpol(size_t length, vector<double> maturity_input, vector<double> mid_price_input) {
	vector<double> delta_T(length);
	delta_T = delta_T_calc(length, maturity_input);
	size_t curve_length = sum(delta_T);
	vector<double> curve(curve_length); 
	vector<double> basis(length);
	size_t pos = 0;

	basis = basis_calc(length, maturity_input, mid_price_input);

	for (size_t j = 0; j < length; j++) { //# of maturities
		for (size_t k = 1; k <= delta_T(j); k++) { //# of days between every maturity
			curve(pos) = basis(j);
			pos++;
		}
	}
	std::cout << "curve: " << curve << std::endl;
}

//Calculates basis point vector.
vector<double> basis_calc(size_t length, vector<double> maturity_input, vector<double> mid_price_input) {
	vector<double> basis(length);
	vector<double> delta_T(length);
	delta_T = delta_T_calc(length, maturity_input);
	vector_range<vector<double>> mid_price0(mid_price_input, range(0, length));
	vector_range<vector<double>> mid_price1(mid_price_input, range(1, length + 1));

	basis = element_div((mid_price1 - mid_price0), delta_T);

	return basis;
}

//Calculates difference in maturity, delta_T
vector<double> delta_T_calc(size_t length, vector<double> maturity_input) {
	vector<double> delta_T(length);
	vector_range<vector<double>> maturity0(maturity_input, range(0, length));
	vector_range<vector<double>> maturity1(maturity_input, range(1, length + 1));

	delta_T = maturity1 - maturity0;

	return delta_T;

}


//for (size_t i = 0; i < curve_length; i++) { //# of rows
//	for (size_t j = 1; j < length; j++) { //# of maturities
//		for (double k = 1; k <= maturity_input(j) - maturity_input(j - 1); k++) { //# of days between every maturity
//			curve(i) = mid_price_input(j);
//		}
//	}
//}
// Run program: Ctrl + F5 or Debug > Start Without Debugging menu
// Debug program: F5 or Debug > Start Debugging menu

// Tips for Getting Started: 
//   1. Use the Solution Explorer window to add/manage files
//   2. Use the Team Explorer window to connect to source control
//   3. Use the Output window to see build output and other messages
//   4. Use the Error List window to view errors
//   5. Go to Project > Add New Item to create new code files, or Project > Add Existing Item to add existing code files to the project
//   6. In the future, to open this project again, go to File > Open > Project and select the .sln file

