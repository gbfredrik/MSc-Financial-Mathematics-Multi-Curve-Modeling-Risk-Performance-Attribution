// RawInterpolTest.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/range/counting_range.hpp>
#include "../../CurveLibrary/sample_handler.h"

using namespace boost::numeric;
//using namespace boost::range;

void price_rawInterpol(double const spot_price, int const maturity_instrument, ublas::vector<int> maturity_data, ublas::vector<double> mid_price_data);
void curve_rawInterpol(ublas::vector<int> maturity_data, ublas::vector<double> mid_price_input);
ublas::vector<double> basis_calc(size_t const length, ublas::vector<double> maturity_input, ublas::vector<double> mid_price_input);
ublas::vector<double> delta_T_calc(size_t const length, ublas::vector<double> maturity_input);

int main()
{	
	double spot_price = 1.1118; //Excel
	int maturity_instrument = 642; //Excel
	ublas::vector<int> maturity_data;
	ublas::vector<double> mid_price_data;

	maturity_data = read_csv_vector<int>("deltaT_EUR.csv");
	std::cout << "deltaTdata: " << maturity_data << std::endl;
	mid_price_data = read_csv_vector<double>("basisFX_EUR.csv");
	std::cout << "basis: " << mid_price_data << std::endl;

	price_rawInterpol(spot_price, maturity_instrument, maturity_data, mid_price_data);
	
	std::cout << "Curve" << std::endl;
	std::cout << "-----" << std::endl;
	curve_rawInterpol(maturity_data, mid_price_data);

}


//Calculates forward price
void price_rawInterpol(double const spot_price, int const maturity_instrument, ublas::vector<int> maturity_data, ublas::vector<double> mid_price_data) {
	double basis_m = 0;
	double forward_price = 0;
	size_t maturity_pos = 0;
	size_t length = mid_price_data.size();
	ublas::vector<double> basis(length);
	ublas::vector<double> delta_T(length);

	for (size_t i = 0; i < length; ++i) {

		if (maturity_instrument < maturity_data(i)) {
			double T = maturity_data(i) - maturity_data(i - 1) - (maturity_data(i) - maturity_instrument);
			basis_m = (mid_price_data(i) - mid_price_data(i - 1)) / ((maturity_data(i) - maturity_data(i - 1)) / 365.0) * T / 365.0; //Hardcoded
			maturity_pos--;
			break;
		}else if (maturity_instrument <= maturity_data(i)) {
			basis_m = 0;
			break;
		}else if (maturity_instrument > maturity_data(length-1)) {
			std::cout << "Cannot price contract. Maturity of instrument is larger than the curve. " << std::endl;
		}
		maturity_pos++;
	}
	basis = basis_calc(maturity_pos + 1, maturity_data/365.0, mid_price_data);
	delta_T = delta_T_calc(maturity_pos + 1, maturity_data/365.0);

	forward_price = spot_price+inner_prod(basis, delta_T)+basis_m;

	std::cout << "-----" <<std::endl;
	std::cout << "forward_price: " << forward_price << std::endl;
	std::cout << "-----" <<std::endl;

}

//Calculates raw interpolation curve. 
void curve_rawInterpol(ublas::vector<int> maturity_data, ublas::vector<double> mid_price_data) {
	size_t length = mid_price_data.size();
	ublas::vector<double> delta_T(length);
	delta_T = delta_T_calc(length, maturity_data);
	size_t curve_length = sum(delta_T);
	ublas::vector<double> curve(curve_length);
	ublas::vector<double> basis(length);
	size_t pos = 0;

	basis = basis_calc(length, maturity_data/365.0, mid_price_data); //Hardcoded

	for (size_t j = 0; j < length-1; ++j) { //# of maturities
		for (size_t k = 0; k < delta_T(j); ++k) { //# of days between every maturity
			curve(pos) = basis(j);
			pos++;
		}
	}
	std::cout << "curve: " << curve << std::endl;
	//write_csv_vector(curve, "RawInterpol.csv");
}

//Calculates basis point vector.
ublas::vector<double> basis_calc(size_t length, ublas::vector<double> maturity_data, ublas::vector<double> mid_price_data) {
	ublas::vector<double> basis(length);
	ublas::vector<double> delta_T(length);
	delta_T = delta_T_calc(length, maturity_data);
	ublas::vector_range<ublas::vector<double>> mid_price0(mid_price_data, ublas::range(0, length-1));
	ublas::vector_range<ublas::vector<double>> mid_price1(mid_price_data, ublas::range(1, length));

	basis = element_div((mid_price1 - mid_price0), delta_T);

	return basis;
}

//Calculates difference in maturity, delta_T
ublas::vector<double> delta_T_calc(size_t length, ublas::vector<double> maturity_input) {
	ublas::vector<double> delta_T(length);
	ublas::vector_range<ublas::vector<double>> maturity0(maturity_input, ublas::range(0, length-1));
	ublas::vector_range<ublas::vector<double>> maturity1(maturity_input, ublas::range(1, length));
	delta_T = (maturity1 - maturity0);
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

