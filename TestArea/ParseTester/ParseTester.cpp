// ParseTester.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>

#include ".../../CurveLibrary/sample_handler.h"
#include <boost/numeric/ublas/matrix.hpp>

int main() {
	boost::numeric::ublas::matrix<double> m{read_csv_matrix()};

	std::cout << "size1: " << m.size1() << std::endl;
	std::cout << "size2: " << m.size2() << std::endl;
}