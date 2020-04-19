/**
	This file contains a main function for testing of 
	the implementations in "/CurveLibrary/sample_handler.h".
*/

#include <iostream>

#include ".../../CurveLibrary/sample_handler.h"
#include <boost/numeric/ublas/matrix.hpp>

int main() {
	boost::numeric::ublas::matrix<double> m{read_txt_matrix("25x10950.txt")};

	std::cout << "size1: " << m.size1() << std::endl;
	std::cout << "size2: " << m.size2() << std::endl;
	
	std::cout << write_txt_matrix(m, "testmatrix.txt");
}
