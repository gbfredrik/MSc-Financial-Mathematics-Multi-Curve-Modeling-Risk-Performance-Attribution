#include "sample_handler.h"
//#include "pch.h"

#include <iostream>
#include <fstream>

#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/matrix.hpp>



void print_hello() {
	std::cout << "hello" << std::endl;

}

boost::numeric::ublas::matrix<double> read_csv_matrix() {
	boost::numeric::ublas::matrix<double> m;
	std::ifstream inf("../../ExampleData/25x10950.txt");
	if (!inf) {
		std::cout << "Failed to open file" << std::endl;
	}

	if (!(inf >> m)) {
		std::cout << "Failed to write contents to matrix" << std::endl;
	}

	inf.close();
	return m;
}
