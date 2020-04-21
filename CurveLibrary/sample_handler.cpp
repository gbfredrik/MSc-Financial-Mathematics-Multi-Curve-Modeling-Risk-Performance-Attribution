#include "sample_handler.h"

#include "pch.h"

#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/matrix.hpp>

#include <fstream>
#include <iostream>
#include <string>




boost::numeric::ublas::matrix<double> read_txt_matrix(std::string const& file_name) {
	boost::numeric::ublas::matrix<double> m;
	std::ifstream inf("../ExampleData/" + file_name);
	if (!inf) {
		std::cout << "Failed to open file" << std::endl;
	}

	if (!(inf >> m)) {
		std::cout << "Failed to write contents to matrix" << std::endl;
	}

	inf.close();
	return m;
}


bool write_txt_matrix(boost::numeric::ublas::matrix<double> const& m, std::string const& file_name) {
	std::ofstream outf;
	outf.open("m_" + file_name); // Appending "m_" to prevent accidental overwriting of important files

	if (!(outf << m)) { // If failed write
		outf.close();
		return false;
	}

	outf.close();
	return true;
}

bool placeholder_ir_measurement_multi(boost::numeric::ublas::matrix<double>& m_rf, 
									  boost::numeric::ublas::matrix<double>& m_tenor) {
	m_rf =  read_txt_matrix("25x10950.txt");
	m_tenor = read_txt_matrix("25x10950.txt");

	bool status{ m_rf.size1() > 0 && m_tenor.size1() > 0 };
	return status;
}
