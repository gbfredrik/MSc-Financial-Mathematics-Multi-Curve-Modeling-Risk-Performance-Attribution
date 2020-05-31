#include "pch.h"
#include "sample_handler.h"

#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/matrix.hpp>

#include <fstream>
#include <iostream>
#include <string>

using namespace boost::numeric;

ublas::matrix<double> read_txt_matrix(std::string const& file_name) {
	ublas::matrix<double> m;
	std::ifstream inf("./MSc Git/MScCurveModeling/Data/" + file_name);

	if (!inf) {
		std::cout << "Failed to open file" << std::endl;
	}

	if (!(inf >> m)) {
		std::cout << "Failed to write contents to matrix" << std::endl;
	}

	inf.close();
	return m;
}


bool write_txt_matrix(ublas::matrix<double> const& m, std::string const& file_name) {
	std::ofstream outf;
	outf.open("./MSc Git/MScCurveModeling/Data/m_" + file_name); // Appending "m_" to prevent accidental overwriting of important files

	if (!(outf << m)) { // If failed write
		outf.close();
		return 0;
	}

	outf.close();
	return 1;
}

bool write_txt_vector(ublas::vector<double> const& m, std::string const& file_name) {
	std::ofstream outf;
	outf.open("./MSc Git/MScCurveModeling/Data/v_" + file_name); // Appending "m_" to prevent accidental overwriting of important files

	if (!(outf << m)) { // If failed write
		outf.close();
		return 0;
	}

	outf.close();
	return 1;
}

ublas::matrix<double> read_csv_matrix(std::string const& file_name) {
	std::ifstream inf;
	inf.open("./MSc Git/MScCurveModeling/Data/" + file_name);

	int rows{ 0 };
	int cols{ 0 };
	char delim;
	inf >> rows >> delim >> cols >> delim;
	ublas::matrix<double> m(rows, cols);
	
	std::string line;
	std::string val;
	for (int i = 0; i < rows; ++i) {
		getline(inf, line);
		std::stringstream s(line);

		for (int j = 0; j < cols; ++j) {
			getline(s, val, ';');
			m(i, j) = std::stod(val);
			//inf >> m(i, j);
			//inf >> delim;
		}
	}

	inf.close();
	return m;
}

ublas::vector<double> read_csv_vector(std::string const& file_name) {
	std::ifstream inf;
	inf.open("./MSc Git/MScCurveModeling/Data/" + file_name);

	int length{ 0 };
	char delim;
	inf >> length;
	ublas::vector<double> v(length);

	for (int i = 0; i < length; ++i) {
		inf >> v(i);
		inf >> delim;
	}

	inf.close();
	return v;
}

bool write_csv_matrix(ublas::matrix<double> const& m, std::string const& file_name) {
	std::ofstream outf;
	outf.open("./MSc Git/MScCurveModeling/Data/m_" + file_name); // Appending "m_" to prevent accidental overwriting of important files

	outf << m.size1() << ";" << m.size2() << std::endl;

	for (size_t r = 0, rows = m.size1(); r < rows; ++r) {
		for (size_t c = 0, cols = m.size2(); c < cols; ++c) {
			outf << m(r, c) << ((c != cols - 1) ? ";" : "");
		}
		outf << std::endl;
	}

	return true;
}

bool write_csv_vector(ublas::vector<double> const& m, std::string const& file_name) {
	std::ofstream outf;
	outf.open("./MSc Git/MScCurveModeling/Data/v_" + file_name); // Appending "m_" to prevent accidental overwriting of important files

	outf << m.size() << std::endl;

	for (size_t i = 0, length = m.size(); i < length; ++i) {
		outf << m(i) << ((i != length - 1) ? ";" : "");
	}
	outf << std::endl;

	return true;
}



bool placeholder_ir_measurement_multi(
	ublas::matrix<double>& m_rf, 
	ublas::matrix<double>& m_tenor
) {
	m_rf = read_csv_matrix("fHist.csv");
	m_tenor = read_csv_matrix("piHist.csv");

	return (m_rf.size1() > 0) && (m_tenor.size2() > 0);
}
