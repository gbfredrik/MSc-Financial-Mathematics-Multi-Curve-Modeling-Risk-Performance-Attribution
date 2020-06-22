#include "pch.h"

#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/matrix.hpp>

#include <fstream>
#include <iostream>
#include <string>

using namespace boost::numeric;

inline ublas::matrix<double> read_txt_matrix(std::string const& file_name) {
	ublas::matrix<double> m;
	std::ifstream inf("./MSc Git/Data/" + file_name);

    if (!inf) {
        std::cout << "Failed to open file" << std::endl;
    }

    if (!(inf >> m)) {
        std::cout << "Failed to read contents to matrix" << std::endl;
    }

    inf.close();
    return m;
}

inline bool write_txt_matrix(ublas::matrix<double> const& m, std::string const& file_name) {
	std::ofstream outf;
	outf.open("./MSc Git/Data/m_" + file_name); // Appending "m_" to prevent accidental overwriting of important files

    if (!(outf << m)) { // If failed write
        outf.close();
        return 0;
    }

    outf.close();
    return 1;
}

inline bool write_txt_vector(ublas::vector<double> const& v, std::string const& file_name) {
	std::ofstream outf;
	outf.open("./MSc Git/Data/v_" + file_name); // Appending "v_" to prevent accidental overwriting of important files

    if (!(outf << v)) { // If failed write
        outf.close();
        return 0;
    }

    outf.close();
    return 1;
}

template<typename T>
ublas::matrix<T> read_csv_matrix(std::string const& file_name) {
	std::ifstream inf;
	inf.open("X:/Examensarbete/Data/" + file_name);

	size_t rows{ 0 };
    size_t cols{ 0 };
	char delim;
	inf >> rows >> delim >> cols;// >> delim;
	if (typeid(T) == typeid(double)) {
		inf >> delim;
	}
	ublas::matrix<T> m(rows, cols);

    std::string line{};
    std::string val{};
    for (size_t i{ 0 }; i < rows; ++i) {
		getline(inf, line);
		std::stringstream s(line);

		if (typeid(T) == typeid(int)) {
            for (size_t j{ 0 }; j < cols; ++j) {
				inf >> m(i, j) >> delim;
			}
		} else if (typeid(T) == typeid(double)) {
            for (size_t j{ 0 }; j < cols; ++j) {
				getline(s, val, ';');
				m(i, j) = std::stod(val);
			}
		}
	}

	inf.close();
	return m;
}

template<typename T>
ublas::vector<T> read_csv_vector(std::string const& file_name) {
	std::ifstream inf;
	inf.open("X:/Examensarbete/Data/" + file_name);

    if (!inf) {
        std::cout << "Cannot open file.\n\n";
    }

	size_t length{};
	char delim{};
	inf >> length;

	if (typeid(T) == typeid(double)) {
		inf >> delim;
	}
	ublas::vector<T> v(length);

    for (size_t i{ 0 }; i < length; ++i) {
		inf >> v(i) >> delim;
	}

    inf.close();
    return v;
}

inline bool write_csv_matrix(ublas::matrix<double> const& m, std::string const& file_name) {
	std::ofstream outf;
	outf.open("X:/Examensarbete/Data/m_" + file_name); // Appending "m_" to prevent accidental overwriting of important files

    outf << m.size1() << ";" << m.size2() << std::endl;

	for (size_t r{ 0 }, rows{ m.size1() }; r < rows; ++r) {
		for (size_t c{ 0 }, cols{ m.size2() }; c < cols; ++c) {
			outf << m(r, c) << ((c != cols - 1) ? ";" : "");
		}
		outf << std::endl;
	}

    return true;
}

inline bool write_csv_vector(ublas::vector<double> const& v, std::string const& file_name) {
	std::ofstream outf;
	outf.open("X:/Examensarbete/Data/v_" + file_name); // Appending "v_" to prevent accidental overwriting of important files

    if (!outf) {
        std::cout << "Failed\n\n";
    }

    outf << v.size() << std::endl;

    for (size_t i{ 0 }, length{ v.size() }; i < length; ++i) {
        outf << v(i) << ((i != length - 1) ? ";" : "");
    }
    outf << std::endl;

    return true;
}
