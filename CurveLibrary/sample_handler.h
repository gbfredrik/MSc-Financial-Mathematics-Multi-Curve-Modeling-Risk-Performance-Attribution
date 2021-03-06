#pragma once

#include <boost/numeric/ublas/matrix.hpp>

#include <string>

boost::numeric::ublas::matrix<double> read_txt_matrix(std::string const& file_name);
bool write_txt_matrix(
    boost::numeric::ublas::matrix<double> const& m, 
    std::string const& file_name
);
bool write_txt_vector(
    boost::numeric::ublas::vector<double> const& v, 
    std::string const& file_name
);

boost::numeric::ublas::matrix<double> read_csv_matrix(std::string const& file_name);
boost::numeric::ublas::vector<double> read_csv_vector(std::string const& file_name);
bool write_csv_matrix(
    boost::numeric::ublas::matrix<double> const& m, 
    std::string const& file_name
);
bool write_csv_vector(
    boost::numeric::ublas::vector<double> const& v, 
    std::string const& file_name
);
