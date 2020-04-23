#pragma once

#include <boost/numeric/ublas/matrix.hpp>

#include <string>

boost::numeric::ublas::matrix<double> read_txt_matrix(std::string const&);
bool write_txt_matrix(boost::numeric::ublas::matrix<double> const&, std::string const&);

bool placeholder_ir_measurement_multi(boost::numeric::ublas::matrix<double>&,
									  boost::numeric::ublas::matrix<double>&);