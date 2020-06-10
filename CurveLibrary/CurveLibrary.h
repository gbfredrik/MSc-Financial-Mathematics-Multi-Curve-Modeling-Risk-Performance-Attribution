// CurveLibrary.h : Declares all DLL functions
#pragma once
#define EXPORT comment(linker, "/EXPORT:" __FUNCTION__ "=" __FUNCDNAME__)

#include <boost/numeric/ublas/matrix.hpp>

LONG __stdcall testXL(int x, int& y);

BOOL __stdcall run_all_multiXL(
	int const eigen_algorithm,
	bool const eval_eigen,
	double* return_norm_errors
);

struct CurveCollection {
	std::string filename;
	boost::numeric::ublas::matrix<double> m_A;
	boost::numeric::ublas::matrix<double> m_diff;
	int k;
	boost::numeric::ublas::matrix<double> m_E;
	boost::numeric::ublas::vector<double> v_Lambda;
	boost::numeric::ublas::matrix<double> m_delta_xi;
	boost::numeric::ublas::vector<double> v_norm_errors;
};

void placeholder_eigen(
	CurveCollection& curve_collection,
	int const eigen_algorithm,
	bool const eval_eigen,
	int const k,
    bool const save = false
);
