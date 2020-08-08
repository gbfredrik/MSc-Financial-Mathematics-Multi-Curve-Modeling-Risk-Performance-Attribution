// CurveLibrary.h : Declares all DLL functions
#pragma once
#define EXPORT comment(linker, "/EXPORT:" __FUNCTION__ "=" __FUNCDNAME__)

#include "CurveCollection.h"

#include <boost/numeric/ublas/matrix.hpp>

LONG __stdcall testXL(int x, int& y);

BOOL __stdcall generate_multi_risk_factorsXL(
    char const* _file_path,
    char const* _file_names,
	int const eigen_algorithm,
	bool const eval_eigen,
	double* return_norm_errors,
	int const curve_length
);

BOOL __stdcall run_all_fxXL(
	int const eigen_algorithm,
	bool const eval_eigen,
	double* return_norm_errors,
	int const curve_length
);

void placeholder_eigen(
	CurveCollection& curve_collection,
	int const eigen_algorithm,
	bool const eval_eigen,
    bool const save = false
);
