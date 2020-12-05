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
	int const curve_length,
    char const* _k_risk_factors,
    double* return_approx_errors
);

BOOL __stdcall run_all_multi_risk_measuresXL(
    char const* _file_path,
    char const* _file_names_curve_sets,
    char const* _file_names_risk_factors,
    bool const backtest,
    double* statistic_m_hypothesis,
    bool pass_hypothesis,
    double* statistic_christoffersen,
    bool pass_christoffersen,
    double* statistic_p_acerbi,
    bool pass_acerbi
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
    bool const save = true
);
