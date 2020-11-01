// CurveLibrary.cpp : Defines the exported functions for the DLL.
//

#include "pch.h"
#include "framework.h"
#include "CurveLibrary.h" // Must be this order, otherwise functions are not exported!

#include "CurveCollection.h"
#include "sample_handler.h"
#include "../MathLibrary/matrixOperations.h"
#include "../MathLibrary/statisticsOperations.h"
#include "../RiskFactorCalculation/FactorCalculation.h"

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>

#include <string>
#include <algorithm>
#include <vector>

#include <WinUser.h>

using namespace boost::numeric;

LONG __stdcall testXL(int x, int& y) {
#pragma EXPORT

	y += 100;
    return x * x;
}
// isBusinessDayVBA(const char* _market, int curDateExcel, int checkDateExcel, int* isBD, BSTR* holidayName)
//{std::string market(_market);


BOOL __stdcall generate_multi_risk_factorsXL(
    char const* _file_path,
    char const* _file_names,
	int const eigen_algorithm,
	bool const eval_eigen,
	double* return_norm_errors,
	int const curve_length,
	char const* _k_risk_factors
) {
#pragma EXPORT

	bool status{ 1 };
    std::string file_names{ _file_names };
    //std::string k_risk_factors{ _k_risk_factors };
    std::stringstream ss_file_names{ file_names };
	std::stringstream ss_k_risk_factors{ std::string{ _k_risk_factors } };
    std::string substr{};
	std::string substr_k{};

	int count_tenor{ std::count(file_names.begin(), file_names.end(), ';') };
	CurveCollection rf;
	std::vector<CurveCollection> tenors(count_tenor);

    // Setup and read files
	rf.file_path = _file_path;
    getline(ss_file_names, substr, ';');
    rf.file_name = substr;

	getline(ss_k_risk_factors, substr_k, ';');
	rf.k = std::stoi(substr_k);
	
	for (CurveCollection& cc : tenors) {
        cc.file_path = _file_path;
        getline(ss_file_names, substr, ';');
		cc.file_name = substr;

		getline(ss_k_risk_factors, substr_k, ';');
        cc.k = std::stoi(substr_k);
	}

	try {
		// Read
		rf.m_A = read_csv_matrix(rf.file_path + rf.file_name);
		for (CurveCollection& cc : tenors) {
			cc.m_A = read_csv_matrix(cc.file_path + cc.file_name);
		}
	} catch (std::exception const&) {
		return -1;
	}

	rf.m_A_trunc = rf.m_A;
	rf.m_A_trunc.resize(rf.m_A_trunc.size1(), curve_length);
	for (CurveCollection& cc : tenors) {
		cc.m_A_trunc = cc.m_A;
		cc.m_A_trunc.resize(cc.m_A_trunc.size1(), curve_length);
	}


	// Risk factor calculation
	placeholder_eigen(rf, eigen_algorithm, eval_eigen);
    for (CurveCollection& cc : tenors) {
        placeholder_eigen(cc, eigen_algorithm, eval_eigen);
    }

	return status;
}

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
) {
#pragma EXPORT

    bool status{ 1 };

    // Parameter estimation


    // Risk measurement
    // for (int date : oos_dates) { 
        // Run simulation 
        // Measure risk 
    //}

    // Backtesting
    if (backtest) {
        NULL;
    }

    // Save results


    return status;
}




BOOL __stdcall run_all_fxXL(
	int const eigen_algorithm,
	bool const eval_eigen,
	double* return_norm_errors,
	int const curve_length
) {
#pragma EXPORT

	bool status{ 1 };
	int count_rf{ 2 };
	int count_tenor{ 2 };

	std::vector<CurveCollection> rf(count_rf);
	std::vector<CurveCollection> tenors(count_tenor);

	rf.at(0).file_name = "FX_SEK_base";
	rf.at(1).file_name = "FX_SEK_term";
	tenors.at(0).file_name = "FX_SEK_fx";
	tenors.at(1).file_name = "FX_SEK_fxAvg";

	try {
		// Read
		rf.at(0).m_A = read_csv_matrix("FX_SEK_base.csv");
		rf.at(1).m_A = read_csv_matrix("FX_SEK_term.csv");
		tenors.at(0).m_A = read_csv_matrix("FX_SEK_fx.csv");
		tenors.at(1).m_A = read_csv_matrix("FX_SEK_fxAvg.csv");
		
	}
	catch (std::exception const&) {
		return 0;
	}
	
	// Truncate to minimum curve lengths
    for (int i{ 0 }; i < count_rf; ++i) {
		//rf.at(i).m_A_trunc = rf.at(i).m_A;
        rf.at(i).k = 6;
		rf.at(i).m_A.resize(rf.at(i).m_A.size1(), curve_length);
	}

	// Truncate to minimum curve lengths
    for (int i{ 0 }; i < count_tenor; ++i) {
		//rf.at(i).m_A_trunc = rf.at(i).m_A;
        tenors.at(i).k = 9;
		tenors.at(i).m_A.resize(tenors.at(i).m_A.size1(), curve_length);
	}
	
	/*
	if (tenors.at(0).m_A.size1() == 0) {
		return 0;
	}
	else
	{
		return 1;
	}
	*/

	placeholder_eigen(rf.at(0), eigen_algorithm, eval_eigen);
	placeholder_eigen(rf.at(1), eigen_algorithm, eval_eigen);
	placeholder_eigen(tenors.at(1), eigen_algorithm, eval_eigen);
	
	return status;
}


void placeholder_eigen(
	CurveCollection& curve_collection, 
	int const eigen_algorithm, 
	bool const eval_eigen,
    bool const save
) {
	curve_collection.m_diff = matrixOperations::abs(matrixOperations::diff_matrix(curve_collection.m_A_trunc));
	//curve_collection.m_diff = matrixOperations::diff_matrix(curve_collection.m_A_trunc);

	// Todo: Clean data here and process:
	//ublas::matrix_range<ublas::matrix<double>> m_diff_clean(
	//	curve_collection.m_diff,
	//	ublas::range(0, min(1500, curve_collection.m_diff.size1())),
	//	ublas::range(0, curve_collection.m_diff.size2())
	//);
	ublas::matrix<double> m_centered{ matrixOperations::center_matrix(curve_collection.m_diff) }; // Previously m_diff_clean

	//curve_collection.m_E = m_rf_E(m_centered.size2(), k);
	//ublas::vector<double> v_Lambda(k);

	enum class PCA_Algo { IRAM, BDCSVD };
	if (PCA_Algo(eigen_algorithm) == PCA_Algo::IRAM && curve_collection.k != 0) {
		FactorCalculation::iram(
			m_centered / sqrt(m_centered.size1() - 1),
            curve_collection.k,
			curve_collection.m_E, 
			curve_collection.v_Lambda
		);
	} else if (PCA_Algo(eigen_algorithm) == PCA_Algo::BDCSVD || curve_collection.k == 0) {
		FactorCalculation::eigen_bdcsvd(
			m_centered / sqrt(m_centered.size1() - 1), 
			curve_collection.k, 
			curve_collection.m_E, 
			curve_collection.v_Lambda
		);
	}

	curve_collection.m_delta_xi = FactorCalculation::compute_risk_factors(
		curve_collection.m_E, 
		curve_collection.m_diff
	);


	if (eval_eigen) { // Evaluate eigendecomposition
		curve_collection.v_norm_errors = FactorCalculation::eig_all_norm_errors(
			statisticsOperations::covm(curve_collection.m_diff), // Previously m_diff_clean
			curve_collection.m_E,
			curve_collection.v_Lambda
		);
	}
    
    if (save) {
        write_csv_matrix(curve_collection.m_E, curve_collection.file_path + "eigvec_" + curve_collection.file_name);
        write_csv_vector(curve_collection.v_Lambda, curve_collection.file_path + "eigval_" + curve_collection.file_name);
    }
}
