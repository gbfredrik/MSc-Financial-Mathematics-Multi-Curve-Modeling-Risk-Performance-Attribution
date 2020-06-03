// CurveLibrary.cpp : Defines the exported functions for the DLL.
//

#include "pch.h"
#include "framework.h"
#include "CurveLibrary.h" // Must be this order, otherwise functions are not exported!

#include "sample_handler.h"
#include "../MathLibrary/matrixOperations.h"
#include "../MathLibrary/statisticsOperations.h"
#include "../RiskFactorCalculation/FactorCalculation.h"

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>

#include <algorithm>
#include <vector>

using namespace boost::numeric;

LONG __stdcall testXL(int x, int& y) {
#pragma EXPORT

	y += 100;

	return x * x;
}

BOOL __stdcall ir_measurement_multiXL(BOOL const save_curves) {
#pragma EXPORT

	bool status{ 0 };

	try	{
		ublas::matrix<double> m_forward_curves_rf;
		ublas::matrix<double> m_forward_curves_tenor;
		
		status = placeholder_ir_measurement_multi(m_forward_curves_rf, m_forward_curves_tenor);

		if (save_curves) { // If computation is successful and saving is toggled
			status = status && write_txt_matrix(m_forward_curves_rf, "forward_curves_rf.txt");
			status = status && write_txt_matrix(m_forward_curves_tenor, "forward_curves_tenor.txt");
		}
	} catch (const std::exception&) {
		return 0;
	}

	return status;
}

BOOL __stdcall run_all_multiXL(
	int const eigen_algorithm,
	bool const eval_eigen,
	double* return_norm_errors
) {
#pragma EXPORT

	bool status{ 1 };
	int count_tenor{ 1 };

	CurveCollection rf;
	std::vector<CurveCollection> tenors(count_tenor);

	rf.filename = "fHist.csv";
	for (CurveCollection& cc : tenors) {
		cc.filename = "piHist.csv";
	}

	try {
		// Read
		rf.m_A = read_csv_matrix("fHist.csv");
		for (CurveCollection& cc : tenors) {
			cc.m_A = read_csv_matrix("piHist.csv");
		}
	} catch (std::exception const&) {
		return -1;
	}

	int k{ 6 };
	placeholder_eigen(rf, eigen_algorithm, eval_eigen, k);
	placeholder_eigen(tenors.at(0), eigen_algorithm, eval_eigen, k);
	
	
	return status;
}

void placeholder_eigen(
	CurveCollection& curve_collection, 
	int const eigen_algorithm, 
	bool const eval_eigen,
	int const k
) {
	curve_collection.m_diff = matrixOperations::diff_matrix(curve_collection.m_A);

	// Todo: Clean data here and process:
	ublas::matrix_range<ublas::matrix<double>> m_diff_clean(
		curve_collection.m_diff,
		ublas::range(0, min(1500, curve_collection.m_diff.size1())),
		ublas::range(0, curve_collection.m_diff.size2())
	);
	ublas::matrix<double> m_centered{ matrixOperations::center_matrix(m_diff_clean) };

	//curve_collection.m_E =  m_rf_E(m_centered.size2(), k);
	//ublas::vector<double> v_Lambda(k);

	enum class PCA_Algo { IRAM, BDCSVD };
	if (PCA_Algo(eigen_algorithm) == PCA_Algo::IRAM) {
		FactorCalculation::iram(
			m_centered / sqrt(m_centered.size1() - 1),
			k, 
			curve_collection.m_E, 
			curve_collection.v_Lambda
		);
	} else if (PCA_Algo(eigen_algorithm) == PCA_Algo::BDCSVD) {
		FactorCalculation::eigen_bdcsvd(
			m_centered / sqrt(m_centered.size1() - 1), 
			k, 
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
			statisticsOperations::covm(m_diff_clean),
			curve_collection.m_E,
			curve_collection.v_Lambda
		);
	}

	write_csv_matrix(curve_collection.m_E, "eigvec_" + curve_collection.filename);
	write_csv_vector(curve_collection.v_Lambda, "eigval_" + curve_collection.filename);
}
