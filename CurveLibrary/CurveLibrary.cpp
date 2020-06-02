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

LONG __stdcall squareXL(int x, int& y) {
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
	bool const compute_curves,
	int const eigen_algorithm,
	bool const eval_eigen,
	double* return_norm_errors
) {
#pragma EXPORT

	bool status{ 1 };
	int count_rf{ 1 };
	int count_tenor{ 1 };

	std::vector<ublas::matrix<double>> m_rf(count_rf);
	std::vector<ublas::matrix<double>> m_tenors(count_tenor);
	try {
		if (compute_curves) { // Compute
			//status = placeholder_ir_measurement_multi(m_rf, m_forward_curves_tenor);
			// TODO: Replace above to compute curve(s)

			// If computation is successful
			for (int i{ 0 }; i < count_rf; ++i) {
				write_csv_matrix(
					m_rf.at(i), 
					"forward_curves_rf" + std::to_string(i) + ".csv"
				);
			}
			for (int i{ 0 }; i < count_tenor; ++i) {
				write_csv_matrix(
					m_tenors.at(i), 
					"forward_curves_tenor" + std::to_string(i) + ".csv"
				);
			}
		} else { // Read
			m_rf.at(0) = read_csv_matrix("fHist.csv");
			m_tenors.at(0) = read_csv_matrix("piHist.csv");
		}
	} catch (std::exception const&) {
		return -1;
	}
	/*
	ublas::matrix<double> m_diff{ matrixOperations::diff_matrix(m_forward_curves_rf) };
	// Todo: Clean data here and process:
	ublas::matrix_range<ublas::matrix<double>> m_diff_clean(
		m_diff, 
		ublas::range(0, min(1500, m_diff.size1())), 
		ublas::range(0, m_diff.size2())
	);
	ublas::matrix<double> m_rf_centered{ matrixOperations::center_matrix(m_diff_clean) };

	int k{ 6 };
	ublas::matrix<double> m_rf_E(m_rf_centered.size2(), k);
	ublas::vector<double> v_rf_Lambda(k);

	enum class PCA_Algo { IRAM, BDCSVD };
	//PCA_Algo pca_algo{ PCA_Algo::IRAM };

	if (PCA_Algo(eigen_algorithm) == PCA_Algo::IRAM) {
		status = status && FactorCalculation::iram(m_rf_centered / sqrt(m_rf_centered.size1() - 1), k, m_rf_E, v_rf_Lambda);
	} else if (PCA_Algo(eigen_algorithm) == PCA_Algo::BDCSVD) {
		status = status && FactorCalculation::eigen_bdcsvd(m_rf_centered / sqrt(m_rf_centered.size1() - 1), k, m_rf_E, v_rf_Lambda);
	}
		
	//status = status && write_csv_matrix(m_rf_E, "rf_vec.csv");
	//status = status && write_csv_vector(v_rf_Lambda, "rf_val.csv");

	ublas::matrix<double> m_rf_delta_xi{ FactorCalculation::compute_risk_factors(m_rf_E, m_diff) };
	//status = status && write_csv_matrix(m_rf_delta_xi, "rf_delta_xi.csv");
		
	if (eval_eigen) { // Evaluate eigendecomposition
		ublas::vector<double> v_errors(
			FactorCalculation::eig_all_norm_errors(
				statisticsOperations::covm(m_diff_clean), 
				m_rf_E, 
				v_rf_Lambda
			)
		);
		for (int i{ 0 }; i < k; ++i) {
			return_norm_errors[i] = v_errors(i);
		}
	}
		
	//...
	*/
	
	return status;
}


