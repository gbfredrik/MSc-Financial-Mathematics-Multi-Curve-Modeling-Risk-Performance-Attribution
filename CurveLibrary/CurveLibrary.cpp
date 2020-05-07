// CurveLibrary.cpp : Defines the exported functions for the DLL.
//

#include "pch.h"
#include "framework.h"
#include "CurveLibrary.h" // Must be this order, otherwise functions are not exported!

#include "sample_handler.h"
#include "../MathLibrary/matrixOperations.h"
#include "../RiskFactorCalculation/FactorCalculation.h"

#include <boost/numeric/ublas/matrix.hpp>

LONG __stdcall squareXL(int x, int &y) {
#pragma EXPORT

	y += 100;
	return x * x;
}

BOOL __stdcall ir_measurement_multiXL(BOOL const save_curves) {
#pragma EXPORT

	bool status{ 0 };

	try	{
		boost::numeric::ublas::matrix<double> m_forward_curves_rf;
		boost::numeric::ublas::matrix<double> m_forward_curves_tenor;
		
		status = placeholder_ir_measurement_multi(m_forward_curves_rf, m_forward_curves_tenor);

		if (status == 1 && save_curves) { // If computation is successful and saving is toggled
			status = write_txt_matrix(m_forward_curves_rf, "forward_curves_rf.txt");
			status = status && write_txt_matrix(m_forward_curves_tenor, "forward_curves_tenor.txt");
		}
	} catch (const std::exception&) {
		return 0;
	}

	return status;
}



BOOL __stdcall run_all_multiXL(BOOL const compute_curves/*, ...*/) {
#pragma EXPORT

	bool status{ 0 };

	try {
		boost::numeric::ublas::matrix<double> m_forward_curves_rf{};
		boost::numeric::ublas::matrix<double> m_forward_curves_tenor{};
		if (compute_curves) { // Compute
			status = placeholder_ir_measurement_multi(m_forward_curves_rf, m_forward_curves_tenor);
			// TODO: Replace above to compute curve(s)

			if (status == 1) { // If computation is successful and saving is toggled
				status = write_txt_matrix(m_forward_curves_rf, "forward_curves_rf.txt");
				status = status && write_txt_matrix(m_forward_curves_tenor, "forward_curves_tenor.txt");
			}
		}
		else { // Read
			m_forward_curves_rf = read_txt_matrix("25x10950.txt");
			m_forward_curves_tenor = read_txt_matrix("25x10950.txt");
		}

		boost::numeric::ublas::matrix<double> m_diff{ matrixOperations::diff_matrix(m_forward_curves_rf) };
		boost::numeric::ublas::matrix<double> m_rf_centered{ matrixOperations::center_matrix(m_forward_curves_rf) };

		int k = 6;
		status = status && write_txt_matrix(m_rf_centered, "rf_centered.txt");
		boost::numeric::ublas::matrix<double> m_rf_E(m_diff.size2(), k);
		boost::numeric::ublas::vector<double> v_rf_Lambda(k);

		status = status && FactorCalculation::iram(m_rf_centered / std::sqrt(m_rf_centered.size1() - 1), k, m_rf_E, v_rf_Lambda);
		//status = status && FactorCalculation::eigen_bdcsvd(m_rf_centered, k, m_rf_E, v_rf_Lambda);
		
		status = status && write_txt_matrix(m_rf_E, "rf_vec.txt");
		status = status && write_txt_vector(v_rf_Lambda, "rf_val.txt");





	} catch (const std::exception&) {
		return 0;
	}

	return status;
}
