// CurveLibrary.cpp : Defines the exported functions for the DLL.
//

#include "pch.h"
#include "framework.h"
#include "CurveLibrary.h" // Must be this order, otherwise functions are not exported!

#include "sample_handler.h"
#include "../MathLibrary/matrixOperations.h"

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

	try	{
		boost::numeric::ublas::matrix<double> m_forward_curves_rf{};
		boost::numeric::ublas::matrix<double> m_forward_curves_tenor{};
		if (compute_curves) { // Compute
			status = placeholder_ir_measurement_multi(m_forward_curves_rf, m_forward_curves_tenor);
			// TODO: Replace above to compute curve(s)

			if (status == 1) { // If computation is successful and saving is toggled
				status = write_txt_matrix(m_forward_curves_rf, "forward_curves_rf.txt");
				status = status && write_txt_matrix(m_forward_curves_tenor, "forward_curves_tenor.txt");
			}
		} else { // Read
			m_forward_curves_rf = read_txt_matrix("25x10950.txt");
			m_forward_curves_tenor = read_txt_matrix("25x10950.txt");
		}

		// ...
		// ...
		boost::numeric::ublas::matrix<double> m_diff{matrixOperations::diff_matrix(m_forward_curves_rf)};
		//status = status && write_txt_matrix(m_diff, "rf_diff.txt");

		int k = 6;
		//boost::numeric::ublas::matrix<double> m_E{ matrixOperations::diff_matrix(m_forward_curves_rf) };


	} catch (const std::exception&) {
		return 0;
	}

	return status;
}
