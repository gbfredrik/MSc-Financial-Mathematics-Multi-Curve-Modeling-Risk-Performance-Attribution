#include "pch.h"
#include "varRed.h"
#include "mex.h"

#include <numeric>

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/io.hpp>



#include "../MathLibrary/matrixOperations.h"
#include "../MathLibrary/rvSim.h"
#include "../MathLibrary/statisticsOperations.h"

using namespace boost::numeric::ublas;

matrix<double> varRed::redVariance(matrix<double> U, std::string varRedType) {
	
	size_t k = U.size1(); // Number of risk factors
	size_t N = U.size2(); // Number of simulations
	matrix<double> V(k, N);

	if (varRedType == "lhsd") {
		V = lhsd_gen(U);
	}
	else if (varRedType == "none") {
		V = U;
	}

	return V;
}


matrix<double> varRed::lhsd_gen(matrix<double> const& U) {
	size_t k = U.size1(); // Number of risk factors
	size_t N = U.size2(); // Number of simulations

	matrix<double> r(k, N);
	matrix<double> V(k, N);


	// Calculate rank statistic of U
	for (size_t i = 0; i < k; i++) {
		row(r, i) = rank(row(U, i));
	}


	for (size_t i = 0; i < k; i++) {
		for (size_t j = 0; j < N; j++) {

			V(i, j) = (r(i, j) - (1.0 / 2.0)) / (N);

		}
	}

	return V;
}

// Sort columns after rank
vector<double> varRed::rank(vector<double> const& U) {
	size_t N = U.size();
	vector<double> r(N);
	vector<double> rowid(N);

	rowid = sort(U);

	for (size_t i = 0; i < N; i++) {
		r(rowid(i)) = i + 1.0;
	}

	return r;

}

// Generate indices used for the rank statistics
vector<double> varRed::sort(vector<double> const& U) {

	vector<size_t> indices(U.size());

	// iota is from <numeric>, C++0x
	std::iota(indices.begin(), indices.end(), 0);

	std::sort(indices.begin(), indices.end(), [&U](std::size_t left, std::size_t right)
		{
			return U[left] < U[right];
		});

	return indices;
}

