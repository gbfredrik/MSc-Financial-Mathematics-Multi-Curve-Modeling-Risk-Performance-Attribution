#include "pch.h"

#include "lhsd.h"

#include <numeric>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/io.hpp>
using namespace boost::numeric::ublas;


matrix<double> lhsd::lhsd_gen(matrix<double> const& U) {
	size_t N = U.size1(); // Number of simulations
	size_t n = U.size2(); // Number of risk factors

	matrix<double> r(N, n);
    matrix<double> V(N, n);

	// Calculate rank statistic of U
	for (int i = 0; i < n; i++) {
		column(r, i) = rank(column(U, i));
	}


	for (int i = 0; i < N; i++) {
		for (int j = 0; j < n; j++) {
			V(i, j) = (r(i, j) - 1 / 2) / N;
		}
	}
	return V;
}

// Sort columns after rank
vector<double> lhsd::rank(vector<double> const& U) {
	size_t N = U.size();
	vector<double> r(N);
	vector<double> rowid(N);
	
	rowid = sort(U);

	for (int i = 0; i < N; i++) {
		r(rowid(i)) = i;
	}

	return r;

}

// Generate indices used for the rank statistics
vector<double> lhsd::sort(vector<double> const& U) {

    vector<size_t> indices(U.size());

    // iota is from <numeric>, C++0x
    std::iota(indices.begin(), indices.end(), 0);

    std::sort(indices.begin(), indices.end(), [&U](std::size_t left, std::size_t right)
        {
            return U[left] < U[right];
        });

    return indices;
}
