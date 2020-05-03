#include "pch.h"
#include "MultipleYieldSim.h"
#include <tuple>
#include "mex.h"

#include <boost/numeric/ublas/matrix.hpp>


using namespace boost::numeric::ublas;

std::tuple<matrix<double>, matrix<double>> MultipleYieldSim::simMultiple(vector<double> eps, matrix<double> Ezero, matrix<double> Etau, vector<double> fZeroCurr, vector<double> fTauCurr, int N) {
	size_t n = fZeroCurr.size();
	matrix<double> fZeroNext(N, n);
	matrix<double> fTauNext(N, n);

	fZeroNext(0, 0) = 1;
	fTauNext(1, 2) = 3;

	return std::make_tuple(fZeroNext, fTauNext);

}

void MultipleYieldSim::testSim(double x, double* y, double* z, mwSize n) {
	mwSize i;
	/* multiply each element y by x */
	for (i = 0; i < n; i++) {
		z[i] = x * y[i];
	}
}