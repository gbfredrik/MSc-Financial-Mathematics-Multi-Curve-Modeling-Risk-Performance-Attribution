#include "pch.h"
#include "mex.h"

#include "pa.h"

#include "../MathLibrary/matrixOperations.h"
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/optional/optional.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/operation.hpp>
#include <boost/numeric/ublas/operation_blocked.hpp>

using namespace boost::numeric::ublas;

void pa::performanceAttribution(vector<double> const& N, vector<double> const& y, vector<matrix<double>> const& E,
	vector<int> k, vector<vector<int>> floatCashFlows, vector<vector<int>> fixCashFlows,
	vector<matrix<double>> const& curveData, vector<int> const& times, vector<int> const& startDate, 
	vector<int> const& endDate, vector<vector<double>>& NPV, vector<vector<double>>& carry, 
	vector<vector<double>>& sumRiskFactors, vector<vector<double>>& epsI, vector<vector<double>>& epsA, 
	vector<vector<double>>& epsP) {
	


	int I = y.size(); // num. of instruments
	int M = curveData.size(); // num. of curves
	int n = curveData(0).size2(); // num. of disc. points
	int dataLength = curveData(0).size1();
	vector<int> instrumentDays(I); // num. of in-sample days for each instrument
	for (int i = 0; i < I; i++) {
		instrumentDays(i) = endDate(i) - startDate(i);
	}
	
	// number of risk factors
	int kZero = k(0);
	int kPi = k(1);
	int kTot = kZero + kPi;
	
	// set risk principal component sizes
	vector<double> dXiZero(kTot);
	vector<double> dXiPi(kTot);
	vector<double> XiBarZero(2 * n);
	vector<double> XiBarPi(2 * n);
	vector<double> XiBarZeroPrev(2 * n);
	vector<double> XiBarPiPrev(2 * n);
	vector<double> XiBardXiZero(2 * n);
	vector<double> XiBardXiPi(2 * n);
	
	// slice E matrix and get E_k
	vector<matrix<double>> E_k(M);
	for (int i = 0; i < M; i++) {
		E_k(i) = matrix<double>(n, k(i));
		E_k(i) = subrange(E(i), 0, 0, n-1, k(i)-1);
	}
	mexPrintf("wat ");
	// set data matrices
	mexPrintf("%d", dataLength);
	mexPrintf(" ");
	mexPrintf("%d", n);
	mexPrintf(" ");

	matrix<double> f(n, dataLength);
	matrix<double> pi(n, dataLength);
	f = trans(curveData(0));
	pi = trans(curveData(1));
	mexPrintf("wat1 ");
	
	// set the integrating matrix and spot data matrices
	matrix<double> A(n, n);
	intMatrix(A);
	mexPrintf("wat2 ");
	
	matrix<double> r(n, dataLength);
	matrix<double> piSpot(n, dataLength);
	mexPrintf("%d", n);
	mexPrintf(" ");
	mexPrintf("%d", dataLength);
	mexPrintf(" ");
	
	axpy_prod(A, f, r);
	mexPrintf("wa3 ");
	axpy_prod(A, pi, piSpot);
	matrix<double> aZero(n, 2 * n);
	matrix<double> aPi(n, 2 * n);
	matrix<double> aZero_k(kZero, kTot);
	matrix<double> aPi_k(kPi, kTot);
	
	seta(A, E, n, n, aZero, "zero");
	mexPrintf("wat4 ");
	/*
	seta(A, E, n, n, aPi, "pi");
	seta(A, E_k, kZero, kPi, aZero_k, "zero");
	seta(A, E_k, kZero, kPi, aPi_k, "pi");
	mexPrintf("wat ");
	// gradient, hessian
	vector<double> g(kTot);
	matrix<double> H(kTot, kTot);
	mexPrintf("wat ");
	// date variables
	vector<double> deltaTj;
	vector<double> deltaTjPrev;
	-int postFirstDate = 0;
	int dateDiff = 0;
	vector<double> floatCashFlowsPrev;
	vector<double> fixCashFlowsPrev; 
	mexPrintf("wat ");
	
	// Loop over contracts
	for (int i = 0; i < I; i++) {

		deltaTj = vector<double>(fixCashFlows(i).size());		
		deltaTj = fixCashFlows(i) / 360;
		

		// Loop over in-sample days
		for (int j = 0; j < instrumentDays(i); j++) {
			if (j == 0) { // if first in-sample day
				setXiBar(E, n, column(f, j), column(pi, j), XiBarZero, "zero");
				setXiBar(E, n, column(f, j), column(pi, j), XiBarPi, "pi");
				carry(i)(j) = 0;
				sumRiskFactors(i)(j) = 0;
				epsI(i)(j) = 0;
				epsA(i)(j) = 0;
				epsP(i)(j) = irsPrice(N(i), y(i), floatCashFlows(i), fixCashFlows(i), deltaTj, column(r, j), column(piSpot, j));
				NPV(i)(j) = epsP(i)(j);
			}
			else {
				dateDiff = times(j) - times(j - 1);
				floatCashFlowsPrev = vector<double>(floatCashFlows(i).size());
				floatCashFlowsPrev = floatCashFlows(i);
				fixCashFlowsPrev = vector<double>(fixCashFlows(i).size());
				fixCashFlowsPrev = fixCashFlows(i);
				deltaTjPrev = vector<double>(deltaTj.size());
				deltaTjPrev = deltaTj;
				
				if (postFirstDate == 0) {
					for (int l = 0; l < floatCashFlows(i).size(); l++) {
						floatCashFlows(i)(l) = floatCashFlows(i)(l) - dateDiff;
					}
					for (int l = 0; l < fixCashFlows(i).size(); l++) {
						fixCashFlows(i)(l) = fixCashFlows(i)(l) - dateDiff;
					}
					if (floatCashFlows(i)(0) == 0) {
						floatCashFlows(i)(0) = 0;
						postFirstDate = 1;
					}
				}
				else if (postFirstDate == 1) { // Todo: fixa fixcashflows också
					for (int l = 1; l < floatCashFlows(i).size(); l++) {
						floatCashFlows(i)(l) = floatCashFlows(i)(l) - dateDiff;
					}
					for (int l = 0; l < fixCashFlows(i).size(); l++) {
						fixCashFlows(i)(l) = fixCashFlows(i)(l) - dateDiff;
					}
					if (floatCashFlows(i)(1) == 0) {
						subrange(floatCashFlows(i), 0, floatCashFlows(i).size() - 2) = subrange(floatCashFlows(i), 1, floatCashFlows(i).size() - 1);
						floatCashFlows(i).resize(floatCashFlows(i).size() - 1, true);
					}
				}

				if (deltaTj.size() != fixCashFlows(i).size()) {
					deltaTj.resize(fixCashFlows(i).size(), false);
				}
				deltaTj = fixCashFlows(i) / 360;

				XiBarZeroPrev = XiBarZero;
				XiBarPiPrev = XiBarZeroPrev;
				setdXi(E, kZero, kPi, column(f, j), column(pi, j), column(f, j - 1), column(pi, j - 1), XiBarZero, "zero");
				setdXi(E, kZero, kPi, column(f, j), column(pi, j), column(f, j - 1), column(pi, j - 1), XiBarPi, "pi");
				setXiBar(E, n, column(f, j), column(pi, j), XiBarZero, "zero");
				setXiBar(E, n, column(f, j), column(pi, j), XiBarPi, "pi");
				XiBardXiZero = XiBarZeroPrev;
				subrange(XiBardXiZero, 0, kZero-1) = subrange(XiBarZeroPrev, 0, kZero-1) + dXiZero;
				XiBardXiPi = XiBarPiPrev;
				subrange(XiBardXiPi, n, n + kPi) = subrange(XiBarPiPrev, n, n + kPi) + dXiPi;

				// ------------------------ GRADIENT AND HESSIAN ------------------
				g = grad(N(i), y(i), floatCashFlows(i), fixCashFlows(i), aZero_k, aPi_k, deltaTj,
					column(f, j), column(pi, j));
				H = hess(N(i), y(i), floatCashFlows(i), fixCashFlows(i), aZero_k, aPi_k, deltaTj,
					column(f, j), column(pi, j));				
				// ------------------------ CARRY. --------------------------------
				carry(i)(j) = irsPriceRiskFactor(N(i), y(i), floatCashFlows(i), fixCashFlows(i), aZero, aPi, deltaTj, XiBarZeroPrev, XiBarPiPrev)
					- irsPriceRiskFactor(N(i), y(i), floatCashFlowsPrev, fixCashFlowsPrev, aZero, aPi, deltaTjPrev, XiBarZeroPrev, XiBarPiPrev);
				// ------------------------ TRUNC. ERR. ----------------------------
				epsI(i)(j) = irsPrice(N(i), y(i), floatCashFlows(i), fixCashFlows(i), deltaTj, column(r, j), column(piSpot, j))
					- irsPriceRiskFactor(N(i), y(i), floatCashFlows(i), fixCashFlows(i), aZero, aPi, deltaTj, XiBardXiZero, XiBardXiPi);
				// ------------------------ TAYL. APPROX. ERR. ---------------------
				epsA(i)(j) = irsPriceRiskFactor(N(i), y(i), floatCashFlows(i), fixCashFlows(i), aZero, aPi, deltaTj, XiBardXiZero, XiBardXiPi)
					- irsPriceRiskFactor(N(i), y(i), floatCashFlows(i), fixCashFlows(i), aZero, aPi, deltaTj, XiBarZeroPrev, XiBarPiPrev)
					- inner_prod(g, dXiZero + dXiPi) - (1/2) * inner_prod(dXiZero + dXiPi, prod(H, dXiZero + dXiPi));
				// ------------------------ RISK FACTORS ---------------------------
				sumRiskFactors(i)(j) = inner_prod(g, dXiZero + dXiPi) + (1 / 2) * inner_prod(dXiZero + dXiPi, prod(H, dXiZero + dXiPi));
				// ------------------------ NPV ------------------------------------
				NPV(i)(j) = irsPrice(N(i), y(i), floatCashFlows(i), fixCashFlows(i), deltaTj, column(r, j), column(piSpot, j))
					- irsPrice(N(i), y(i), floatCashFlowsPrev, fixCashFlowsPrev, deltaTjPrev, column(r, j - 1), column(piSpot, j - 1));
			}
		}
	}
	*/
}

void pa::intMatrix(matrix<double>& A) {
	mexPrintf("hej ");
	int n = A.size1();

	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			if (j <= i) {
				A(i, j) = 1 / (i + 1);
			}
			else {
				A(i, j) = 0;
			}
		}
	}
	mexPrintf("då ");
}

vector<double> pa::grad(double N, double y, vector<int> floatCashFlows, vector<int> fixCashFlows, matrix<double> aZero,
	matrix<double> aPi, vector<double> deltaTj, vector<double> r, vector<double> pi) {
	
	int numCashFix = fixCashFlows.size();
	int numCashFloat = floatCashFlows.size();

	int kTot = aZero.size2();

	vector<double> fixL(kTot, 0);
	// Calc fix leg
	for (int i = 0; i < numCashFix; i++) {
		fixL = fixL + deltaTj(i) * ((floatCashFlows(0) / 365 * row(aZero, floatCashFlows(0)) - fixCashFlows(i) / 365 * row(aZero, fixCashFlows(i)))
			* exp(floatCashFlows(0) / 365 * r(floatCashFlows(0)) - fixCashFlows(i) / 365 * r(fixCashFlows(i))));

	}
	fixL = N * y * fixL;

	// Calc float leg
	vector<double> floatL(kTot, 0);
	for (int i = 0; i < numCashFloat - 1; i++) {
		floatL = floatL + (floatCashFlows(i + 1) / 365 * row(aPi, floatCashFlows(i + 1))
			- floatCashFlows(i) / 365 * (row(aZero, floatCashFlows(i)) + row(aPi, floatCashFlows(i)))
			+ floatCashFlows(0) / 365 * row(aZero, floatCashFlows(0)))
			* exp(floatCashFlows(i + 1) / 365 * pi(floatCashFlows(i + 1))
				- floatCashFlows(i) / 365 * (r(floatCashFlows(i)) + pi(floatCashFlows(i)))
				+ floatCashFlows(0) / 365 * r(floatCashFlows(0)))
			- (floatCashFlows(0) / 365 * row(aZero, floatCashFlows(0))
				- floatCashFlows(i + 1) / 365 * row(aZero, floatCashFlows(i + 1)))
			* exp(floatCashFlows(0) / 365 * r(floatCashFlows(0))
				- floatCashFlows(i + 1) / 365 * r(floatCashFlows(i + 1)));
	}
	
	floatL = N * floatL;

	return floatL - fixL;
}

matrix<double> pa::hess(double N, double y, vector<int> floatCashFlows, vector<int> fixCashFlows, matrix<double> aZero,
	matrix<double> aPi, vector<double> deltaTj, vector<double> r, vector<double> pi) {

	int numCashFix = fixCashFlows.size();
	int numCashFloat = floatCashFlows.size();

	int kTot = aZero.size2();

	matrix<double> fixL(kTot, kTot, 0);
	// Calc fix leg
	for (int i = 0; i < numCashFix; i++) {
		fixL = fixL + deltaTj(i) * (outer_prod((floatCashFlows(0) / 365 * row(aZero, floatCashFlows(0)) - fixCashFlows(i) / 365 * row(aZero, fixCashFlows(i)))
			, (floatCashFlows(0) / 365 * row(aZero, floatCashFlows(0)) - fixCashFlows(i) / 365 * row(aZero, fixCashFlows(i))))
			* exp(floatCashFlows(0) / 365 * r(floatCashFlows(0)) - fixCashFlows(i) / 365 * r(fixCashFlows(i))));

	}
	fixL = N * y * fixL;

	// Calc float leg
	matrix<double> floatL(kTot, kTot, 0);
	for (int i = 0; i < numCashFloat - 1; i++) {
		floatL = floatL + outer_prod((floatCashFlows(i + 1) / 365 * row(aPi, floatCashFlows(i + 1)) 
			- floatCashFlows(i) / 365 * (row(aZero, floatCashFlows(i)) + row(aPi, floatCashFlows(i))) + floatCashFlows(0) / 365 * row(aZero, floatCashFlows(0)))
			, (floatCashFlows(i + 1) / 365 * row(aPi, floatCashFlows(i + 1)) - floatCashFlows(i) / 365 * (row(aZero, floatCashFlows(i)) + row(aPi, floatCashFlows(i)))
			+ floatCashFlows(0) / 365 * row(aZero, floatCashFlows(0))))
			* exp(floatCashFlows(i + 1) / 365 * pi(floatCashFlows(i + 1)) - floatCashFlows(i) / 365 * (r(floatCashFlows(i)) + pi(floatCashFlows(i)))
				+ floatCashFlows(0) / 365 * r(floatCashFlows(0)))
			- outer_prod((floatCashFlows(0) / 365 * row(aZero, floatCashFlows(0)) - floatCashFlows(i + 1) / 365 * row(aZero, floatCashFlows(i + 1)))
			, (floatCashFlows(0) / 365 * row(aZero, floatCashFlows(0)) - floatCashFlows(i + 1) / 365 * row(aZero, floatCashFlows(i + 1))))
			* exp(floatCashFlows(0) / 365 * r(floatCashFlows(0)) - floatCashFlows(i + 1) / 365 * r(floatCashFlows(i + 1)));

	}
	floatL = N * floatL;

	return floatL - fixL;

}

void pa::setdXi(vector<matrix<double>> const& E, int kZero, int kPi, vector<double> const& f, vector<double> const& pi, 
	vector<double> const& fPrev, vector<double> const& piPrev, vector<double>& XiBar, std::string type) {
	
	vector<double> XiBarZero(kZero);
	XiBarZero = prod(trans(E(0)), f - fPrev);
	vector<double> XiBarPi(kPi);
	XiBarPi = prod(trans(E(1)), pi - piPrev);
	
	if (type == "zero") {
		vector<double> zeroVector(kPi, 0);
		subrange(XiBar, 0, kZero - 1) = XiBarZero;
		subrange(XiBar, kZero, kZero + kPi - 1) = zeroVector;
	}
	else if (type == "pi") {
		vector<double> zeroVector(kZero, 0);
		subrange(XiBar, 0, kZero - 1) = zeroVector;
		subrange(XiBar, kZero, kZero + kPi - 1) = XiBarPi;
	}
}

void pa::setXiBar(vector<matrix<double>> const& E, int n, vector<double> const& f, vector<double> const& pi,
	vector<double>& XiBar, std::string type) {
	
	vector<double> XiBarZero(n);
	XiBarZero = prod(trans(E(0)), f);
	vector<double> XiBarPi(n);
	XiBarPi = prod(trans(E(1)), pi);
	vector<double> zeroVector(n, 0);
	

	if (type == "zero") {
		subrange(XiBar, 0, n - 1) = XiBarZero;
		subrange(XiBar, n, 2 * n - 1) = zeroVector;
	}
	else if (type == "pi") {
		subrange(XiBar, 0, n - 1) = zeroVector;
		subrange(XiBar, n, 2 * n - 1) = XiBarPi;
	}
}

double pa::irsPrice(double N, double y, vector<int> const& floatCashFlows, vector<int> const& fixCashFlows,
	vector<double> const& deltaTj, vector<double> const& r, vector<double> const& pi) {
	
	int numCashFix = fixCashFlows.size();
	int numCashFloat = floatCashFlows.size();

	// Fix leg
	double fixL = 0;
	
	for (int i = 0; i < numCashFix; i++) {
		fixL = fixL + deltaTj(i) * exp(floatCashFlows(0) / 365 * r(floatCashFlows(0)) - fixCashFlows(i) / 365 * r(fixCashFlows(i)));

	}
	fixL = N * y * fixL;

	// Float leg
	double floatL = 0;
	for (int i = 0; i < numCashFloat - 1; i++) {
		floatL = floatL + exp(floatCashFlows(i + 1) / 365 * pi(floatCashFlows(i + 1)) - floatCashFlows(i) / 365 * (r(floatCashFlows(i)) + pi(floatCashFlows(i)))
			+ floatCashFlows(0) / 365 * r(floatCashFlows(0)))
			- exp(floatCashFlows(0) / 365 * r(floatCashFlows(0)) - floatCashFlows(i + 1) / 365 * r(floatCashFlows(i + 1)));
	}
	floatL = N * floatL;

	return floatL - fixL;
}

double pa::irsPriceRiskFactor(double N, double y, vector<int> const& floatCashFlows, vector<int> const& fixCashFlows,
	matrix<double> const& aZero, matrix<double> aPi, vector<double> const& deltaTj, vector<double> const& XiZero,
	vector<double> XiPi) {

	int numCashFix = fixCashFlows.size();
	int numCashFloat = floatCashFlows.size();

	double fixL = 0;
	// Calc fix leg
	for (int i = 0; i < numCashFix; i++) {
		fixL = fixL + deltaTj(i) * exp(floatCashFlows(0) / 365 * inner_prod(row(aZero, floatCashFlows(0)), XiZero)
			- fixCashFlows(i) / 365 * inner_prod(row(aZero, fixCashFlows(i)), XiZero));
	}
	fixL = N * y * fixL;
	
	// Calc float leg
	double floatL = 0;
	for (int i = 0; i < numCashFloat - 1; i++) {
		floatL = floatL + exp(floatCashFlows(i + 1) / 365 * inner_prod(row(aPi, floatCashFlows(i + 1)), XiPi)
			- floatCashFlows(i) / 365 * (inner_prod(row(aZero, floatCashFlows(i)), XiZero) + inner_prod(row(aPi, floatCashFlows(i)), XiPi))
			+ floatCashFlows(0) / 365 * inner_prod(row(aZero, floatCashFlows(0)), XiZero))
			- exp(floatCashFlows(0) / 365 * inner_prod(row(aZero, floatCashFlows(0)), XiZero)
				- floatCashFlows(i + 1) / 365 * inner_prod(row(aZero, floatCashFlows(i + 1)), XiZero));
	
	
	}
	floatL = N * floatL;

	return floatL - fixL;

}

void pa::seta(matrix<double> const& A, vector<matrix<double>> const& E, int kZero, int kPi, matrix<double>& a, std::string type) {

	int n = A.size1();

	if (type == "zero") {
		matrix<double> aTemp1(n, kZero);
		axpy_prod(A, E(0), aTemp1);
		matrix<double> aTemp2(n, kPi, 0);

		subrange(a, 0, 0, n - 1, kZero-1) = aTemp1;
		subrange(a, n, kZero, n, kZero + kPi - 1) = aTemp2;

	}
	else if (type == "pi") {
		matrix<double> aTemp3(n, kZero, 0);
		matrix<double> aTemp4(n, kPi);
		axpy_prod(A, E(1), aTemp4);

		subrange(a, 0, 0, n - 1, kZero - 1) = aTemp3;
		subrange(a, n, kZero, n, kZero + kPi - 1) = aTemp4;
	}

	
}

