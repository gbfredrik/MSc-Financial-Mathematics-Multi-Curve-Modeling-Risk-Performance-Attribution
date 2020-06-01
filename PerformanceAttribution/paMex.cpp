#include "pch.h"
#include "mex.h"

#include "pa.h"
#include <boost/numeric/ublas/matrix.hpp>

using namespace boost::numeric::ublas;

/* The gateway function */
void mexFunction(int nlhs, mxArray* plhs[],
	int nrhs, const mxArray* prhs[]) {

	// Get sizes
	size_t I = mxGetNumberOfElements(prhs[1]); // num. of instruments
	size_t n = mxGetM(mxGetField(prhs[2], 0, mxGetFieldNameByNumber(prhs[2], 0))); // Number of disc points
	size_t M = mxGetNumberOfElements(prhs[2]); // num. curves
	vector<size_t> numFloat(I);
	vector<size_t> numFix(I);

	for (size_t i = 0; i < I; i++) {
		numFloat(i) = mxGetNumberOfElements(mxGetField(prhs[4], 0, mxGetFieldNameByNumber(prhs[4], i)));
		numFix(i) = mxGetNumberOfElements(mxGetField(prhs[5], 0, mxGetFieldNameByNumber(prhs[5], i)));
	}

	// Get N, y
	double* NMex;
	double* yMex;
	vector<double> y(I);
	vector<double> N(I);
	NMex = mxGetPr(prhs[0]);
	yMex = mxGetPr(prhs[1]);
	for (size_t i = 0; i < I; i++) {
		y(i) = yMex[i];
		N(i) = NMex[i];
	}

	// Get E
	vector<double*> EMex(M);
	vector<matrix<double>> E(M);
	for (size_t i = 0; i < n; i++) {
		EMex(i) = mxGetPr(mxGetField(prhs[2], 0, mxGetFieldNameByNumber(prhs[2], i)));
		E(i) = matrix<double>(n, n);
	}
	for (size_t i = 0; i < n; i++) {
		for (size_t j = 0; j < n; j++) {
			for (size_t l = 0; l < n; l++) {
				E(i)(j, l) = EMex(i)[j + l * n];
			}
		}
	}

	// Get kZero, kTau
	double* kMex;
	vector<size_t> k(M); // Number of risk factors for each curve
	kMex = mxGetPr(prhs[3]);
	for (size_t i = 0; i < M; i++) {
		k(i) = kMex[i];
	}
	
	// Get floatCashFlows, fixCashFlows
	vector<double*> floatCashFlowsMex;
	vector<double*> fixCashFlowsMex;
	vector<vector<double>> floatCashFlows(I);
	vector<vector<double>> fixCashFlows(I);
	for (int i = 0; i < I; i++) {
		floatCashFlowsMex(i) = mxGetPr(mxGetField(prhs[4], 0, mxGetFieldNameByNumber(prhs[4], i)));
		floatCashFlows(i) = vector<double>(numFloat(i));
		fixCashFlowsMex(i) = mxGetPr(mxGetField(prhs[5], 0, mxGetFieldNameByNumber(prhs[5], i)));
		fixCashFlows(i) = vector<double>(numFix(i));
		for (size_t j = 0; j < numFloat(i); j++) {
			floatCashFlows(i)(j) = floatCashFlowsMex(i)[j];
		}
		for (size_t j = 0; j < numFix(i); j++) {
			fixCashFlows(i)(j) = fixCashFlowsMex(i)[j];
		}
	}

	// Get curveData
	vector<double*> curveDataMex(M);
	vector<matrix<double>> curveData(M);
	

	



	/* call the computational routine */
	pa::performanceAttribution(N, y, E, k, floatCashFlows, fixCashFlows, curveData, 
		times, startDate, endDate NPV, carry, sumRiskFactors, epsI, epsA, epsP);








}