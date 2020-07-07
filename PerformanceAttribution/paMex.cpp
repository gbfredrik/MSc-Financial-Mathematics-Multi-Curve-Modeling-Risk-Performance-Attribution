#include "pch.h"
#include "mex.h"

#include "pa.h"

#include "../MathLibrary/matrixOperations.h"
#include <boost/numeric/ublas/matrix.hpp>

using namespace boost::numeric::ublas;

/* The gateway function */
void mexFunction(int nlhs, mxArray* plhs[],
	int nrhs, const mxArray* prhs[]) {

	// Get sizes
	int I = mxGetNumberOfElements(prhs[1]); // num. of instruments
	int n = mxGetM(mxGetField(prhs[2], 0, mxGetFieldNameByNumber(prhs[2], 0))); // Number of disc points
	int days = mxGetM(mxGetField(prhs[6], 0, mxGetFieldNameByNumber(prhs[6], 0))); // Number of in sample days
	int M = mxGetNumberOfFields(prhs[2]); // number of curves
	mexPrintf("%d", M);
	mexPrintf(" ");
	vector<int> numFloat(I); // number of floating payments for each instrument
	vector<int> numFix(I); // number of fix payments for each instrument
	for (int i = 0; i < I; i++) {
		numFloat(i) = mxGetNumberOfElements(mxGetField(prhs[4], 0, mxGetFieldNameByNumber(prhs[4], i)));
		numFix(i) = mxGetNumberOfElements(mxGetField(prhs[5], 0, mxGetFieldNameByNumber(prhs[5], i)));
	}
	
	// Get N, y
	double* NMex;
	double* yMex;
	vector<double> y(I); // save each instrument´s yield
	vector<double> N(I); // save each instrument´s nominal amount
	NMex = mxGetPr(prhs[0]);
	yMex = mxGetPr(prhs[1]);
	for (int i = 0; i < I; i++) {
		y(i) = yMex[i];
		N(i) = NMex[i];
	}
	
	// Get E
	vector<double*> EMex(M);
	vector<matrix<double>> E(M);
	for (int i = 0; i < M; i++) {
		EMex(i) = mxGetPr(mxGetField(prhs[2], 0, mxGetFieldNameByNumber(prhs[2], i)));
		E(i) = matrix<double>(n, n);
	}
	for (int i = 0; i < M; i++) {
		for (int j = 0; j < n; j++) {
			for (int l = 0; l < n; l++) {
				E(i)(j, l) = EMex(i)[j + l * n];
			}
		}
	}

	// Get kZero, kTau
	double* kMex;
	vector<int> k(M); // Number of risk factors for each curve
	kMex = mxGetPr(prhs[3]);
	for (int i = 0; i < M; i++) {
		k(i) = kMex[i];
	}
	
	// Get floatCashFlows, fixCashFlows
	vector<double*> floatCashFlowsMex(I);
	vector<double*> fixCashFlowsMex(I);
	vector<vector<int>> floatCashFlows(I);
	vector<vector<int>> fixCashFlows(I);
	for (int i = 0; i < I; i++) {
		floatCashFlowsMex(i) = mxGetPr(mxGetField(prhs[4], 0, mxGetFieldNameByNumber(prhs[4], i)));
		floatCashFlows(i) = vector<int>(numFloat(i));
		fixCashFlowsMex(i) = mxGetPr(mxGetField(prhs[5], 0, mxGetFieldNameByNumber(prhs[5], i)));
		fixCashFlows(i) = vector<int>(numFix(i));
		for (int j = 0; j < numFloat(i); j++) {
			floatCashFlows(i)(j) = floatCashFlowsMex(i)[j];
		}
		for (int j = 0; j < numFix(i); j++) {
			fixCashFlows(i)(j) = fixCashFlowsMex(i)[j];
		}
	}

	// Get curveData
	vector<double*> curveDataMex(M);
	vector<matrix<double>> curveData(M);
	for (int i = 0; i < M; i++) {
		curveDataMex(i) = mxGetPr(mxGetField(prhs[6], 0, mxGetFieldNameByNumber(prhs[6], i)));
		curveData(i) = matrix<double>(days, n);
		for (int j = 0; j < days; j++) {
			for (int l = 0; l < n; l++) {
				curveData(i)(j, l) = curveDataMex(i)[j + l * days];
			}
		}
	}
	
	// Get times
	double* timesMex;
	vector<int> times(days);
	timesMex = mxGetPr(prhs[7]);
	for (int i = 0; i < days; i++) {
		times(i) = timesMex[i];
	}

	// Get startDate, endDate
	double* startDateISMex;
	double* endDateISMex;
	vector<int> startDateIS(I);
	vector<int> endDateIS(I);
	startDateISMex = mxGetPr(prhs[8]);
	endDateISMex = mxGetPr(prhs[9]);
	int maxEndDate = 0;
	for (int i = 0; i < I; i++) {
		startDateIS(i) = startDateISMex[i];
		endDateIS(i) = endDateISMex[i];
		if (endDateIS(i) > maxEndDate) {
			maxEndDate = endDateIS(i);
		}
	}


	// Set output
	double* NPVMex;
	double* carryMex;
	double* sumRiskFactorsMex;
	double* epsIMex;
	double* epsAMex;
	double* epsPMex;
	vector<vector<double>> NPV(I);
	vector<vector<double>> carry(I);
	vector<vector<double>> sumRiskFactors(I);
	vector<vector<double>> epsI(I);
	vector<vector<double>> epsA(I);
	vector<vector<double>> epsP(I);
	for (int i = 0; i < nlhs; i++) {
		plhs[i] = mxCreateDoubleMatrix(maxEndDate - startDateIS(0), I, mxREAL);
	}
	NPVMex = mxGetPr(plhs[0]);
	carryMex = mxGetPr(plhs[1]);
	sumRiskFactorsMex = mxGetPr(plhs[2]);
	epsIMex = mxGetPr(plhs[3]);
	epsAMex = mxGetPr(plhs[4]);
	epsPMex = mxGetPr(plhs[5]);
	for (int i = 0; i < I; i++) {
		NPV(i) = vector<double>(endDateIS(i) - startDateIS(i));
		carry(i) = vector<double>(endDateIS(i) - startDateIS(i));
		sumRiskFactors(i) = vector<double>(endDateIS(i) - startDateIS(i));
		epsI(i) = vector<double>(endDateIS(i) - startDateIS(i));
		epsA(i) = vector<double>(endDateIS(i) - startDateIS(i));
		epsP(i) = vector<double>(endDateIS(i) - startDateIS(i));
	}


	/* call the computational routine */
	pa::performanceAttribution(N, y, E, k, floatCashFlows, fixCashFlows, curveData, 
		times, startDateIS, endDateIS, NPV, carry, sumRiskFactors, epsI, epsA, epsP);
	
	/*Convert result*/
	for (int i = 0; i < I; i++) {
		for (int j = 0; j < maxEndDate - startDateIS(0); j++) {
			NPV(i)(j) = NPVMex[i + j * I];
			carry(i)(j) = carryMex[i + j * I];
			sumRiskFactors(i)(j) = sumRiskFactorsMex[i + j * I];
			epsI(i)(j) = epsIMex[i + j * I];
			epsA(i)(j) = epsAMex[i + j * I];
			epsP(i)(j) = epsPMex[i + j * I];
			if (j >= endDateIS(i) - startDateIS(0)) {
				NPV(i)(j) = 0;
				carry(i)(j) = 0;
				sumRiskFactors(i)(j) = 0;
				epsI(i)(j) = 0;
				epsA(i)(j) = 0;
				epsP(i)(j) = 0;
			}
		}
	}
}