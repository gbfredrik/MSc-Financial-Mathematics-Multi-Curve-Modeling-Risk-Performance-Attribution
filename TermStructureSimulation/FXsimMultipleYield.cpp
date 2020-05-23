#include "pch.h"
#include "mex.h"
#include "MultipleYieldSim.h"
#include <boost/numeric/ublas/matrix.hpp>
#include <string>
#include "../MathLibrary/rvSim.h"

using namespace boost::numeric::ublas;

/* The gateway function */
void mexFunction(int nlhs, mxArray* plhs[],
	int nrhs, const mxArray* prhs[]) {

	// Get sizes
	size_t n = mxGetM(mxGetField(prhs[0], 0, mxGetFieldNameByNumber(prhs[0], 0))); // Number of disc points
	size_t m = mxGetM(mxGetField(prhs[6], 0, mxGetFieldNameByNumber(prhs[6], 0))); // Number of historical dates
	size_t M = nlhs; // Number of curves
	vector<size_t> k(M); // Number of risk factors for each curve
	for (size_t i = 0; i < M; i++) {
		k(i) = mxGetN(mxGetField(prhs[0], 0, mxGetFieldNameByNumber(prhs[0], i)));
	}

	// Get E, rho
	vector<double*> EMex(M);
	vector<double*> rhoMex(M);
	vector<matrix<double>> E(M);
	vector<matrix<double>> rho(M);
	for (size_t i = 0; i < M; i++) {
		EMex(i) = mxGetPr(mxGetField(prhs[0], 0, mxGetFieldNameByNumber(prhs[0], i)));	
		rhoMex(i) = mxGetPr(mxGetField(prhs[1], 0, mxGetFieldNameByNumber(prhs[1], i)));
		E(i) = matrix<double>(n, k(i));
		rho(i) = matrix<double>(k(i), k(i));
	}

	for (size_t i = 0; i < M; i++) {
		for (size_t j = 0; j < n; j++) {
			for (size_t l = 0; l < k(i); l++) {
				E(i)(j, l) = EMex(i)[j + l * n];
			}
		}
	}

	for (size_t i = 0; i < M; i++) {
		for (size_t j = 0; j < k(i); j++) {
			for (size_t l = 0; l < k(i); l++) {
				rho(i)(j, l) = rhoMex(i)[j + l * k(i)];
			}
		}
	}

	// Get mu, omega, alpha, beta, gamma
	vector<double*> muMex(M);
	vector<double*> omegaMex(M);
	vector<double*> alphaMex(M);
	vector<double*> betaMex(M);
	vector<double*> gammaMex(M);
	vector<double*> dfCMex(M);
	vector<double*> dfMMex(M);
	vector<double*> xiHatMex(M);
	vector<vector<double>> mu(M);
	vector<vector<double>> omega(M);
	vector<vector<double>> alpha(M);
	vector<vector<double>> beta(M);
	vector<vector<double>> gamma(M);
	vector<vector<double>> dfC(M);
	vector<vector<double>> dfM(M);
	vector<vector<double>> xiHat(M);
	for (size_t i = 0; i < M; i++) {
		muMex(i) = mxGetPr(mxGetField(prhs[2], 0, mxGetFieldNameByNumber(prhs[2], i)));
		omegaMex(i) = mxGetPr(mxGetField(prhs[3], 0, mxGetFieldNameByNumber(prhs[3], i)));
		alphaMex(i) = mxGetPr(mxGetField(prhs[4], 0, mxGetFieldNameByNumber(prhs[4], i)));
		betaMex(i) = mxGetPr(mxGetField(prhs[5], 0, mxGetFieldNameByNumber(prhs[5], i)));
		gammaMex(i) = mxGetPr(mxGetField(prhs[13], 0, mxGetFieldNameByNumber(prhs[13], i)));
		xiHatMex(i) = mxGetPr(mxGetField(prhs[15], 0, mxGetFieldNameByNumber(prhs[15], i)));
		dfCMex(i) = mxGetPr(mxGetField(prhs[16], 0, mxGetFieldNameByNumber(prhs[16], i)));
		dfMMex(i) = mxGetPr(mxGetField(prhs[17], 0, mxGetFieldNameByNumber(prhs[17], i)));
		mu(i) = vector<double>(k(i));
		omega(i) = vector<double>(k(i));
		alpha(i) = vector<double>(k(i));
		beta(i) = vector<double>(k(i));
		gamma(i) = vector<double>(k(i));
		xiHat(i) = vector<double>(k(i));
		dfC(i) = vector<double>(k(i));
		dfM(i) = vector<double>(k(i));
	}
	for (size_t i = 0; i < M; i++) {
		for (size_t j = 0; j < k(i); j++) {
			mu(i)(j) = muMex(i)[j];
			omega(i)(j) = omegaMex(i)[j];
			alpha(i)(j) = alphaMex(i)[j];
			beta(i)(j) = betaMex(i)[j];
			gamma(i)(j) = gammaMex(i)[j];
			xiHat(i)(j) = xiHatMex(i)[j];
			dfC(i)(j) = dfCMex(i)[j];
			dfM(i)(j) = dfMMex(i)[j];
		}
	}

	// Get hist
	vector<double*> histMex(M);
	vector<matrix<double>> hist(M);
	for (size_t i = 0; i < M; i++) {
		histMex(i) = mxGetPr(mxGetField(prhs[6], 0, mxGetFieldNameByNumber(prhs[6], i)));
		hist(i) = matrix<double>(m, n);
	}
	for (size_t i = 0; i < M; i++) {
		for (size_t j = 0; j < m; j++) {
			for (size_t k = 0; k < n; k++) {
				hist(i)(j, k) = histMex(i)[j + k * m];
			}
		}
	}

	// Get marginal, copula, varRedType, kappa
	vector<std::string> marginal(M);
	vector<std::string> copula(M);
	vector<std::string> varRedType(M);
	double* kappaMex;
	vector<double> kappa(M);
	kappaMex = mxGetPr(prhs[14]);
	for (size_t i = 0; i < M; i++) {
		marginal(i) = mxArrayToString(mxGetField(prhs[7], 0, mxGetFieldNameByNumber(prhs[7], i)));
		copula(i) = mxArrayToString(mxGetField(prhs[8], 0, mxGetFieldNameByNumber(prhs[8], i)));
		varRedType(i) = mxArrayToString(mxGetField(prhs[9], 0, mxGetFieldNameByNumber(prhs[9], i)));
		kappa(i) = kappaMex[i];
	}

	// Get d, N, problemType
	size_t N;
	size_t d;
	d = mxGetScalar(prhs[10]);
	N = mxGetScalar(prhs[11]);

	// Get fRes (create the output matrix and pointer)
	vector<double*> fResMex(M);
	vector<matrix<double>> fRes(M, matrix<double>(n, N));
	for (size_t i = 0; i < M; i++) {
		plhs[i] = mxCreateDoubleMatrix(n, N, mxREAL);
		fResMex(i) = mxGetPr(plhs[i]);
	}

	/*
	double* randomMex;
	plhs[2] = mxCreateDoubleMatrix(6, 2000, mxREAL);
	randomMex = mxGetPr(plhs[2]);
	*/

	/* call the computational routine */
	MultipleYieldSim::simMultipleFull(E, rho, mu, omega, alpha, beta, hist,
		marginal, copula, varRedType, d, N, fRes, gamma, kappa,
		xiHat, dfC, dfM);

	/*Convert result*/
	for (size_t i = 0; i < M; i++) {
		for (size_t j = 0; j < n; j++) {
			for (size_t k = 0; k < N; k++) {
				fResMex(i)[j + k * n] = fRes(i)(j, k);
			}
		}
	}
	/*
	for (size_t i = 0; i < k(0); i++) {
		for (size_t j = 0; j < N; j++) {
			randomMex[i + j * k(0)] = test(i, j);
		}
	}
	*/
}
