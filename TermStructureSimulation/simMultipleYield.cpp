#include "pch.h"
#include "mex.h"
#include "MultipleYieldSim.h"
#include <boost/numeric/ublas/matrix.hpp>
#include <string>

using namespace boost::numeric::ublas;


void simConverter(double* EZeroMex, double* ETauMex, double* fZeroMex, double* piMex, 
	double* kappaMex, double* xiHatMex, int d, double* fZeroResMex, double* 
	fTauResMex, int m, int kZero, int kTau, int N, double* omegaMex, double* alphaMex,
	double* betaMex, double* fHistMex, double* piHistMex, int days,
	std::string marginalZeroMex, std::string marginalTauMex, std::string copulaZeroMex, std::string copulaTauMex,
	std::string varRedTypeZeroMex, std::string varRedTypeTauMex, double* muMex, double* dfMex) {
	

	/*Create new variables corresponding to the funtion that is being tested*/
	vector<matrix<double>> E(2);
	E(0) = matrix<double>(m, kZero);
	E(1) = matrix<double>(m, kTau);


	vector<double> fZero(m);
	matrix<double> pi(m, 1);
	vector<double> kappa(kZero);
	matrix<double> xiHat(kZero, 2);
	vector<std::string> marginal(2);
	vector<std::string> copula(2);
	vector<std::string> varRedType(2);
	vector<matrix<double>> fRes(2, matrix<double>(m, N));
	vector<matrix<double>> hist(2, matrix<double>(days, m));

	vector<double> omega(kZero);
	vector<double> beta(kZero);
	vector<double> alpha(kZero);
	matrix<double> mu(kZero, 2);
	matrix<double> df(kZero, 2);

	/*Convert input variables*/
	for (int i = 0; i < m; i++) {
		fZero(i) = fZeroMex[i];
		pi(i, 0) = piMex[i];
	}

	for (int i = 0; i < m; i++) {		
		for (int j = 0; j < kZero; j++) {
			E(0)(i, j) = EZeroMex[i + j*m];
		}
	}
	for (int i = 0; i < m; i++) {
		for (int j = 0; j < kTau; j++) {
			E(1)(i, j) = ETauMex[i + j * m];
		}
	}
	
	for (int i = 0; i < kZero; i++) {
		for (int j = 0; j < 2; j++) {
			xiHat(i, j) = xiHatMex[i + j * kZero];
			mu(i, j) = muMex[i + j * kZero];
			df(i, j) = dfMex[i + j * kZero];
		}
	}

	for (int i = 0; i < kZero; i++) {
		kappa(i) = kappaMex[i];
		omega(i) = omegaMex[i];
		alpha(i) = alphaMex[i];
		beta(i) = betaMex[i];
	}

	for (int j = 0; j < days; j++) {
		for (int k = 0; k < m; k++) {
			hist(0)(j, k) = fHistMex[j + k * days];
			hist(1)(j, k) = piHistMex[j + k * days];
		}
	}

	marginal(0) = marginalZeroMex;
	marginal(1) = marginalZeroMex;
	copula(0) = copulaZeroMex;
	copula(1) = copulaTauMex;
	varRedType(0) = varRedTypeZeroMex;
	varRedType(1) = varRedTypeTauMex;


	/*Call function*/
	MultipleYieldSim::simMultipleFull(E, fZero, pi, kappa, xiHat, omega, 
		alpha, beta, hist, marginal, copula, varRedType, 
		mu, df, d, N, fRes);

	/*Convert result*/
	for (int i = 0; i < m; i++) {
		for (int j = 0; j < N; j++) {
			fZeroResMex[i + j * m] = fRes(0)(i, j);
			fTauResMex[i + j * m] = fRes(1)(i, j);
		}
	} 
}


/* The gateway function */
void mexFunction(int nlhs, mxArray* plhs[],
	int nrhs, const mxArray* prhs[]) {

	int N = 0;

	int d = 0;
	double* EZero;			            
	double* ETau;
	double* fZero;
	double* pi;
	double* kappa;
	double* xiHat;
	double* omega;
	double* alpha;
	double* beta;
	double* fHist;
	double* piHist;
	double* mu;
	double* df;


	std::string marginalZero;
	std::string copulaZero;
	std::string varRedTypeZero;

	std::string marginalTau;
	std::string copulaTau;
	std::string varRedTypeTau;

	double *fZeroRes;
	double *fTauRes;
	
	int days = mxGetN(prhs[11]);
	int m = mxGetM(prhs[2]);
	int kZero = mxGetN(prhs[2]);
	int kTau = mxGetN(prhs[3]);


	/* get the values of input scalars  */
	d = mxGetScalar(prhs[0]);
	N = mxGetScalar(prhs[1]);

	/* get the values of input matrices  */
	EZero = mxGetPr(prhs[2]);
	ETau = mxGetPr(prhs[3]);
	fZero = mxGetPr(prhs[4]);
	pi = mxGetPr(prhs[5]);
	kappa = mxGetPr(prhs[6]);
	xiHat = mxGetPr(prhs[7]);
	omega = mxGetPr(prhs[8]);
	alpha = mxGetPr(prhs[9]);
	beta = mxGetPr(prhs[10]);
	fHist = mxGetPr(prhs[11]);
	piHist = mxGetPr(prhs[12]);
	mu = mxGetPr(prhs[13]);
	df = mxGetPr(prhs[14]);

	marginalZero = mxArrayToString(prhs[15]);
	marginalTau = mxArrayToString(prhs[16]);
	copulaZero = mxArrayToString(prhs[17]);
	copulaTau = mxArrayToString(prhs[18]);
	varRedTypeZero = mxArrayToString(prhs[19]);
	varRedTypeTau = mxArrayToString(prhs[20]);

	/* create the output matrix */	
	plhs[0] = mxCreateDoubleMatrix(m, N, mxREAL);
	plhs[1] = mxCreateDoubleMatrix(m, N, mxREAL);

	/* get a pointer to the real data in the output matrix */
	fZeroRes = mxGetPr(plhs[0]);
	fTauRes = mxGetPr(plhs[1]);
	
	
	/* call the computational routine */
	simConverter(EZero, ETau, fZero, pi, kappa, xiHat, d, fZeroRes, 
		fTauRes, m, kZero, kTau, N, omega, alpha, beta, fHist, piHist, days,
		marginalZero, marginalTau, copulaZero, copulaTau,
		varRedTypeZero, varRedTypeTau, mu, df);

}

