#include "mex.h"
#include "MultipleYieldSim.h"
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/utility/value_init.hpp>
#include <boost/assign.hpp>
using namespace boost::numeric::ublas;


void simConverter(double* EZeroMex, double* ETauMex, double* fZeroMex, double* piMex, double* kappaMex, double* xiHatMex, int d, double* fZeroResMex, double* fTauResMex, int m, int n, int N, int k) {
	
	/*Create new variables corresponding to the funtion that is being tested*/
		vector<matrix<double>> E(2, matrix<double>(m, n));
	
	vector<double> fZero(m);
	matrix<double> pi(m, 1);
	vector<double> kappa(k);
	matrix<double> xiHat(k, 2);
	vector<matrix<double>> fRes(2, matrix<double>(m, N));


	
	/*Convert input variables*/
	for (int i = 0; i < m; i++) {
		fZero(i) = fZeroMex[i];
		pi(i, 1) = piMex[i];
	}

	/*
	for (int i = 0; i < m; i++) {		
		for (int j = 0; j < n; j++) {
			E(0)(i, j) = EZeroMex[j + i*n];
			E(1)(i, j) = ETauMex[j + i*n];
		}
	}
	
	mexPrintf("1");
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < 2; i++) {
			xiHat(i, j) = xiHatMex[j + i * 2];
		}
	}
	mexPrintf("2");
	
	mexPrintf("3");
	for (int i = 0; i < n; i++) {
		kappa(i) = kappaMex[i];
	}
	mexPrintf("4");
	*/

	
	/*Call function*/
	mexPrintf("hehe");
	MultipleYieldSim::simMultipleFull(E, fZero, pi, kappa, xiHat, d, fRes);
	mexPrintf("japp");

	/*Convert result*/
	for (int i = 0; i < m; i++) {
		for (int j = 0; j < N; j++) {
			fZeroResMex[j + i*N] = fRes(0)(i, j);
			fTauResMex[j + i*N] = fRes(1)(i, j);
		}
	}
	mexPrintf("yeboii");

}


/* The gateway function */
void mexFunction(int nlhs, mxArray* plhs[],
	int nrhs, const mxArray* prhs[]) {

	int N = 1;
	int k = 6;

	int d = 0;
	double *EZero;			            
	double *ETau;
	double *fZero;
	double *pi;
	double *kappa;
	double *xiHat;

	double* fZeroRes;
	double* fTauRes;
	
	int m = mxGetM(prhs[1]);
	int n = mxGetN(prhs[1]);

	/* get the values of input scalars  */
	d = mxGetScalar(prhs[0]);

	/* get the values of input matrices  */
	EZero = mxGetPr(prhs[1]);
	ETau = mxGetPr(prhs[2]);
	fZero = mxGetPr(prhs[3]);
	pi = mxGetPr(prhs[4]);
	kappa = mxGetPr(prhs[5]);
	xiHat = mxGetPr(prhs[6]);

	/* create the output matrix */	
	plhs[0] = mxCreateDoubleMatrix(m, N, mxREAL);
	plhs[1] = mxCreateDoubleMatrix(m, N, mxREAL);

	/* get a pointer to the real data in the output matrix */
	fZeroRes = mxGetPr(plhs[0]);
	fTauRes = mxGetPr(plhs[1]);
	
	/* call the computational routine */
	simConverter(EZero, ETau, fZero, pi, kappa, xiHat, d, fZeroRes, fTauRes, m, n, N, k);

}