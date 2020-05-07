#include "pch.h"
#include "mex.h"
#include "MultipleYieldSim.h"
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/utility/value_init.hpp>
#include <boost/assign.hpp>
using namespace boost::numeric::ublas;


void simConverter(double* EZeroMex, double* ETauMex, double* fZeroMex, double* piMex, double* kappaMex, double* xiHatMex, int d, double* fZeroResMex, double* fTauResMex, int m, int k, int N) {
	

	/*Create new variables corresponding to the funtion that is being tested*/
	vector<matrix<double>> E(2, matrix<double>(m, k));
	
	vector<double> fZero(m);
	matrix<double> pi(m, 1);
	vector<double> kappa(k);
	matrix<double> xiHat(k, 2);
	vector<matrix<double>> fRes(2, matrix<double>(m, N));


	
	/*Convert input variables*/
	for (int i = 0; i < m; i++) {
		fZero(i) = fZeroMex[i];
		pi(i, 0) = piMex[i];
	}

	for (int i = 0; i < m; i++) {		
		for (int j = 0; j < k; j++) {
			E(0)(i, j) = EZeroMex[i + j*m];
			E(1)(i, j) = ETauMex[i + j*m];
		}
	}

	for (int i = 0; i < k; i++) {
		for (int j = 0; j < 2; j++) {
			xiHat(i, j) = xiHatMex[j + i * 2];
		}
	}

	for (int i = 0; i < k; i++) {
		kappa(i) = kappaMex[i];
	}
	

	
	/*Call function*/
	MultipleYieldSim::simMultipleFull(E, fZero, pi, kappa, xiHat, d, N, fRes);

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
	double *EZero;			            
	double *ETau;
	double *fZero;
	double *pi;
	double *kappa;
	double *xiHat;

	double* fZeroRes;
	double* fTauRes;
	
	int m = mxGetM(prhs[2]);
	int k = mxGetN(prhs[2]);

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

	/* create the output matrix */	
	plhs[0] = mxCreateDoubleMatrix(m, N, mxREAL);
	plhs[1] = mxCreateDoubleMatrix(m, N, mxREAL);

	/* get a pointer to the real data in the output matrix */
	fZeroRes = mxGetPr(plhs[0]);
	fTauRes = mxGetPr(plhs[1]);
	
	
	/* call the computational routine */
	simConverter(EZero, ETau, fZero, pi, kappa, xiHat, d, fZeroRes, fTauRes, m, k, N);

}