#include "Gaussian_Copula.h"


Gaussian_Copula::Gaussian_Copula(matrix<double> series) : Distribution(series) {
	time_series = series;
}

void Gaussian_Copula::getSeries() {
	std::cout << time_series << "\n";
}

double Gaussian_Copula::function_value(vector<double> const& x) {
	std::cout << "in FN Gaussian copula" << "\n\n";
	double sum = 0;
	int n = time_series.size2(); //Number of riskfaktors
	matrix<double> P = buildP(x);
	std::cout << "P = " << P << "\n\n";
	double det_P = matrixOperations::ublasToMatrixXd(P).determinant();
	matrix<double> P_inv = matrixOperations::matrixXdToUblas(matrixOperations::ublasToMatrixXd(P).inverse());

	for (size_t i = 0; i < time_series.size1(); ++i) {

		matrix_row<matrix<double> > U_row(time_series, i);

		//Get N_inv(U)
		vector<double> N_inv(n);
		for (size_t j = 0; j < N_inv.size(); ++j) {
			boost::math::normal norm = boost::math::normal::normal_distribution();
			N_inv(j) = quantile(norm, U_row(j));
		}

		//Help products and summations
		identity_matrix<double> I(n);
		vector<double> nPI = prod(N_inv, P_inv-I);
		double nPIn = inner_prod(nPI, N_inv);
				
		sum = sum - 0.5 * nPIn - 0.5 * log(det_P);
	}

	return -sum;
}


vector<double> Gaussian_Copula::calcGradients(vector<double> const& x) {
	vector<double> gradients(x.size());
	int n = time_series.size2(); //Number of riskfaktors
	double sum = 0;
	zero_vector<double> zeroVec(n * n);
	vector<double> dFdP = zeroVec;
	double dFdnu = 0;

	//Get rho gradients as a vector
	for (size_t i = 0; i < time_series.size1(); ++i) {

		//Get T_inv(U)
		matrix_row<matrix<double> > U_row(time_series, i);
		//std::cout << "U_row = " << U_row << "\n\n";
		vector<double> N_inv(n);
		for (size_t j = 0; j < N_inv.size(); ++j) {
			boost::math::normal norm = boost::math::normal::normal_distribution();
			N_inv(j) = quantile(norm, U_row(j));
		}
		//std::cout << "T_inv = " << T_inv << "\n\n";

		matrix<double> P = buildP(x);
		//std::cout << "P = " << P << "\n\n";
		matrix<double> P_inv = matrixOperations::matrixXdToUblas(matrixOperations::ublasToMatrixXd(P).inverse());

		//Get P_inv as a vector of the columns
		vector<double> vec_Pinv = matrixToVector(P_inv);

		//Kroneckers product of two vectors
		vector<double> vKron = kronOfVectors(prod(N_inv, P_inv), prod(N_inv, P_inv));

		//Get gradients for time_series(i)
		vector<double> dFdP_temp = 0.5 * vKron - 0.5 * vec_Pinv;
		dFdP = dFdP + dFdP_temp;
	}

	std::cout << "dFdP = " << dFdP << "\n\n";

	//Get rho gradients as matrix
	matrix<double> dFdP_mat = vectorToMatrix(dFdP);
	std::cout << "dFdP_mat = " << dFdP_mat << "\n\n";

	//Get optimization parameters from rho matrix
	vector<double> dfdParams = getElements(dFdP_mat);

	return -dfdParams;
}


matrix<double> Gaussian_Copula::vectorToMatrix(vector<double> const& vec) {
	int n = time_series.size2();
	matrix<double> resMatrix(n, n);
	int counter = 0;

	for (size_t j = 0; j < n; ++j) {
		for (size_t i = 0; i < n; ++i) {
			resMatrix(i, j) = vec(counter);
			counter = counter + 1;
		}
	}

	return resMatrix;
}

matrix<double> Gaussian_Copula::buildP(vector<double> const& x) {
	int n = time_series.size2();
	matrix<double> P(n, n);
	int counter = 0;	//keep track of fetched elements

	for (size_t i = 0; i < n; ++i) {
		for (size_t j = 0; j < i + 1; ++j) {
			if (i == j) {
				P(i, j) = 1;
			}
			else {
				P(i, j) = x(counter);
				P(j, i) = x(counter);
				counter = counter + 1;
			}
		}
	}
	return P;
}


vector<double> Gaussian_Copula::getElements(matrix<double> const& P) {
	int n = time_series.size2();
	vector<double> resVec((n - 1) * n * 0.5); // Vector with optimization parameters in P
	int counter = 0;	//keep track of fetched elements

	for (size_t i = 0; i < n; ++i) {
		for (size_t j = 0; j < i; ++j) {
			resVec(counter) = P(i, j);
			counter = counter + 1;
		}
	}
	return resVec;
}

vector<double> Gaussian_Copula::kronOfVectors(vector<double> const& v1, vector<double> const& v2) {
	vector<double> vKron(v1.size() * v2.size());

	int counter = 0;

	for (size_t i = 0; i < v1.size(); ++i) {
		for (size_t j = 0; j < v2.size(); ++j) {
			vKron(counter) = v1(i) * v2(j);
			counter = counter + 1;
		}
	}

	return vKron;
}

vector<double> Gaussian_Copula::matrixToVector(matrix<double> const& matrix) {
	int n = matrix.size1();
	vector<double> resVec(n * n);
	int counter = 0;

	for (size_t j = 0; j < n; ++j) {
		for (size_t i = 0; i < n; ++i) {
			resVec(counter) = matrix(i, j);
			counter = counter + 1;
		}
	}
	return resVec;
}

vector<double> Gaussian_Copula::calcNumGradients(vector<double> const& x) {

	double epsilon = 2.2 * pow(10, -16);
	vector<double> increment(sqrt(epsilon) * x);
	vector<double> num_gradients(x.size());

	/*
	vector<double> x_0diff = x;
	x_0diff(0) += increment(0);
	vector<double> x_1diff = x;
	x_1diff(1) += increment(1);
	vector<double> x_2diff = x;
	x_2diff(2) += increment(2);
	vector<double> x_3diff = x;
	x_3diff(3) += increment(3);
	vector<double> x_4diff = x;
	x_4diff(4) += increment(4);

	num_gradients(0) = (function_value(x_0diff) - function_value(x)) / (x_0diff(0) - x(0));
	num_gradients(1) = (function_value(x_1diff) - function_value(x)) / (x_1diff(1) - x(1));
	num_gradients(2) = (function_value(x_2diff) - function_value(x)) / (x_2diff(2) - x(2));
	num_gradients(3) = (function_value(x_3diff) - function_value(x)) / (x_3diff(3) - x(3));
	num_gradients(4) = (function_value(x_4diff) - function_value(x)) / (x_4diff(4) - x(4));

	*/

	return num_gradients;
}

double Gaussian_Copula::calcStepSize(vector<double> const& x, vector<double> const& d) {

	double a = 1;
	//bool accepted = false;
	////Kontrollera att rho är positiv definit, dvs minsta egenvärdet är positivt, annars halvera steglängden.

	//while (!accepted) {
	//	accepted = true;

	//	for (size_t i = 0; i < x.size() - 1; ++i) {
	//		if (x(i) + a * d(i) < 0 || x(i) + a * d(i) > 1) {
	//			accepted = false;
	//			a = 0.5 * a;
	//			break;
	//		}
	//	}

	//	if (x(x.size() - 1) + a * d(d.size() - 1) <= 2) {
	//		accepted = false;
	//		a = 0.5 * a;
	//	}

	//	matrix<double> Pnext = buildP(x + a * d);
	//	double minEigenvalue = ();

	//	if (minEigenvalue <= 0) {
	//		accepted = false;
	//	}
	//}
	return a;
}
