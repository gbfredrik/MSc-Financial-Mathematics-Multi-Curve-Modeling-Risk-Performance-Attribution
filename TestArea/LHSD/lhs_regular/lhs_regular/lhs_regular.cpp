#include <iostream>

//Boost packages for numeric represenation
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/io.hpp>

//Boost packages for maths
#include <boost/math/tools/bivariate_statistics.hpp>

//Package to generate random variables
#include <random>

using namespace boost::numeric::ublas;

matrix<double> GC_Sim(matrix<double> corr, int n);
vector<double> gen_normal(double m, double s, int n);
matrix<double> corr(matrix<double> input);
matrix<double> gen_test(int rows, int cols);

int main() {
    
    //gen_normal(0, 1, 10);

	//Generate test matrix
	matrix<double> test(3, 10);
	test = gen_test(3, 10);
	
	//Calculate correlation matrix
	matrix<double> corr(3, 10);
	corr = corr(test);
	
	//GC_Sim(corr(test), 20);
	
    
}


matrix<double> GC_Sim(matrix<double> corr, int n) {
    int dim = corr.size1(); // Number of risk factors
    matrix<double> U(dim, n); // Return matrix
    
    //Generate a matrix with random numbers
    matrix<double> y(dim, n);

    //for (auto values_iter = gen_normal(0,1,n).begin(); values_iter != gen_normal(0,1,n).end(); ++values_iter, ++iter1) 
    
    //std::copy(gen_normal(0, 1, n).begin(), gen_normal(0, 1, n).end(), y.begin1());
    
	/*
	for (int i = 0; i < n; ++i) {
		std::cout << y << std::endl;
	}
	*/
    return U;
}


//Generate n normal variables with mean m and standard deviation s
vector<double> gen_normal(double m, double s, int n) {
    std::default_random_engine generator;
    std::normal_distribution<double> distribution(m, s);
    vector<double> rand(n);
        
    for (int i = 0; i < n; ++i) {
        double number = distribution(generator);
        rand[i] = number;
        //std::cout << rand[i] << std::endl;
    }
    
    return rand;
}

matrix<double> corr(matrix<double> input) {
	matrix<double> corr;



	return corr;
}

matrix<double> gen_test(int rows, int cols) {
	std::random_device rd;
	std::mt19937 mt(rd());
	std::uniform_real_distribution<double> dist(1.0, 10.0);
	matrix<double> test(3, 10);
	for (int i = 0; i < test.size1; i++) {
		for (int j = 0; j < test.size2; j++) {
			test(i, j) = dist(mt);
		}
	}

	return test;
}


