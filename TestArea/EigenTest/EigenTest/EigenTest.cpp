// EigenTest.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include <chrono>

#include <Eigen/Dense>
#include <Eigen/SVD>

using namespace Eigen;

// https://phylogeny.uconn.edu/tutorial-v2/part-1-ide-project-v2/setting-up-the-eigen-library-v2/
int main() {
	MatrixXd mat = MatrixXd::Random(5, 5);

	std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

	HouseholderQR<MatrixXd> qr(mat.rows(), mat.cols());
	qr.compute(mat);

	std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
	std::cout << "Time difference 1 = " << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << "[ms]" << std::endl;
	begin = std::chrono::steady_clock::now();
	
	//BDCSVD<MatrixXd> svd(mat.transpose() * mat, Eigen::ComputeThinU | Eigen::ComputeThinV);

	end = std::chrono::steady_clock::now();
	//std::cout << svd.singularValues()(127) << std::endl;
	std::cout << "Time difference 2 = " << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << "[ms]" << std::endl;


	//std::cout << qr.rows() << "x" << qr.cols() << std::endl;
}

// Run program: Ctrl + F5 or Debug > Start Without Debugging menu
// Debug program: F5 or Debug > Start Debugging menu

// Tips for Getting Started: 
//   1. Use the Solution Explorer window to add/manage files
//   2. Use the Team Explorer window to connect to source control
//   3. Use the Output window to see build output and other messages
//   4. Use the Error List window to view errors
//   5. Go to Project > Add New Item to create new code files, or Project > Add Existing Item to add existing code files to the project
//   6. In the future, to open this project again, go to File > Open > Project and select the .sln file
