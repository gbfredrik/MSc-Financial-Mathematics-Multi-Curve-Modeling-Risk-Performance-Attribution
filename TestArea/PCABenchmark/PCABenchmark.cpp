#include "bench/BenchTimer.h"
#include <Eigen/Dense>

#include <iostream>
#include <map>
#include <vector>
#include <string>
#include <sstream>

using namespace Eigen;

std::map<std::string, Array<double, 1, 8, DontAlign|RowMajor>> results;
std::vector<std::string> labels;
std::vector<Array2i> sizes;

template<typename Solver, typename MatrixType>
EIGEN_DONT_INLINE
void compute_norm_equation(Solver& solver, MatrixType const& A) {
	if (A.rows() != A.cols()) {
		solver.compute(A.transpose() * A);
	} else {
		solver.compute(A);
	}
}

template<typename Solver, typename MatrixType>
EIGEN_DONT_INLINE
void compute(Solver& solver, MatrixType const& A) {
	solver.compute(A);
}

template<typename Scalar, int Size>
void bench(int id, int rows, int size = Size) {
	typedef Matrix<Scalar, Dynamic, Size> Mat;
	typedef Matrix<Scalar, Dynamic, Dynamic> MatDyn;
	typedef Matrix<Scalar, Size, Size> MatSquare;
	Mat A(rows, size);
	A.setRandom();

	if (rows == size) {
		A = A * A.adjoint();
	}

	BenchTimer t_llt, t_ldlt, t_lu, t_fplu, t_qr, t_cpqr, t_cod, t_fpqr, t_jsvd, t_bdcsvd;

	int svd_opt = ComputeThinU|ComputeThinV;
  
	int tries = 3;
	int rep = 1000 / size;
	if (rep == 0) {
		rep = 1;
	}
	//   rep = rep*rep;
  
	LLT<MatSquare> llt(size);
	LDLT<MatSquare> ldlt(size);
	PartialPivLU<MatSquare> lu(size);
	FullPivLU<MatSquare> fplu(size, size);
	HouseholderQR<Mat> qr(A.rows(), A.cols());
	ColPivHouseholderQR<Mat> cpqr(A.rows(), A.cols());
	CompleteOrthogonalDecomposition<Mat> cod(A.rows(), A.cols());
	FullPivHouseholderQR<Mat> fpqr(A.rows(), A.cols());
	JacobiSVD<MatDyn> jsvd(A.rows(), A.cols());
	BDCSVD<MatDyn> bdcsvd(A.rows(), A.cols());
  
	BENCH(t_llt, tries, rep, compute_norm_equation(llt, A));
	BENCH(t_ldlt, tries, rep, compute_norm_equation(ldlt, A));
	BENCH(t_lu, tries, rep, compute_norm_equation(lu, A));
	if (size <= 1000) {
		BENCH(t_fplu, tries, rep, compute_norm_equation(fplu, A));
	}
	BENCH(t_qr, tries, rep, compute(qr, A));
	BENCH(t_cpqr, tries, rep, compute(cpqr, A));
	BENCH(t_cod, tries, rep, compute(cod, A));
	if (size * rows <= 10000000) {
		BENCH(t_fpqr, tries, rep, compute(fpqr, A));
	}
	if (size < 500) {// JacobiSVD is really too slow for too large matrices
		BENCH(t_jsvd, tries, rep, jsvd.compute(A, svd_opt));
	}
	//   if(size*rows<=20000000)
	BENCH(t_bdcsvd, tries, rep, bdcsvd.compute(A, svd_opt));
  
	results["LLT"][id] = t_llt.best();
	results["LDLT"][id] = t_ldlt.best();
	results["PartialPivLU"][id] = t_lu.best();
	results["FullPivLU"][id] = t_fplu.best();
	results["HouseholderQR"][id] = t_qr.best();
	results["ColPivHouseholderQR"][id] = t_cpqr.best();
	results["CompleteOrthogonalDecomposition"][id] = t_cod.best();
	results["FullPivHouseholderQR"][id] = t_fpqr.best();
	results["JacobiSVD"][id] = t_jsvd.best();
	results["BDCSVD"][id] = t_bdcsvd.best();
}

int main() {
	labels.push_back("LLT");
	labels.push_back("LDLT");
	labels.push_back("PartialPivLU");
	labels.push_back("FullPivLU");
	labels.push_back("HouseholderQR");
	labels.push_back("ColPivHouseholderQR");
	labels.push_back("CompleteOrthogonalDecomposition");
	labels.push_back("FullPivHouseholderQR");
	labels.push_back("JacobiSVD");
	labels.push_back("BDCSVD");

	for (size_t i{ 0 }; i < labels.size(); ++i) {
		results[labels[i]].fill(-1);
	}

	const int small = 8;
	//sizes.push_back(Array2i(small,small));
	//sizes.push_back(Array2i(100,100));
	sizes.push_back(Array2i(1500,1500));
	//sizes.push_back(Array2i(4000,4000));
	//sizes.push_back(Array2i(10000,small));
	//sizes.push_back(Array2i(10000,100));
	sizes.push_back(Array2i(10950,1500));
	//sizes.push_back(Array2i(10000,4000));

	using namespace std;

	for (size_t k{ 0 }; k < sizes.size(); ++k) {
		cout << sizes[k](0) << "x" << sizes[k](1) << "...\n";
		bench<double, Dynamic>(k, sizes[k](0), sizes[k](1));
	}
	
	cout.width(32);
	cout << "solver/size";
	cout << "  ";
	for (size_t k{ 0 }; k < sizes.size(); ++k) {
		std::stringstream ss;
		ss << sizes[k](0) << "x" << sizes[k](1);
		cout.width(10);
		cout << ss.str();
		cout << " ";
	}
	cout << endl;

	for (size_t i{ 0 }; i < labels.size(); ++i) {
		cout.width(32);
		cout << labels[i];
		cout << "  ";
		ArrayXd r = (results[labels[i]] * 100000.f).floor() / 100.f;
		for (size_t k{ 0 }; k < sizes.size(); ++k) {
			cout.width(10);
			if (r(k) >= 1e6) {
				cout << "-";
			} else {
				cout << r(k);
			}
			cout << " ";
		}
		cout << endl;
	}

	// HTML output
	/*cout << "<table class=\"manual\">" << endl;
	cout << "<tr><th>solver/size</th>" << endl;
	for(int k=0; k<sizes.size(); ++k)
	cout << "  <th>" << sizes[k](0) << "x" << sizes[k](1) << "</th>";
	cout << "</tr>" << endl;
	for(int i=0; i<labels.size(); ++i)
	{
	cout << "<tr";
	if(i%2==1) cout << " class=\"alt\"";
	cout << "><td>" << labels[i] << "</td>";
	ArrayXd r = (results[labels[i]]*100000.f).floor()/100.f;
	for(int k=0; k<sizes.size(); ++k)
	{
		if(r(k)>=1e6) cout << "<td>-</td>";
		else
		{
		cout << "<td>" << r(k);
		if(i>0)
			cout << " (x" << numext::round(10.f*results[labels[i]](k)/results["LLT"](k))/10.f << ")";
		if(i<4 && sizes[k](0)!=sizes[k](1))
			cout << " <sup><a href=\"#note_ls\">*</a></sup>";
		cout << "</td>";
		}
	}
	cout << "</tr>" << endl;
	}
	cout << "</table>" << endl;*/

	//cout << "llt                             (ms)  " << (results["llt"]*1000.).format(fmt) << "\n";
	//cout << "ldlt                             (%)  " << (results["ldlt"]/results["llt"]).format(fmt) << "\n";
	//cout << "partialpivlu                     (%)  " << (results["partialpivlu"]/results["llt"]).format(fmt) << "\n";
	//cout << "fullpivlu                        (%)  " << (results["fullpivlu"]/results["llt"]).format(fmt) << "\n";
	//cout << "householderqr                    (%)  " << (results["householderqr"]/results["llt"]).format(fmt) << "\n";
	//cout << "colpivhouseholderqr              (%)  " << (results["colpivhouseholderqr"]/results["llt"]).format(fmt) << "\n";
	//cout << "completeorthogonaldecomposition  (%)  " << (results["completeorthogonaldecomposition"]/results["llt"]).format(fmt) << "\n";
	//cout << "fullpivhouseholderqr             (%)  " << (results["fullpivhouseholderqr"]/results["llt"]).format(fmt) << "\n";
	//cout << "jacobisvd                        (%)  " << (results["jacobisvd"]/results["llt"]).format(fmt) << "\n";
	//cout << "bdcsvd                           (%)  " << (results["bdcsvd"]/results["llt"]).format(fmt) << "\n";
}
