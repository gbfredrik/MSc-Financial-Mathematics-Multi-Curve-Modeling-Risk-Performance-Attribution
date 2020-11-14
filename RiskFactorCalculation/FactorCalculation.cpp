#include "FactorCalculation.h"

#include "../MathLibrary/matrixOperations.h"

#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/QR>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <Spectra/MatOp/DenseSymMatProd.h>
#include <Spectra/SymEigsSolver.h>
#include <rsvd/Prelude.hpp>

#include <iostream>
#include <Windows.h>
#include <WinUser.h>

#include "sample_handler.h"

using namespace boost::numeric;

bool FactorCalculation::iram(
    ublas::matrix<double> const& input, 
    int const k, 
    ublas::matrix<double>& m_E, 
    ublas::vector<double>& v_Lambda
) {
    // Utilizes the SymEigsSolver (IRAM) from the Spectra library
    using namespace Spectra;

    // Transform ublas matrix into Eigen::matrix
    Eigen::MatrixXd D{ matrixOperations::ublasToMatrixXd(input) };
    
    if (D.rows() < D.cols()) {
        Eigen::HouseholderQR<Eigen::MatrixXd> qr{ D.transpose() };
        // Can use FullPivHouseholderQR for greater precision
        //Eigen::FullPivHouseholderQR<Eigen::MatrixXd> qr(D.transpose());

        // Define the positive definite RTR matrix
        Eigen::MatrixXd thinQ{ qr.householderQ() * Eigen::MatrixXd::Identity(D.cols(), D.rows()) }; // qr.householderQ()
        Eigen::MatrixXd R_temp{ qr.matrixQR().triangularView<Eigen::Upper>() };
        Eigen::MatrixXd thinR{ R_temp.topRows(R_temp.cols()) };
        Eigen::MatrixXd RRT{ thinR * thinR.transpose() };

        // Construct matrix operation objects
        DenseSymMatProd<double> op{ RRT }; // RRT

        // Construct eigen solver object, requesting the largest k eigenvalues in magnitude
        SymEigsSolver<double, LARGEST_ALGE, DenseSymMatProd<double>> eigs(&op, k, 2 * k + 1);

        // Initialize and compute
        eigs.init();
        int nconv{ eigs.compute() };

        // Retrieve results
        if (eigs.info() == SUCCESSFUL) {
            // Return eigenpairs
            m_E = matrixOperations::matrixXdToUblas(thinQ * eigs.eigenvectors());
            v_Lambda = matrixOperations::vectorXdToUblas(eigs.eigenvalues()); // OK
        }
    } else {
        // Define the positive definite C matrix (covariance)
        Eigen::MatrixXd C{ D.transpose() * D };

        // Construct matrix operation objects
        DenseSymMatProd<double> op(C);

        // Construct eigen solver object, requesting the largest k eigenvalues in magnitude
        SymEigsSolver<double, LARGEST_ALGE, DenseSymMatProd<double>> eigs(&op, k, 2 * k + 1);

        eigs.init();
        int nconv{ eigs.compute() };

        // Retrieve results
        if (eigs.info() == SUCCESSFUL) {
            // Return eigenpairs
            m_E = matrixOperations::matrixXdToUblas(eigs.eigenvectors());
            v_Lambda = matrixOperations::vectorXdToUblas(eigs.eigenvalues()); // OK
        }
    }
    
    return v_Lambda.size() > 0;
}

bool FactorCalculation::eigen_bdcsvd(
    ublas::matrix<double> const& input, 
    int const k, 
    ublas::matrix<double>& m_E, 
    ublas::vector<double>& v_Lambda
) {
    // Utilizes BDCSVD from Eigen library
    int svd_opt{ Eigen::ComputeThinU | Eigen::ComputeThinV };

    Eigen::MatrixXd H{ matrixOperations::ublasToMatrixXd(input) };

    if (H.rows() < H.cols()) {
        Eigen::HouseholderQR<Eigen::MatrixXd> qr(H.transpose()); // FullPivHouseholderQR<Matrix<double, Dynamic, Size>> fpqr(A.rows(), A.cols());
        Eigen::MatrixXd thinQ{ qr.householderQ() * Eigen::MatrixXd::Identity(H.cols(), H.rows()) };
        Eigen::MatrixXd R_temp{ qr.matrixQR().triangularView<Eigen::Upper>() };
        Eigen::MatrixXd thinR{ R_temp.topRows(R_temp.cols()) };
        Eigen::MatrixXd RRT{ thinR * thinR.transpose() };

        Eigen::BDCSVD<Eigen::MatrixXd> bdcsvd(RRT, svd_opt);

        // Return eigenpairs
        if (k == 0) {
            m_E = matrixOperations::matrixXdToUblas(thinQ * bdcsvd.matrixV());
            v_Lambda = matrixOperations::vectorXdToUblas(bdcsvd.singularValues()); // OK
        } else {
            m_E = matrixOperations::matrixXdToUblas(thinQ * bdcsvd.matrixV().block(0, 0, bdcsvd.matrixV().rows(), k));
            v_Lambda = matrixOperations::vectorXdToUblas(bdcsvd.singularValues().head(k)); // OK
            // See: https://stackoverflow.com/questions/34373757/piece-wise-square-of-vector-piece-wise-product-of-two-vectors-in-c-eigen
        }
    } else {
        Eigen::BDCSVD<Eigen::MatrixXd> bdcsvd(H.transpose() * H, svd_opt);
        
        if (k == 0) {
            m_E = matrixOperations::matrixXdToUblas(bdcsvd.matrixV());
            v_Lambda = matrixOperations::vectorXdToUblas(bdcsvd.singularValues().array().square()); // OK
        } else {
            m_E = matrixOperations::matrixXdToUblas(bdcsvd.matrixV().block(0, 0, bdcsvd.matrixV().rows(), k));
            v_Lambda = matrixOperations::vectorXdToUblas(bdcsvd.singularValues().head(k).array().square()); // OK
        }
    }
        
    return v_Lambda.size() > 0;
}

bool FactorCalculation::eigen_rsvd(
    ublas::matrix<double> const& input,
    int const k,
    ublas::matrix<double>& m_E,
    ublas::vector<double>& v_Lambda
) {
    // Utilizes Randomized SVD from Eigen extension
    int svd_opt{ Eigen::ComputeThinU | Eigen::ComputeThinV };

    Eigen::MatrixXd H{ matrixOperations::ublasToMatrixXd(input) };
    const Eigen::Index reducedRank{ k };

    // Initialize PRNG for the Eigen random matrix generation
    std::srand(static_cast<unsigned int>(777));
    std::mt19937_64 randomEngine{};
    randomEngine.seed(777);

    if (H.rows() < H.cols()) {
        Eigen::HouseholderQR<Eigen::MatrixXd> qr(H.transpose()); // FullPivHouseholderQR<Matrix<double, Dynamic, Size>> fpqr(A.rows(), A.cols());
        Eigen::MatrixXd thinQ{ qr.householderQ() * Eigen::MatrixXd::Identity(H.cols(), H.rows()) };
        Eigen::MatrixXd R_temp{ qr.matrixQR().triangularView<Eigen::Upper>() };
        Eigen::MatrixXd thinR{ R_temp.topRows(R_temp.cols()) };
        Eigen::MatrixXd RRT{ thinR * thinR.transpose() };

        // Randomized SVD
        Rsvd::RandomizedSvd<Eigen::MatrixXd, std::mt19937_64, Rsvd::SubspaceIterationConditioner::Lu>
            rsvd(randomEngine);
        rsvd.compute(RRT, reducedRank);

        m_E = matrixOperations::matrixXdToUblas(thinQ * rsvd.matrixV());
        v_Lambda = matrixOperations::vectorXdToUblas(rsvd.singularValues());
    } else {
        Rsvd::RandomizedSvd<Eigen::MatrixXd, std::mt19937_64, Rsvd::SubspaceIterationConditioner::Lu>
            rsvd(randomEngine);
        rsvd.compute(H, reducedRank);

        // Retrieve results
        m_E = matrixOperations::matrixXdToUblas(rsvd.matrixV());
        v_Lambda = matrixOperations::vectorXdToUblas(rsvd.singularValues().array().square());
    }

    //const Eigen::MatrixXd rsvdApprox{ rsvd.matrixU() *
    //    rsvd.singularValues().asDiagonal() *
    //    rsvd.matrixV().adjoint() };
    //const auto rsvdErr{ Rsvd::relativeFrobeniusNormError(x, rsvdApprox) };
    //std::cout << "Randomized SVD reconstruction error: " << rsvdErr << std::endl;

    return v_Lambda.size() > 0;
}


ublas::matrix<double> FactorCalculation::compute_risk_factors(
    ublas::matrix<double> const& m_E_k, 
    ublas::matrix<double> const& m_delta_f
) {
    return prod(trans(m_E_k), trans(m_delta_f));
}

double FactorCalculation::smallest_eigval(ublas::matrix<double> const& input) {
    // Utilizes the SymEigsSolver (IRAM) from the Spectra library
    using namespace Spectra;
    int k{ 1 };
    // Transform ublas matrix into Eigen::matrix
    Eigen::MatrixXd C{ matrixOperations::ublasToMatrixXd(input) };

    // Construct matrix operation objects
    DenseSymMatProd<double> op(C);

    // Construct eigen solver object, requesting the largest k eigenvalues in magnitude
    SymEigsSolver<double, SMALLEST_ALGE, DenseSymMatProd<double>> eigs(&op, k, 2 * k + 1);

    eigs.init();
    int nconv{ eigs.compute() };

    // Retrieve results
    if (eigs.info() == SUCCESSFUL) {
        // Return smallest eigenvalue
        return eigs.eigenvalues()[0];
    }
    
    MessageBoxA(
        NULL,
        "Failed to find smallest eigenvalue.",
        "Status",
        MB_OK);
    return -1.0;
}

double FactorCalculation::eig_norm_error(
    ublas::matrix<double> const& m_A, 
    ublas::vector<double> const& v_x, 
    double const lambda
) {
    return norm_2(prod(m_A, v_x) - lambda * v_x);
}

ublas::vector<double> FactorCalculation::eig_all_norm_errors(
    ublas::matrix<double> const& m_A,
    ublas::matrix<double> const& m_x,
    ublas::vector<double> const& v_Lambda
) {
    size_t len{ v_Lambda.size() };
    ublas::vector<double> v_errors(len);

    for (size_t i{ 0 }; i < len; ++i) {
        v_errors(i) = eig_norm_error(m_A, column(m_x, i), v_Lambda(i));
    }

    return v_errors;
}

ublas::matrix<double> FactorCalculation::clean_data(ublas::matrix<double> const& m) {
    // TODO: Implement

    return ublas::matrix<double>();
}
