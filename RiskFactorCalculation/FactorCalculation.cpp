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
    ublas::vector<double>& v_Lambda,
    double& approximation_error,
    ublas::vector<double>& v_norm_errors
) {
    // Utilizes the SymEigsSolver (IRAM) from the Spectra library
    using namespace Spectra;

    // Transform ublas matrix into Eigen::matrix
    Eigen::MatrixXd Data{ matrixOperations::ublasToMatrixXd(input) };
    Eigen::MatrixXd U{ };
    Eigen::VectorXd D{ };
    Eigen::MatrixXd V{ };

    if (Data.rows() < Data.cols()) {
        Eigen::HouseholderQR<Eigen::MatrixXd> qr{ Data.transpose() };
        // Can use FullPivHouseholderQR for greater precision
        //Eigen::FullPivHouseholderQR<Eigen::MatrixXd> qr(Data.transpose());

        // Define the positive definite RTR matrix
        Eigen::MatrixXd thinQ{ qr.householderQ() * Eigen::MatrixXd::Identity(Data.cols(), Data.rows()) }; // qr.householderQ()
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
            U = thinQ * eigs.eigenvectors();
            V = thinQ * eigs.eigenvectors();
            D = eigs.eigenvalues().array().sqrt();
        }
    } else {
        // Define the positive definite C matrix (covariance)
        Eigen::MatrixXd C{ Data.transpose() * Data };

        // Construct matrix operation objects
        DenseSymMatProd<double> op(C);

        // Construct eigen solver object, requesting the largest k eigenvalues in magnitude
        SymEigsSolver<double, LARGEST_ALGE, DenseSymMatProd<double>> eigs(&op, k, 2 * k + 1);

        eigs.init();
        int nconv{ eigs.compute() };

        // Retrieve results
        if (eigs.info() == SUCCESSFUL) {
            U = eigs.eigenvectors();
            V = eigs.eigenvectors();
            D = eigs.eigenvalues().array().sqrt();
        }
    }

    // Assign results to referenced variables
    m_E = matrixOperations::matrixXdToUblas(V);
    v_Lambda = matrixOperations::vectorXdToUblas(D.array().square());
    if (approximation_error == 1.0) {
        Eigen::MatrixXd C{ Data.transpose() * Data };
        Eigen::MatrixXd rsvdApprox{ U * D.asDiagonal() * D.asDiagonal() * V.adjoint() };
        approximation_error = Rsvd::relativeFrobeniusNormError(C, rsvdApprox);
        
        v_norm_errors = eig_all_norm_errors(matrixOperations::matrixXdToUblas(C), m_E, v_Lambda);
    }
    
    return v_Lambda.size() > 0;
}

bool FactorCalculation::eigen_bdcsvd(
    ublas::matrix<double> const& input, 
    int const k, 
    ublas::matrix<double>& m_E, 
    ublas::vector<double>& v_Lambda,
    double& approximation_error,
    ublas::vector<double>& v_norm_errors
) {
    // Utilizes BDCSVD from Eigen library
    int svd_opt{ Eigen::ComputeThinU | Eigen::ComputeThinV };

    Eigen::MatrixXd H{ matrixOperations::ublasToMatrixXd(input) };
    Eigen::MatrixXd U{ };
    Eigen::VectorXd D{ };
    Eigen::MatrixXd V{ };

    if (H.rows() < H.cols()) {
        Eigen::HouseholderQR<Eigen::MatrixXd> qr(H.transpose()); // FullPivHouseholderQR<Matrix<double, Dynamic, Size>> fpqr(A.rows(), A.cols());
        Eigen::MatrixXd thinQ{ qr.householderQ() * Eigen::MatrixXd::Identity(H.cols(), H.rows()) };
        Eigen::MatrixXd R_temp{ qr.matrixQR().triangularView<Eigen::Upper>() };
        Eigen::MatrixXd thinR{ R_temp.topRows(R_temp.cols()) };
        Eigen::MatrixXd RRT{ thinR * thinR.transpose() };

        Eigen::BDCSVD<Eigen::MatrixXd> bdcsvd(RRT, svd_opt);

        // Retrieve results
        if (k == 0) {
            U = thinQ * bdcsvd.matrixU();
            V = thinQ * bdcsvd.matrixV();
            D = bdcsvd.singularValues().array().sqrt();
        } else {
            U = thinQ * bdcsvd.matrixU().block(0, 0, bdcsvd.matrixU().rows(), k);
            V = thinQ * bdcsvd.matrixV().block(0, 0, bdcsvd.matrixV().rows(), k);
            D = bdcsvd.singularValues().head(k).array().sqrt();
            // See: https://stackoverflow.com/questions/34373757/piece-wise-square-of-vector-piece-wise-product-of-two-vectors-in-c-eigen
        }
    } else {
        Eigen::BDCSVD<Eigen::MatrixXd> bdcsvd(H.transpose() * H, svd_opt);
        
        if (k == 0) {
            U = bdcsvd.matrixU();
            V = bdcsvd.matrixV();
            D = bdcsvd.singularValues();
        } else {
            U = bdcsvd.matrixU().block(0, 0, bdcsvd.matrixU().rows(), k);
            V = bdcsvd.matrixV().block(0, 0, bdcsvd.matrixV().rows(), k);
            D = bdcsvd.singularValues().head(k);
        }
    }
    
    // Assign results to referenced variables
    m_E = matrixOperations::matrixXdToUblas(V);
    v_Lambda = matrixOperations::vectorXdToUblas(D.array().square());
    if (approximation_error == 1.0) {
        Eigen::MatrixXd C{ H.transpose() * H };
        Eigen::MatrixXd rsvdApprox{ U * D.asDiagonal() * D.asDiagonal() * V.adjoint() };
        approximation_error = Rsvd::relativeFrobeniusNormError(C, rsvdApprox);
        
        v_norm_errors = eig_all_norm_errors(matrixOperations::matrixXdToUblas(C), m_E, v_Lambda);
    }


    return v_Lambda.size() > 0;
}

bool FactorCalculation::eigen_rsvd(
    ublas::matrix<double> const& input,
    int const k,
    ublas::matrix<double>& m_E,
    ublas::vector<double>& v_Lambda,
    double& approximation_error,
    ublas::vector<double>& v_norm_errors
) {
    // Utilizes Randomized SVD from Eigen extension
    int svd_opt{ Eigen::ComputeThinU | Eigen::ComputeThinV };

    Eigen::MatrixXd H{ matrixOperations::ublasToMatrixXd(input) };
    const Eigen::Index reducedRank{ k };
    Eigen::MatrixXd U{ };
    Eigen::VectorXd D{ };
    Eigen::MatrixXd V{ };

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

        // Retrieve results
        U = thinQ * rsvd.matrixU();
        V = thinQ * rsvd.matrixV();
        D = rsvd.singularValues().array().sqrt();
    } else {
        Rsvd::RandomizedSvd<Eigen::MatrixXd, std::mt19937_64, Rsvd::SubspaceIterationConditioner::Lu>
            rsvd(randomEngine);
        rsvd.compute(H, reducedRank);

        // Retrieve results
        U = rsvd.matrixU();
        V = rsvd.matrixV();
        D = rsvd.singularValues();
    }
    
    // Assign results to referenced variables
    m_E = matrixOperations::matrixXdToUblas(V);
    v_Lambda = matrixOperations::vectorXdToUblas(D.array().square());
    if (approximation_error == 1.0) {
        Eigen::MatrixXd C{ H.transpose() * H };
        Eigen::MatrixXd rsvdApprox{ svd_approximation(U, D, V) };
        approximation_error = Rsvd::relativeFrobeniusNormError(C, rsvdApprox);
        
        v_norm_errors = eig_all_norm_errors(matrixOperations::matrixXdToUblas(C), m_E, v_Lambda);
    }

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

Eigen::MatrixXd FactorCalculation::svd_approximation(
    Eigen::MatrixXd const& m_U,
    Eigen::MatrixXd const& v_D,
    Eigen::MatrixXd const& m_V
) {
    return m_U * v_D.asDiagonal() * v_D.asDiagonal() * m_V.transpose();
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

double FactorCalculation::relativeFrobeniusNormError(
    Eigen::MatrixXd const& m_original, 
    Eigen::MatrixXd const& m_approximation
) {
    const auto differenceNorm{ (m_approximation - m_original).norm() };
    const auto referenceNorm{ m_original.norm() };

    assert(referenceNorm > 0);

    return differenceNorm / referenceNorm;
}

ublas::matrix<double> FactorCalculation::clean_data(ublas::matrix<double> const& m) {
    // TODO: Implement

    return ublas::matrix<double>();
}
