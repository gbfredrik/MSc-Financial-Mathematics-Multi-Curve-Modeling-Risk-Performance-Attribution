#pragma once

#include <boost/numeric/ublas/matrix.hpp>

#include <string>

struct CurveCollection {
    // General type data
    std::string file_path;
    std::string file_name;

    // Raw curve data
    boost::numeric::ublas::matrix<double> m_A;
    boost::numeric::ublas::matrix<double> m_diff;
    boost::numeric::ublas::matrix<double> m_A_trunc;
    boost::numeric::ublas::matrix<double> m_Avg;
    boost::numeric::ublas::matrix<double> m_Avg_short;

    // Risk factor data
    int k{ 6 };
    boost::numeric::ublas::matrix<double> m_E;
    boost::numeric::ublas::vector<double> v_Lambda;
    boost::numeric::ublas::matrix<double> m_delta_xi;
    boost::numeric::ublas::vector<double> v_norm_errors;
    double approximation_error{ 1.0 };

    // Parameter fitting data
    boost::numeric::ublas::vector<double> optimal_parameter_set;
    double optimal_f_val{ 0.0 };

    // Simulation data


    // Probably not: Risk measurement data


    // Probably not: Performance attribution data

};
