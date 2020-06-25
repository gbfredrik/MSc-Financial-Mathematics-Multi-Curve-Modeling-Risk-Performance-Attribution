#pragma once
#include "Distribution.h"

#include <boost/numeric/ublas/matrix.hpp>

class T_Copula : public Distribution {
public:
    T_Copula(boost::numeric::ublas::matrix<double> series);
    double function_value(boost::numeric::ublas::vector<double> const& x);
    boost::numeric::ublas::vector<double> calcGradients(boost::numeric::ublas::vector<double> const& x);
    boost::numeric::ublas::vector<double> calcNumGradients(boost::numeric::ublas::vector<double> const& x);
    void getSeries();
    double calcStepSize(
        boost::numeric::ublas::vector<double> const& x,
        boost::numeric::ublas::vector<double> const& d
    );

    //private:
    boost::numeric::ublas::matrix<double> time_series;
    boost::numeric::ublas::matrix<double> buildP(boost::numeric::ublas::vector<double> const& x);
    boost::numeric::ublas::vector<double> matrixToVector(boost::numeric::ublas::matrix<double> const& matrix);
    boost::numeric::ublas::matrix<double> vectorToMatrix(boost::numeric::ublas::vector<double> const& vec);
    boost::numeric::ublas::vector<double> getElements(boost::numeric::ublas::matrix<double> const& matrix);
    boost::numeric::ublas::vector<double> kronOfVectors(boost::numeric::ublas::vector<double> const& v1, boost::numeric::ublas::vector<double> const& v2);
    double dGamma(double t);
};