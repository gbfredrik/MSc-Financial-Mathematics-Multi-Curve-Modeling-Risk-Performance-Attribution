#include "Distribution.h"

#include <boost/numeric/ublas/matrix.hpp>

//Rosenbrock function: f(x,y) = (1-x)^2 + 100(y-x^2)^2 

class bfgs {
public:
    static boost::numeric::ublas::vector<double> minimize(
        boost::numeric::ublas::vector<double> x,
        boost::numeric::ublas::matrix<double> H_inv,
        int const max_iter,
        double const epsilon,
        Distribution* dist
    );

    //static Distribution dist;
    //static int function_type;
    //static vector<double> calcGradients(vector<double> x);
    //static double calcStepSize(vector<double> x, vector<double> d, Distribution* dist);
    //static double f(vector<double> x);
    //static double rosenbrock(vector<double> x);
};
