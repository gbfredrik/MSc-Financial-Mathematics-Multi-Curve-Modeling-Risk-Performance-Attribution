#pragma once

#include <boost/numeric/ublas/matrix.hpp>

class Distribution {
public:
    Distribution(boost::numeric::ublas::matrix<double> time_series);

	virtual boost::numeric::ublas::matrix<double> calcNumHessian(boost::numeric::ublas::vector<double> const& x);
	virtual boost::numeric::ublas::vector<double> calcNumGradients(boost::numeric::ublas::vector<double> const& x);
	virtual boost::numeric::ublas::vector<double> calcGradients(boost::numeric::ublas::vector<double> const& x);
	virtual double function_value(boost::numeric::ublas::vector<double> const& x);
	virtual double calcStepSize(
        boost::numeric::ublas::vector<double> const& x, 
        boost::numeric::ublas::vector<double> const& d
    );
	~Distribution(void);

private:
    boost::numeric::ublas::vector<double> time_series;
};
