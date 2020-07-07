#pragma once

#include <boost/numeric/ublas/matrix.hpp>

class Distribution {
public:
    Distribution(boost::numeric::ublas::matrix<double> time_series);

	virtual boost::numeric::ublas::matrix<double> calc_num_hessian(boost::numeric::ublas::vector<double> const& x);
	virtual boost::numeric::ublas::vector<double> calc_num_gradients(boost::numeric::ublas::vector<double> const& x);
	virtual boost::numeric::ublas::vector<double> calc_gradients(boost::numeric::ublas::vector<double> const& x);
	virtual double function_value(boost::numeric::ublas::vector<double> const& x);
	virtual double calc_step_size(
        boost::numeric::ublas::vector<double> const& x, 
        boost::numeric::ublas::vector<double> const& d
    );
	~Distribution(void);

private:
    boost::numeric::ublas::vector<double> time_series;
};
