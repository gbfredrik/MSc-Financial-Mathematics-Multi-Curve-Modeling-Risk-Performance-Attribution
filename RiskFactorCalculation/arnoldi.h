#ifndef ARNOLDI
#define ARNOLDI

#include <boost/numeric/ublas/matrix.hpp>
#include <Eigen/Core>


class arnoldi {
private:

public:
	static std::tuple<boost::numeric::ublas::matrix<double>, boost::numeric::ublas::vector<double>> iram(boost::numeric::ublas::matrix<double> const& input, int k);
	

};

#endif