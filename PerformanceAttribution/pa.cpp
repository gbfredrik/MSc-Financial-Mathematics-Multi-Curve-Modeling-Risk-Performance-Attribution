#include "pch.h"
#include "mex.h"

#include "pa.h"

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/optional/optional.hpp>

using namespace boost::numeric::ublas;



void pa::performanceAttribution(vector<double> const& N, vector<double> const& y, vector<matrix<double>> const& E,
	vector<size_t> k, vector<vector<double>> const& floatCashFlows, vector<vector<double>> const& fixCashFlows, 
	vector<matrix<double>> const& curveData, vector<double> const& times, int startDate, int endDate, vector<double>& NPV, 
	vector<double>& carry, vector<double>& sumRiskFactors, vector<double>& epsI, vector<double>& epsA, vector<double>& epsP) {



}