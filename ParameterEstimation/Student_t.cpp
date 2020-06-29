#include "Student_t.h"

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/variance.hpp>
#include <boost/bind.hpp>

#include <numeric>
#include <cmath>

using namespace boost::numeric;

Student_t::Student_t(ublas::matrix<double> series) : Distribution(series) {
    ublas::matrix_column<ublas::matrix<double>> x(series, 0);
    time_series = x;

    ublas::vector<double> garch_vec(time_series.size() + 1);
    m_GARCH_vec = garch_vec;

    //Calc variance for timeseries
    boost::accumulators::accumulator_set<double, boost::accumulators::stats<boost::accumulators::tag::variance>> acc;
    for_each(time_series.begin(), time_series.end(), boost::bind<void>(boost::ref(acc), _1));

    //Set variance as first element
    m_GARCH_vec(0) = boost::accumulators::variance(acc);
}

void Student_t::update_GARCH_vec(ublas::vector<double> const& x) {
    ublas::vector<double> garch_vec(time_series.size() + 1);

    for (size_t i{ 0 }, n{ time_series.size() }; i < n; ++i) {
        m_GARCH_vec(i + 1) = x(0) + x(1) * pow(time_series(i), 2) + x(2) * m_GARCH_vec(i);
    }
}

double Student_t::function_value(ublas::vector<double> const& x) {
    update_GARCH_vec(x);

    double sum{ 0.0 };
    double constant{ lgamma((x(4) + 1) / 2) - lgamma(x(4) / 2) - 0.5 * log(M_PI)
        - 0.5 * log(x(4)) };

    for (size_t i{ 0 }, n{ m_GARCH_vec.size() - 1}; i < n; ++i) {
        sum = sum + constant - 0.5 * log(m_GARCH_vec(i)) - 0.5 * (x(4) + 1) * log(1 + pow(time_series(i) - x(3), 2) / (x(4) * m_GARCH_vec(i)));
    } 

    return -sum;
}

ublas::vector<double> Student_t::calcGradients(ublas::vector<double> const& x) {
    update_GARCH_vec(x);

    double dw{ 0.0 };
    double da{ 0.0 };
    double db{ 0.0 };
    double dmu{ 0.0 };
    double ddf{ 0.0 };

    ublas::vector<double> inner_w{ derivative_w(x) };
    ublas::vector<double> inner_a{ derivative_a(x) };
    ublas::vector<double> inner_b{ derivative_b(x) };

    double temp{ 0.0 };
    for (size_t i{ 0 }, n{ m_GARCH_vec.size() - 1 }; i < n; ++i) {
        temp = 0.5 * (1 / m_GARCH_vec(i) + (x(4) + 1) * ((m_GARCH_vec(i) * x(4))
            / (m_GARCH_vec(i) * x(4) + pow(time_series(i) - x(3), 2)))
            * (-pow(time_series(i) - x(3), 2)) / (pow(m_GARCH_vec(i), 2) * x(4)));

        dmu += - (x(4) + 1) * (time_series(i) - x(3)) / (m_GARCH_vec(i) * x(4) + pow(time_series(i) - x(3), 2));

        dw += temp * inner_w(i);
        da += temp * inner_a(i);
        db += temp * inner_b(i);
    }

    ublas::vector<double> gradients(x.size());
    gradients(0) = dw;
    gradients(1) = da;
    gradients(2) = db;
    gradients(3) = dmu;
    gradients(4) = calcNumGradients(x)(4);

    return gradients;
}

ublas::vector<double> Student_t::calcNumGradients(ublas::vector<double> const& x) {
    double epsilon{ 2.2 * pow(10, -16) };
    ublas::vector<double> increment{ sqrt(epsilon) * x };
    ublas::vector<double> num_gradients(x.size());

    ublas::vector<double> x_0diff{ x };
    x_0diff(0) += increment(0);
    ublas::vector<double> x_1diff{ x };
    x_1diff(1) += increment(1);
    ublas::vector<double> x_2diff{ x };
    x_2diff(2) += increment(2);
    ublas::vector<double> x_3diff{ x };
    x_3diff(3) += increment(3);
    ublas::vector<double> x_4diff{ x };
    x_4diff(4) += increment(4);

    double f_val{ function_value(x) };
    num_gradients(0) = (function_value(x_0diff) - f_val) / (x_0diff(0) - x(0));
    num_gradients(1) = (function_value(x_1diff) - f_val) / (x_1diff(1) - x(1));
    num_gradients(2) = (function_value(x_2diff) - f_val) / (x_2diff(2) - x(2));
    num_gradients(3) = (function_value(x_3diff) - f_val) / (x_3diff(3) - x(3));
    num_gradients(4) = (function_value(x_4diff) - f_val) / (x_4diff(4) - x(4));

    return num_gradients;
}

double Student_t::calcStepSize(
    ublas::vector<double> const& x, 
    ublas::vector<double> const& d
) {
    double a{ 1.0 };
    double c1{ pow(10, -4) };
    double c2{ 0.9 };

    while (x(0) + a * d(0) < 0 
        || x(1) + a * d(1) < 0 
        || x(2) + a * d(2) < 0 
        || x(1) + a * d(1) + x(2) + a * d(2) >= 1 
        || x(4) + a * d(4) <= 2
    ) {
        a *= 0.5;
        
        if (a == 0) {
            break;
        }
    }

    while (function_value(x + a * d) > function_value(x) + c1 * a * inner_prod(calcGradients(x), d)) {
        a *= 0.5;

        if (a == 0) {
            break;
        }
    }
    
    return a;
}

ublas::vector<double> Student_t::derivative_w(ublas::vector<double> const& x) {
    ublas::vector<double> inner_dw(time_series.size());
    inner_dw(0) = 0;

    for (size_t i{ 1 }, n{ inner_dw.size() }; i < n; ++i) {
        inner_dw(i) = 1 + x(2) * inner_dw(i - 1);
    }

    return inner_dw;
}

ublas::vector<double> Student_t::derivative_a(ublas::vector<double> const& x) {
    ublas::vector<double> inner_da(time_series.size());
    inner_da(0) = 0;

    for (size_t i{ 1 }, n{ inner_da.size() }; i < n; ++i) {
        inner_da(i) = pow(time_series(i - 1), 2) + x(2) * inner_da(i - 1);
    }

    return inner_da;
}

ublas::vector<double> Student_t::derivative_b(ublas::vector<double> const& x) {
    ublas::vector<double> inner_db(time_series.size());
    inner_db(0) = 0;

    for (size_t i{ 1 }, n{ inner_db.size() }; i < n; ++i) {
        inner_db(i) = m_GARCH_vec(i - 1) + x(2) * inner_db(i - 1);
    }

    return inner_db;
}
