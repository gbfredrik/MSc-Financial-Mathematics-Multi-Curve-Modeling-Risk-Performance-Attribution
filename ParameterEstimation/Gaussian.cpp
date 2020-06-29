#include "Gaussian.h"

#include "../MathLibrary/matrixOperations.h"

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/variance.hpp>
#include <boost/bind.hpp>

#include <cmath>

using namespace boost::numeric;

Gaussian::Gaussian(ublas::matrix<double> series) : Distribution(series) {
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

//Update garch vector with new parameters
void Gaussian::update_GARCH_vec(ublas::vector<double> const& x) {  // datum växer med index
    for (size_t i{ 0 }, n{ time_series.size() }; i < n; ++i) {
        m_GARCH_vec(i + 1) = x(0) + x(1) * pow(time_series(i), 2) + x(2) * m_GARCH_vec(i);
    }
}

//Calculate function value given parameters x
double Gaussian::function_value(ublas::vector<double> const& x) {
    update_GARCH_vec(x);

    double sum{ 0.0 };
    for (size_t i{ 0 }, n{ m_GARCH_vec.size() - 1 }; i < n; ++i) {
        sum += log(2 * M_PI * m_GARCH_vec(i)) + pow(time_series(i) - x(3), 2) / (m_GARCH_vec(i));
    }

    return 0.5 * sum;
}

//Calculate gradients, x = [omega alpha beta mu]
ublas::vector<double> Gaussian::calcGradients(ublas::vector<double> const&  x) {
    update_GARCH_vec(x);
    
    double dw{ 0.0 };
    double da{ 0.0 };
    double db{ 0.0 };
    double dmu{ 0.0 };

    ublas::vector<double> inner_w{ derivative_w(x) };
    ublas::vector<double> inner_a{ derivative_a(x) };
    ublas::vector<double> inner_b{ derivative_b(x) };

    double temp{ 0.0 };
    for (size_t i{ 0 }, n{ m_GARCH_vec.size() - 1 }; i < n; ++i) {
        temp = 0.5 * (m_GARCH_vec(i) - pow(time_series(i) - x(3), 2)) / (pow(m_GARCH_vec(i), 2));

        dw += temp * inner_w(i);
        da += temp * inner_a(i);
        db += temp * inner_b(i);
        dmu -= (time_series(i) - x(3)) / m_GARCH_vec(i);
    }

    ublas::vector<double> gradients(4);
    gradients(0) = dw;
    gradients(1) = da;
    gradients(2) = db;
    gradients(3) = dmu;

    return gradients;
}

ublas::vector<double> Gaussian::calcNumGradients(ublas::vector<double> const& x) {
    double epsilon{ 2.2 * pow(10, -16) };
    ublas::vector<double> increment{ sqrt(epsilon) * x };
    ublas::vector<double> num_gradients(4);

    ublas::vector<double> x_0diff{ x };
    x_0diff(0) += increment(0);
    ublas::vector<double> x_1diff{ x };
    x_1diff(1) += increment(1);
    ublas::vector<double> x_2diff{ x };
    x_2diff(2) += increment(2);
    ublas::vector<double> x_3diff{ x };
    x_3diff(3) += increment(3);

    double f_val{ function_value(x) };
    num_gradients(0) = (function_value(x_0diff) - f_val) / increment(0);
    num_gradients(1) = (function_value(x_1diff) - f_val) / increment(1);
    num_gradients(2) = (function_value(x_2diff) - f_val) / increment(2);
    num_gradients(3) = (function_value(x_3diff) - f_val) / increment(3);

    return num_gradients;
}

// Hessian: https://v8doc.sas.com/sashtml/ormp/chap5/sect28.htm
//
ublas::matrix<double> Gaussian::calcNumHessian(ublas::vector<double> const& x) {
    //ublas::vector<double> num_gradients(4);
    ublas::vector<double> h(4);
    for (size_t i{ 0 }, n{ h.size() }; i < n; ++i) {
        h(i) = 0;
    }

    double epsilon{ 4.8 * pow(10, -6) };
    
    ublas::vector<double> inc(4);
    inc(0) = 0.0000000000001;
    inc(1) = 0.0000000000001;
    inc(2) = 0.0000000000001;
    inc(3) = 0.0000000000001;
    
    ublas::vector<double> hw{ h };
    hw(0) += inc(0);
    ublas::vector<double> ha{ h };
    ha(1) += inc(1);
    ublas::vector<double> hb{ h };
    hb(2) += inc(2);
    ublas::vector<double> hmu{ h };
    hmu(3) += inc(3);

    ublas::vector<double> grad_x_add_hw{ calcGradients(x + hw) };
    ublas::vector<double> grad_x_add_ha{ calcGradients(x + ha) };
    ublas::vector<double> grad_x_add_hb{ calcGradients(x + hb) };
    ublas::vector<double> grad_x_add_hmu{ calcGradients(x + hmu) };

    ublas::vector<double> grad_x_sub_hw{ calcGradients(x - hw) };
    ublas::vector<double> grad_x_sub_ha{ calcGradients(x - ha) };
    ublas::vector<double> grad_x_sub_hb{ calcGradients(x - hb) };
    ublas::vector<double> grad_x_sub_hmu{ calcGradients(x - hmu) };

    double dw2{ (grad_x_add_hw(0) - grad_x_sub_hw(0)) / (4.0 * inc(0)) 
        + (grad_x_add_hw(0) - grad_x_sub_hw(0)) / (4.0 * inc(0)) };
    double dwa{ (grad_x_add_ha(0) - grad_x_sub_ha(0)) / (4.0 * inc(1)) 
        + (grad_x_add_hw(1) - grad_x_sub_hw(1)) / (4.0 * inc(0)) };
    double dwb{ (grad_x_add_hb(0) - grad_x_sub_hb(0)) / (4.0 * inc(2)) 
        + (grad_x_add_hw(2) - grad_x_sub_hw(2)) / (4.0 * inc(0)) };
    double dwmu{ (grad_x_add_hmu(0) - grad_x_sub_hmu(0)) / (4.0 * inc(3)) 
        + (grad_x_add_hw(3) - grad_x_sub_hw(3)) / (4.0 * inc(0)) };
    double dab{ (grad_x_add_hb(1) - grad_x_sub_hb(1)) / (4.0 * inc(2)) 
        + (grad_x_add_ha(2) - grad_x_sub_ha(2)) / (4.0 * inc(1)) };
    double damu{ (grad_x_add_hmu(1) - grad_x_sub_hmu(1)) / (4.0 * inc(3)) 
        + (grad_x_add_ha(3) - grad_x_sub_ha(3)) / (4.0 * inc(1)) };
    double dbmu{ (grad_x_add_hmu(2) - grad_x_sub_hmu(2)) / (4.0 * inc(3)) 
        + (grad_x_add_hb(3) - grad_x_sub_hb(3)) / (4.0 * inc(2)) };
    double da2{ (grad_x_add_ha(1) - grad_x_sub_ha(1)) / (4.0 * inc(1)) 
        + (grad_x_add_ha(1) - grad_x_sub_ha(1)) / (4.0 * inc(1)) };
    double db2{ (grad_x_add_hb(2) - grad_x_sub_hb(2)) / (4.0 * inc(2)) 
        + (grad_x_add_hb(2) - grad_x_sub_hb(2)) / (4.0 * inc(2)) };
    double dmu2{ (grad_x_add_hmu(3) - grad_x_sub_hmu(3)) / (4.0 * inc(3)) 
        + (grad_x_add_hmu(3) - grad_x_sub_hmu(3)) / (4.0 * inc(3)) };

    ublas::matrix<double> Hessian(4, 4);
    Hessian(0, 0) = dw2;
    Hessian(0, 1) = dwa;
    Hessian(0, 2) = dwb;
    Hessian(0, 3) = dwmu;
    Hessian(1, 0) = dwa;
    Hessian(1, 1) = da2;
    Hessian(1, 2) = dab;
    Hessian(1, 3) = damu;
    Hessian(2, 0) = dwb;
    Hessian(2, 1) = dab;
    Hessian(2, 2) = db2;
    Hessian(2, 3) = dbmu;
    Hessian(3, 0) = dwmu;
    Hessian(3, 1) = damu;
    Hessian(3, 2) = dbmu;
    Hessian(3, 3) = dmu2;

    return matrixOperations::matrix_inv(Hessian);
}

double Gaussian::calcStepSize(
    ublas::vector<double> const& x, 
    ublas::vector<double> const& d
) {
    double a{ 1.0 };
    double c1{ pow(10, -3) };
    double c2{ 0.9 };

    while (x(0) + a * d(0) < 0 
        || x(1) + a * d(1) < 0 
        || x(2) + a * d(2) < 0 
        || x(1) + a * d(1) + x(2) + a * d(2) >= 1
    ) {
        a *= 0.5;

        if (a == 0) {
            break;
        }
    }

    while (function_value(x + a * d) > function_value(x))
    {
        a *= 0.5;
    }

    return a;
}

//Calculate vector with derivatives of garch with respect to omega. dv/dw.
ublas::vector<double> Gaussian::derivative_w(ublas::vector<double> const& x) {
    ublas::vector<double> inner_dw(time_series.size());
    inner_dw(0) = 0;

    for (size_t i{ 1 }, n{ inner_dw.size() }; i < n; ++i) {
        inner_dw(i) = 1 + x(2) * inner_dw(i - 1);
    }

    return inner_dw;
}

//Calculate vector with derivatives of garch with respect to alpha. dv/da.
ublas::vector<double> Gaussian::derivative_a(ublas::vector<double> const& x) {
    ublas::vector<double> inner_da(time_series.size());
    inner_da(0) = 0;

    for (size_t i{ 1 }, n{ inner_da.size() }; i < n; ++i) {
        inner_da(i) = pow(time_series(i - 1), 2) + x(2) * inner_da(i - 1);
    }

    return inner_da;
}

//Calculate vector with derivatives of garch with respect to beta. dv/db
ublas::vector<double> Gaussian::derivative_b(ublas::vector<double> const& x) {
    ublas::vector<double> inner_db(time_series.size());
    inner_db(0) = 0;

    for (size_t i{ 1 }, n{ inner_db.size() }; i < n; ++i) {
        inner_db(i) = m_GARCH_vec(i - 1) + x(2) * inner_db(i - 1);
    }

    return inner_db;
}
