#include "bfgs.h"

#include "Distribution.h"

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>

using namespace boost::numeric;

ublas::vector<double> bfgs::minimize(
    ublas::vector<double> x,
    ublas::matrix<double> H_inv,
    int const max_iter,
    double const epsilon,
    Distribution* dist
) {
    size_t n{ x.size() };
    ublas::matrix<double> help_prod(n, n);
    ublas::vector<double> d(n);
    ublas::vector<double> x_new(n);
    ublas::vector<double> y(n);
    ublas::vector<double> s(n);
    ublas::vector<double> gradient_vec(n);
    double alpha{};
    double scale_H{};
    double l{};
    int k{ 0 };
    //Init hessian inverse matrix
    ublas::identity_matrix<double> I(n);

    gradient_vec = dist->calcGradients(x);

    while (ublas::norm_2(gradient_vec) > epsilon && k < max_iter) {
        //gradient_vec = dist->calcGradients(x);
        d = -prod(H_inv, gradient_vec);
        alpha = dist->calcStepSize(x, d);
        s = alpha * d;
        x_new = x + s;
        //x_new = x + alpha * d;
        //std::cout << std::setprecision(16) << "Nya s   = " << s << std::endl;
        //std::cout << std::setprecision(16) << "Gamla s = " << x_new - x << std::endl << std::endl;

        ublas::vector<double> gradient_vec_new{ dist->calcGradients(x_new) };
        y = gradient_vec_new - gradient_vec;
        //s = x_new - x; // Gammal lösning men verkar sämre
        scale_H = inner_prod(y, s);
        l = 1.0 / scale_H;

        if (scale_H == 0) {
            std::cout << "Breaking since H has invalid values.\n\n";
            break;
        }

        help_prod = prod((I - l * outer_prod(s, y)), H_inv);
        H_inv = prod(help_prod, (I - l * outer_prod(y, s))) + l * outer_prod(s, s);
        x = x_new;
        gradient_vec = gradient_vec_new;
        ++k;
    }

    ublas::vector<double> results(n + 1);

    for (size_t i{ 0 }; i < results.size() - 1; ++i) {
        results(i) = x(i);
    }
    results(n) = dist->function_value(x);

    return results;
}
