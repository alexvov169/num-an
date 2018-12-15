#include <iostream>
#include <fstream>
#include <cmath>

#include "simpson.h"
#include "sole.h"
#include "polynomials.h"

using namespace std;

template<typename ValTy, typename IntTy, typename PolynomsTy, typename ApprFunc>
class normal_form_impl : public sole<ValTy, IntTy> {
public:
    normal_form_impl(ValTy a, ValTy b,
                     IntTy form_degree,
                     ValTy int_relative_error,
                     PolynomsTy polynomial_basis,
                     ApprFunc f):
        sole<ValTy, IntTy>(form_degree),
        a(a), b(b),
        int_relative_error(int_relative_error),
        polynomial_basis(polynomial_basis), form_degree(form_degree),
        f(f)
    {
        build(form_degree);
        this->complete_pivoting();
    }

    normal_form_impl& operator++() {
        ++form_degree;
        build(form_degree);
        this->complete_pivoting();
        return *this;
    }

    ValTy get_value(const ValTy& x) {
        ValTy result = 0;
        IntTy i = 0;
        for (const auto& root : this->get_roots()) {
            result += polynomial_basis(i, x) * root;
            ++i;
        }
        return result;
    }

    const IntTy& get_form_degree() const { return form_degree; }

protected:
    void build(IntTy n) {
        this->resize(n);
        this->roots.resize(n);
        for (IntTy i = 0; i < n; ++i) {
            (*this)[i].resize(n+1);
        }
        int j, i;

        for (j = 0; j < n; j++) {
            auto simpson1 = simpson_integrate(a, b, [&](double x) {
                double tmp = polynomial_basis(j, x);
                return tmp * tmp;
            });
            (*this)[j][j] = simpson1.relative_runge(int_relative_error);
            for (i = j + 1; i < n; i++) {
                auto simpson = simpson_integrate(a, b, [&](double x) {
                    return polynomial_basis(j, x) * polynomial_basis(i, x);
                });
                (*this)[j][i] = (*this)[i][j] = simpson.relative_runge(int_relative_error);
            }
            auto simpson = simpson_integrate(a, b, [&](double x) {
                return polynomial_basis(j, x) * f(x);
            });
            (*this)[j][n] = simpson.relative_runge(int_relative_error);
        }
    }

private:
    ValTy a;
    ValTy b;
    ValTy int_relative_error;
    PolynomsTy polynomial_basis;
    IntTy form_degree;
    ApprFunc f;
};

template<typename ValTy, typename IntTy, typename ApprFunc, typename PolynomsTy>
class lsd_approximator_impl {
public:
    lsd_approximator_impl(ValTy a, ValTy b, ValTy error,
                          IntTy npoints, IntTy pol_degree,
                          ValTy int_relative_error,
                          PolynomsTy polynomial,
                          ApprFunc f):
        error(error),
        a(a), b(b),
        npoints(npoints),
        normal_form(a, b, pol_degree, int_relative_error, polynomial, f),
        f(f)
    {
    }

    ValTy operator()() {
        ValTy x;
        ValTy h = (b - a) / npoints;
        ValTy approximation, exact;
        ValTy deviation_sum, delta;
        ValTy deviation;
        do {
            deviation_sum = 0;
            for (IntTy i = 0; i <= npoints; ++i) {
                x = a + h * i;
                exact = f(x);
                approximation = normal_form.get_value(x);
                delta = abs(approximation - exact);
                deviation_sum = delta * delta + deviation_sum;
                //printf("%f %f %f\n", x, exact, approximation);
                //cout << "N = " << N << " exact " << exact << " apprx " << result_n  << endl;
            }
            deviation = sqrt(deviation_sum / (normal_form.get_form_degree() + 1));
            cout << "N = " << normal_form.get_form_degree() << "; deviation = " << deviation << endl;
            ++normal_form;
        } while (deviation > error);

        ofstream approximation_out;
        approximation_out.open("points.csv");
        for (IntTy i = 0; i <= npoints; ++i) {
            x = a + h * i;
            approximation = normal_form.get_value(x);
            approximation_out << x << ';' << approximation << endl;
        }
        return approximation;
    }
private:
    ValTy error;
    ValTy a;
    ValTy b;
    IntTy npoints;
    normal_form_impl<ValTy, IntTy, PolynomsTy, ApprFunc> normal_form;

    ApprFunc f;
};

template<typename ValTy, typename IntTy, typename ApprFunc, typename PolynomsTy>
lsd_approximator_impl<ValTy, IntTy, ApprFunc, PolynomsTy>
make_lsd_approximator(ValTy a, ValTy b, ValTy error,
         IntTy npoints, IntTy pol_degree,
         ValTy int_relative_error,
         PolynomsTy polynomial,
         ApprFunc f) {
    return { a, b, error, npoints, pol_degree, int_relative_error, polynomial, f };
}

int main() {
    auto f = [&](double x) { return 1 / x - 0.1 * x * x * sin(2 * x); };

	double a = 2, b = 11;
    int degree = 10;
    int npoints = 90;
    double epsilon = 0.01;
    //auto approximator = make_lsd_approximator(a, b, epsilon, npoints, degree, 0.001, chebyshov_polynomial<double, int>, f);
    auto approximator = make_lsd_approximator(a, b, epsilon, npoints, degree, 0.001, algebraic_polynomial<double, int>, f);
    approximator();

	return 0;
}
