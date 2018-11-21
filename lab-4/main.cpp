#include <iostream>
#include <iomanip>
#include "trapezium.h"

using namespace std;

static long double f(long double x) { return 10 * log(x); }
static long double F(long double x) {
    return 10 * x * (log(x)-1); // the primitive function is miswritten in variant 6
}
static long double ddf(long double x) {
    return -10 / x / x;
}

template<typename ValTy, typename Integrator>
ValTy print_analytical_table(ostream& stream, Integrator integrate, ValTy epsilon, ValTy max_ddf) {
    stream << "analytic step" << endl;
    stream << setw(10) << "epsilon"
         << setw(20) << "step"
         << setw(20) << "exact"
         << setw(20) << "delta"
         << endl;
    auto exact_value = integrate.exact(),
         trapezium = integrate.analytical(epsilon, max_ddf),
         delta = abs(exact_value - trapezium);
    stream << setw(10) << setprecision(0) << epsilon
         << setw(20) << setprecision(10) << integrate.get_step()
         << setw(20) << exact_value
         << setw(20) << delta
         << endl;
    return delta;
}

template<typename ValTy, typename Integrator>
void print_runge_table(ostream& stream, Integrator integrate, ValTy epsilon) {
    stream << "runge step" << endl;
    stream << setw(10) << "epsilon"
         << setw(20) << "step"
         << setw(20) << "delta"
         << endl;
    auto exact_value = integrate.exact();
    auto trapezium = integrate.runge(epsilon, int(1/sqrt(epsilon)));
    stream << setw(10) << setprecision(10) << epsilon
         << setw(20) << setprecision(10) << integrate.get_step()
         << setw(20) << abs(exact_value - trapezium)
         << endl;
}

int main() {
    long double a = 1, b = 15;
    auto integrate = trapezium_integrate(f, F, a, b);
    auto epsilon_start = static_cast<long double>(1e-2),
            epsilon_mult = static_cast<long double>(1e-3), epsilon = epsilon_start;
    int nrows = 4;

    for (int i = nrows; i--;) {
        print_runge_table(cout, integrate, print_analytical_table(cout, integrate, epsilon, ddf(a)));
        epsilon *= epsilon_mult;
    }

    return 0;
}
