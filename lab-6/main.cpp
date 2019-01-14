#include <iostream>
#include <iomanip>
#include <cmath>

using namespace std;

template<typename ValTy, typename IntTy, typename FuncTy, typename ExactTy>
class runge_kutta {
    ValTy error;
    ValTy a;
    ValTy b;
    IntTy init_nsteps;
    IntTy nsteps;
    ValTy y0;
    ValTy t0;
    FuncTy dt;
    ExactTy exact;
    const ValTy r = 4;
    const ValTy R_denom = 15;
public:

    runge_kutta(ValTy eps, ValTy a, ValTy b, IntTy nsteps, ValTy y0, ValTy t0, FuncTy dt, ExactTy exact):
        a(a), b(b), init_nsteps(nsteps), nsteps(nsteps), y0(y0), t0(t0), dt(dt), exact(exact) {
        reach_error(eps);
    }


    const ValTy& get_error() const { return error; }
    const IntTy& get_nsteps() const { return nsteps; }
    void set_r_root_step() {
        nsteps = init_nsteps = 1+ (b - a) / sqrt(sqrt(error));
    }

    void double_count() {
        IntTy n = init_nsteps;
        ValTy delta;
        ValTy y_h, y_h_2 = eval(n);
        do {
            y_h = y_h_2;
            n *= 2;
            y_h_2 = eval(n);
            delta = abs(y_h_2 - y_h) / R_denom;
        } while (delta > error);
        nsteps = n;
        error = delta;
    }
    void next(const IntTy& i, const ValTy& h, const ValTy& x0, ValTy& xk, ValTy& yk, ValTy& tk) {
        ValTy k1 = h * dt(xk, yk, tk);
        ValTy k2 = h * dt(xk + h / 2, yk + k1 / 2, tk);
        ValTy k3 = h * dt(xk + h / 2, yk + k2 / 2, tk);
        ValTy k4 = h * dt(xk, yk + k3, tk);
        ValTy delta_tk = 1. / 6 * (k1 + 2 * k2 + 2 * k3 + k4);
        tk = delta_tk + tk;
        yk = h * tk + yk;
        xk = h * i + x0;
    }
    template<typename FTy>
    ValTy eval(IntTy n, FTy callable) {
        const ValTy& x0 = a;
        ValTy h = (b - a) / n;

        ValTy xk = x0, yk = y0, tk = t0;
        callable(xk, yk);
        for (IntTy i = 1; i <= n; ++i) {
            next(i, h, x0, xk, yk, tk);
            callable(xk, yk);
        }
        return yk;
    }
protected:

    ValTy eval(IntTy n) {
        const ValTy& x0 = a;
        ValTy h = (b - a) / n;

        ValTy yk = y0, tk = t0;
        for (IntTy i = 1; i <= n; ++i) {
            ValTy xk = h * i + x0;
            ValTy k1 = h * dt(xk, yk, tk);
            ValTy k2 = h * dt(xk + h / 2, yk + k1 / 2, tk);
            ValTy k3 = h * dt(xk + h / 2, yk + k2 / 2, tk);
            ValTy k4 = h * dt(xk, yk + k3, tk);
            ValTy delta_tk = 1. / 6 * (k1 + 2 * k2 + 2 * k3 + k4);
            tk = delta_tk + tk;
            yk = h * tk + yk;
        }
        return yk;
    }
    void reach_error(ValTy eps) {
        ValTy delta, approx_val, exact_val = exact(b);
        IntTy n = init_nsteps;
        do {
            approx_val = eval(n);
            delta = abs(approx_val - exact_val);
            //cout << approx_val << ' ' << exact_val << ' ' << delta << endl;
            ++n;
            //n *= 2;
        } while (delta > eps);
        error = delta;
        nsteps = n;
    }
};

template<typename ValTy, typename IntTy, typename FuncTy, typename ExactTy>
runge_kutta<ValTy, IntTy, FuncTy, ExactTy>
make_runge_kutta(ValTy eps, ValTy a, ValTy b, IntTy n, ValTy y0, ValTy t0, FuncTy dt, ExactTy exact) {
    return { eps, a, b, n, y0, t0, dt, exact };
}


int main() {
    // Variant 6
    double a = 0, b = 1;
    double y0 = 1, dy0 = 0;
    // ddy = x * dy + y + 1
    // t = dy; dy = t
    // dt = x * t + y + 1
    auto dt = [](double x, double y, double t) {
        return x * t + y + 1;
    };
    auto exact = [](double x) {
        return 2 * exp(x*x / 2) - 1;
    };
    int nsteps = 10;
    cout << "Initial number of steps: " << nsteps << endl;
    auto rk = make_runge_kutta(0.01, a, b, nsteps, y0, dy0, dt, exact);
    cout << "Reached error: " << rk.get_error() << " with n = " << rk.get_nsteps() << endl;
    cout << "Double count: " << endl;
    rk.set_r_root_step();
    cout << "Initial number of steps: " << rk.get_nsteps() << endl;
    rk.double_count();
    cout << "Reached |y_h - y_h_2|/15: " << rk.get_error() << " with n = " << rk.get_nsteps() << endl;

    cout << setw(10) << "x" << " " << setw(10) << "y" << endl;
    rk.eval(10, [](const double& x, const double& y) {
        cout << setw(10) << x << " " << setw(10) << y << endl;
    });

    return 0;
}
