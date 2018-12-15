#include <iostream>
#include <vector>
#include <functional>
#include <cmath>
#include <memory>

#include "simpson.h"

using namespace std;


template<typename ValTy>
class polynomial {
    ValTy x;
public:
    polynomial(ValTy x) : x(x) {}
    const ValTy& get_x() const { return x; }
};

template<typename ValTy>
class algebraic_polynomials : public polynomial<ValTy> {
public:
    algebraic_polynomials(ValTy x) : polynomial<ValTy>(x) {}

    class iterator {
        const algebraic_polynomials *that;
        ValTy current;
    public:
        iterator(const algebraic_polynomials *that, ValTy initial):
            that(that), current(initial) {}
        const ValTy& operator*() const { return current; }
        iterator& operator++() {
            current = that->get_x() * current;
            return *this;
        }
    };

    iterator begin() const {
        return { this, 1 };
    }
};


template<typename ValTy>
algebraic_polynomials<ValTy>
make_algebraic_polynomials(ValTy x) {
    return {x};
}
/*
double chebyshov_polynomial(int n, double x) {
    auto t = make_chebyshov_polynomials(x);
    auto i = t.begin();
    while (n--) {
        ++i;
    }
    return *i;
}*/

double algebraic_polynomial(int n, double x) {
    auto t = make_algebraic_polynomials(x);
    auto i = t.begin();
    while (n--) {
        ++i;
    }
    return *i;
}

double chebyshev(int n, double x) {
	double Tn1, Tn = x, Tn_1 = 1;

	if (n == 0) return Tn_1;

	int i = 1;
	while (i < n) {
		Tn1 = 2 * x * Tn - Tn_1;
		Tn_1 = Tn;
		Tn = Tn1;
		++i;
	}

	return Tn;
}

vector<double> scheme_selection_main_element(vector<vector<double> >& matrix, int N){
    int i, j, line, row, k, line_new = 0;
    double R, max;
    vector<double> M(N), res(N);
    vector<vector<double> > new_matrix = matrix;

    for (k = 0; k < N - 1; k++) {
        max = 0;
        for (i = 0; i < N; i++) {
            for (j = 0; j < N; j++) {
                if (fabs(matrix[i][j]) > max) {
                    max = fabs(matrix[i][j]);
                    line = i;
                    row = j;
                }
            }
        }

        for (j = 0; j < N + 1; j++)
            new_matrix[line_new][j] = matrix[line][j];
        line_new++;

        for (i = 0; i < N; i++)
            M[i] = -(matrix[i][row] / max);

        for (i = 0; i < N; i++)
            for (j = 0; j < N + 1; j++)
                if (i != line)
                    matrix[i][j] += matrix[line][j] * M[i];

        for (j = 0; j < N + 1; j++)
            matrix[line][j] = 0;

        for (i = 0; i < N; i++)
            matrix[i][row] = 0;
    }

    for (i = 0; i < N; i++)
        if (matrix[i][N] != 0)
            line = i;

    for (j = 0; j < N + 1; j++)
        new_matrix[line_new][j] = matrix[line][j];

    for (j = 0; j < N; j++)
        res[j] = 0;

    for (i = N - 1; i >= 0; i--){
        R = new_matrix[i][N];
        for (j = 0; j < N; j++)
            if ((new_matrix[i][j] != 0) && (res[j] != 0))
                R -= new_matrix[i][j] * res[j];
        for (j = 0; j < N; j++)
            if ((new_matrix[i][j] != 0) && (res[j] == 0))
                res[j] = R / new_matrix[i][j];
    }

    //cout << "res = " << res << endl;
    return res;
}

template<typename IndexType>
class zero_leader_exception_impl : public std::exception {
public:
    zero_leader_exception_impl(IndexType i) : i(i) {}
    virtual const char *what() const noexcept override { return "zero leader"; }
    IndexType i;
};
template<typename IndexType>
zero_leader_exception_impl<IndexType>
zero_leader_exception(IndexType i) {
    return zero_leader_exception_impl<IndexType>(i);
}

template<typename ValTy, typename IntTy>
class sole : protected vector<vector<ValTy> > {
public:
    sole(IntTy n) : vector<vector<ValTy> >(n) {}

    const vector<ValTy>& get_roots() const { return roots; }

    void complete_pivoting1() {
        const IntTy N = this->size();
        int i, j, line, row, k, line_new = 0;
        double R, max;
        vector<double> M(N);
        vector<vector<double> > new_atrix = *this;

        for (k = 0; k < N - 1; k++) {
            max = 0;
            for (i = 0; i < N; i++) {
                for (j = 0; j < N; j++) {
                    if (fabs((*this)[i][j]) > max) {
                        max = fabs((*this)[i][j]);
                        line = i;
                        row = j;
                    }
                }
            }

            for (j = 0; j < N + 1; j++)
                new_atrix[line_new][j] = (*this)[line][j];
            line_new++;

            for (i = 0; i < N; i++)
                M[i] = -((*this)[i][row] / max);

            for (i = 0; i < N; i++)
                for (j = 0; j < N + 1; j++)
                    if (i != line)
                        (*this)[i][j] += (*this)[line][j] * M[i];

            for (j = 0; j < N + 1; j++)
                (*this)[line][j] = 0;

            for (i = 0; i < N; i++)
                (*this)[i][row] = 0;
        }

        for (i = 0; i < N; i++)
            if ((*this)[i][N] != 0)
                line = i;

        for (j = 0; j < N + 1; j++)
            new_atrix[line_new][j] = (*this)[line][j];

        for (j = 0; j < N; j++)
            roots[j] = 0;

        for (i = N - 1; i >= 0; i--){
            R = new_atrix[i][N];
            for (j = 0; j < N; j++)
                if ((new_atrix[i][j] != 0) && (roots[j] != 0))
                    R -= new_atrix[i][j] * roots[j];
            for (j = 0; j < N; j++)
                if ((new_atrix[i][j] != 0) && (roots[j] == 0))
                    roots[j] = R / new_atrix[i][j];
        }
    }

    void complete_pivoting() {
        const IntTy n = this->size();
        vector<IntTy> indices(n);
        for (IntTy m = 0; m < n; ++m) {
            ValTy maximum = (*this)[m][m];
            IntTy imax = m;
            IntTy jmax = m;
            for (IntTy i = m+1; i < n; ++i) {
                for (IntTy j = m+1; j < n; ++j) {
                    if (maximum < (*this)[i][j]) {
                        maximum = (*this)[i][j];
                        imax = i;
                        jmax = j;
                    }
                }
            }
            indices[m] = jmax;

            if (m != imax) {
                swap((*this)[m], (*this)[imax]);
            }
            if (m != jmax) {
                for (IntTy i = 0; i < n; ++i) {
                    swap((*this)[i][m], (*this)[i][jmax]);
                }
            }
        }
        gaussian_elimination();
        for (IntTy i = 0; i < n; ++i) {
            swap(roots[i], roots[indices[i]]);
        }
    }
    void gaussian_elimination() {
        roots.resize(this->size());
        auto& m = *this;
                const int n = this->size();

                for (int i = 0; i < n; ++i) {
                    if (m[i][i] != 0) {
                        const double leader = m[i][i];
                        for (int j = 0; j < n + 1; ++j) {
                            m[i][j] /= leader;
                        }
                        for (int k = i + 1; k < n; ++k) {
                            const double under_leader = m[k][i];
                            for (int j = 0; j < n + 1; j++) {
                                m[k][j] -= m[i][j] * under_leader;
                            }
                        }
                    } else {
                        throw zero_leader_exception(i);
                    }
                }

                for (int i = n - 1; i >= 0; --i) {
                    for (int ki = i; ki > 0; --ki) {
                        const double above_leader = m[ki - 1][i];
                        for (int j = 0; j <= n; ++j) {
                            m[ki - 1][j] -= m[i][j] * above_leader;
                        }
                    }
                }


                for (int i = 0; i < n; ++i) {
                    roots[i] = m[i][n];
                }
    }
protected:
    vector<ValTy> roots;
};

template<typename ValTy, typename IntTy, typename PolynomsTy, typename ApprFunc>
class normal_form_impl : public sole<ValTy, IntTy> {
public:
    normal_form_impl(ValTy a, ValTy b,
                     IntTy form_degree,
                     ValTy int_relative_error,
                     PolynomsTy polynomial_basis,
                     ApprFunc f):
        a(a), b(b),
        sole<ValTy, IntTy>(form_degree),
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
                printf("%.6f %.6f %.6f\n", x, exact, approximation);
                //cout << "N = " << N << " exact " << exact << " apprx " << result_n  << endl;
            }
            deviation = sqrt(deviation_sum / (normal_form.get_form_degree() + 1));
            cout << "N = " << normal_form.get_form_degree() << " dev = " << deviation << endl;
            ++normal_form;
        } while (deviation > error);

        return approximation;
    }
private:
    ValTy error;
    ValTy a;
    ValTy b;
    IntTy npoints;
    ApprFunc f;
    //vector<
    normal_form_impl<ValTy, IntTy, PolynomsTy, ApprFunc> normal_form;
    //_buffer;
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
    auto f = [&](double x) {return 1 / x - 0.1 * x * x * sin(2 * x); };

	double a = 2, b = 11;
    int N = 10;
    int n = 9;
    double epsilon = 0.01;
    auto approximator = make_lsd_approximator(a, b, epsilon, n, N, 0.00001, chebyshev, f);
    //auto approximator = make_lsd_approximator(a, b, epsilon, n, N, 0.01, algebraic_polynomial, f);
    approximator();

	return 0;
}
