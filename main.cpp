#include <iostream>
#include <iomanip>
#include <math.h>
#include "ranges.h"

using namespace std;

typedef double value_type;
typedef int integer_type;

class cosh_approximator {
private:
    integer_type nterms;
    value_type remainder_term;
public:
    const integer_type& get_nterms() const { return nterms; }
    const integer_type& get_remainder_term() const { return remainder_term; }
    value_type operator()(value_type argument, value_type error_value) {
        value_type current_term = 1, sum = 1;
        integer_type denominator = 2;

        nterms = 0;
        do {
            current_term = argument * argument * current_term / denominator / (denominator - 1);
            denominator += 2;
            sum = current_term + sum;
            ++nterms;
        } while (!(current_term < error_value));

        remainder_term = current_term;

        return sum;
    }

    value_type operator()(value_type argument, integer_type nterms) {
        value_type current_term = 1, sum = 1;
        integer_type denominator = 2;

        for (auto i : range(nterms)) {
            current_term = argument * argument * current_term / denominator / (denominator - 1);
            denominator += 2;
            sum = current_term + sum;
            i = i;
        }

        remainder_term = current_term;

        return sum;
    }
};

template<typename ApproximatorType>
class tables {
private:
    approximator_type approximator;
public:
    typedef ApproximatorType approximator_type;
    tables(approximator_type approximator, ) : approximator(approximator) {}

    void build_error_table(value_type x,
                           value_type step,
                           value_type eps_begin,
                           value_type eps_last) {
        integer_type nterms;
        value_type exact, approximated, remainder_term;
        cout
            << setw(20) << "eps"
            << setw(20) << "n"
            << setw(20) << "remainder term"
            << setw(20) << "delta"
            << endl;
        for (auto error = eps_begin; (error - eps_last) > 1e-14; error *= step) {
            exact = cosh(x);
            approximated = cosh_impl(x, error, nterms, remainder_term);
            cout << scientific << right
                << setw(20) << error
                << setw(20) << nterms
                << setw(20) << remainder_term
                << setw(20) << exact - approximated
                << endl;
        }

    }
    void build_argument_table(value_type a,
                              value_type b,
                              integer_type nterms,
                              integer_type nsteps) {
        value_type exact, approximated, remainder_term;
        cout
            << setw(20) << "x"
            << setw(20) << "remainder term"
            << setw(20) << "delta"
            << endl;
        for (const auto x : range(a, b, nsteps)) {
            exact = cosh(x);
            approximated = cosh_impl(x, nterms, remainder_term);
            cout << scientific << right
                << setw(20) << x
                << setw(20) << remainder_term
                << setw(20) << exact - approximated
                << endl;
        }
    }
private:
    value_type initial_error;
    value_type final_error;
    value_type error_step;
    value_type initial_value;
    value_type final_value;
    value_type value_step;

    friend std::ostream& operator<<(std::ostream&, tables);
};

std::ostream& operator<<(std::ostream& stream, tables object) {
    auto x = (object.initial_value + object.final_value) / 2;
    build_error_table(x, object.value_step, object.initial_error, object.final_error);
    stream << endl;
    build_argument_table(object.initial_value, object.final_value, 5, 10);
    return stream;
}

int main() {
    const value_type a = -10.8,
            b = 11.9,
            x = (a + b) / 2,
            step = 1e-3, eps_begin = 1e-2, eps_last = 1e-14;

    return 0;
}
