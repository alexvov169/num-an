#include <iostream>
#include <iomanip>
#include <math.h>
#include "ranges.h"

using namespace std;

typedef double value_type;
typedef int integer_type;

template<typename ApproximatorType>
class tables {
public:
    typedef ApproximatorType approximator_type;
    tables(approximator_type approximator,
           value_type second_table_error,
           value_type initial_error,
           value_type final_error,
           value_type error_step,
           value_type first_table_value,
           value_type initial_value,
           value_type final_value,
           integer_type nvalues):
        second_table_error(second_table_error),
        initial_error(initial_error),
        final_error(final_error),
        error_step(error_step),
        first_table_value(first_table_value),
        initial_value(initial_value),
        final_value(final_value),
        nvalues(nvalues),
        approximator(approximator) {}
private:
    void print_first_table(std::ostream& stream) {
        stream
            << setw(20) << "eps"
            << setw(20) << "n"
            << setw(20) << "remainder term"
            << setw(20) << "delta"
            << endl;
        value_type approximated;
        for (const auto& error : rrange(initial_error, second_table_error,
                                        [this](const value_type& x) { return x * error_step; },
                                        greater<double>())) {
            approximated = approximator(first_table_value, error); // must be called first
            stream << scientific << right
                << setw(20) << error
                << setw(20) << approximator.get_nterms()
                << setw(20) << approximator.get_remainder_term()
                << setw(20) << approximator.get_exact(first_table_value) - approximated
                << endl;
        }
        approximated = approximator(first_table_value, second_table_error); // must be called first
        stream << scientific << right
            << setw(20) << second_table_error
            << setw(20) << approximator.get_nterms()
            << setw(20) << approximator.get_remainder_term()
            << setw(20) << approximator.get_exact(first_table_value) - approximated
            << endl;
        second_table_nterms = approximator.get_nterms();
    }
    void print_second_table(std::ostream& stream) {
        stream
            << setw(20) << "x"
            << setw(20) << "remainder term"
            << setw(20) << "delta"
            << endl;
        value_type approximated;
        for (const auto& x : range(initial_value, final_value, nvalues)) {
            approximated = approximator(x, second_table_nterms);
            stream << scientific << right
                << setw(20) << x
                << setw(20) << approximator.get_remainder_term()
                << setw(20) << approximator.get_exact(x) - approximated
                << endl;
        }
    }
private:
    value_type second_table_error;
    value_type initial_error;
    value_type final_error;
    value_type error_step;
    value_type first_table_value;
    value_type initial_value;
    value_type final_value;
    integer_type nvalues;
    approximator_type approximator;

    integer_type second_table_nterms;
    template<typename T>
    friend std::ostream& operator<<(std::ostream&, tables<T>);
};

template<typename T>
inline auto make_tables(T approximator,
                        value_type second_table_error,
                        value_type initial_error,
                        value_type final_error,
                        value_type error_step,
                        value_type first_table_value,
                        value_type initial_value,
                        value_type final_value,
                        integer_type nvalues) {
    return tables<T>(approximator,
                     second_table_error, initial_error, final_error, error_step,
                     first_table_value, initial_value, final_value, nvalues);
}

template<typename T>
std::ostream& operator<<(std::ostream& stream, tables<T> object) {
    object.print_first_table(stream);
    stream << endl;
    object.print_second_table(stream);
    return stream;
}



class cosh_approximator {
public:
    value_type get_exact(value_type x) const { return cosh(x); }
    const integer_type& get_nterms() const { return nterms; }
    const integer_type& get_remainder_term() const { return this->remainder_term; }
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
            i = i; // keeps compiler silent
        }

        remainder_term = current_term;

        return sum;
    }
private:
    integer_type nterms;
    value_type remainder_term;
};



int main() {
    const value_type a = -10.8,
            b = 11.9,
            n = 10,
            x = (a + b) / 2,
            ed = 1e-3, e0 = 1e-2, e1 = 1e-14, e2 = 1e-8;

    cout << make_tables(cosh_approximator(),
                        e2, e0, e1, ed,
                        x, a, b, n) << endl;
    return 0;
}
