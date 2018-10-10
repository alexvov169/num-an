#include <iostream>
#include <iomanip>
#include <math.h>
#include "ranges.h"
#include "tables.h"

using namespace std;

typedef double value_type;
typedef int integer_type;

class cosh_approximator {
public:
	value_type operator()(value_type argument, value_type error_value) const {
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

	value_type operator()(value_type argument, integer_type nterms) const {
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
protected:
	value_type get_exact(value_type x) const { return cosh(x); }
	const integer_type& get_nterms() const { return nterms; }
	const value_type& get_remainder_term() const { return this->remainder_term; }
private:
	mutable integer_type nterms;
	mutable value_type remainder_term;
	friend class tables<cosh_approximator, value_type, integer_type>;
};


int main() {
    const value_type a = -10.8,
            b = 11.9,
            x = (a + b) / 2,
            ed = 1e-3, e0 = 1e-2, e1 = 1e-14, e2 = 1e-8;

    cout << make_tables(cosh_approximator(),
                        e2, e0, e1, ed,
                        x, a, b, 10) << endl;
    return 0;
}
