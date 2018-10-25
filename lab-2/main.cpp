#define _USE_MATH_DEFINES

#include <iostream>
#include <iomanip>
#include <tuple>
#include <cmath>
#include "ranges.h"
#include "isolator.h"

using namespace std;

typedef double value_type;
typedef value_type error_type;
typedef int integer_type;
typedef edges<value_type> edges_type;

template<typename FunctionType, typename DerivativeType, typename IntegerType, typename RootValueType, typename DeltaType>
class root_approximator {
public:
    root_approximator(FunctionType f, DerivativeType df, DerivativeType ddf) : f(f), df(df), ddf(ddf) {}

    const FunctionType& get_f() const { return f; }

    const DerivativeType& get_df() const { return df; }

	const DerivativeType& get_ddf() const { return ddf; }

    const RootValueType& get_root_value() const { return root_value; }

    const IntegerType& get_nterms() const { return nterms; }

    const DeltaType& get_delta() const { return delta; }

private:
    FunctionType f;
    DerivativeType df;
	DerivativeType ddf;
protected:
    RootValueType root_value;
    IntegerType nterms;
    DeltaType delta;
};


template<typename FunctionType,
         typename DerivativeType,
         typename IntegerType,
         typename RootValueType,
         typename DeltaType = RootValueType>
class iteration_approximator : public root_approximator<FunctionType, DerivativeType, IntegerType, RootValueType, DeltaType> {
public:
    iteration_approximator(FunctionType f, DerivativeType df, DerivativeType ddf) :
        root_approximator<FunctionType, DerivativeType, IntegerType, RootValueType, DeltaType>(f, df, ddf) {}

    void operator()(const edges_type& root_range, const error_type& error) {
        const auto& a = root_range.left(), b = root_range.right();
        auto m1 = std::abs(this->get_df()(a)),
            M1 = std::abs(this->get_df()(b));

        if (m1 > M1) std::swap(m1, M1);

        auto xn = root_range.middle(), x0 = xn,
            lambda = this->get_df()(root_range.middle()) > 0 ? 1 / M1 : -1 / M1, q = 1 - m1 / M1;
        this->nterms = 0;
        do {
            x0 = xn;
            xn = x0 - lambda * this->get_f()(x0);
            ++this->nterms;
        } while (std::abs(x0 - xn) > error * (1 - q) / q);

        this->delta = std::abs(x0 - xn) * q / (1 - q);
        this->root_value = xn;
    }
};


template<typename FunctionType, typename DerivativeType, typename IntegerType, typename RootValueType, typename DeltaType = RootValueType>
class secant_approximator : public root_approximator<FunctionType, DerivativeType, IntegerType, RootValueType, DeltaType> {
public:
    secant_approximator(FunctionType f, DerivativeType df, DerivativeType ddf) :
        root_approximator<FunctionType, DerivativeType, IntegerType, RootValueType, DeltaType>(f, df, ddf) {}

    void operator()(const edges_type& root_range, const error_type& error) {
		const auto& a = root_range.left(), b = root_range.right();
		value_type m1 = std::abs(this->get_df()(a)),
			M1 = std::abs(this->get_df()(b)),
			m = m1 > M1 ? M1 : m1, c, x0;

		if (this->get_f()(a) * this->get_ddf()(a) > 0) {
			c = a;
			x0 = b;
		}
		else {
			c = b;
			x0 = a;
		}

		auto xn = x0;
		this->nterms = 0;
		do {
			x0 = xn;
			xn = x0 - this->get_f()(x0) / (this->get_f()(x0) - this->get_f()(c)) * (x0 - c);
			++this->nterms;
		} while (error < std::abs(this->get_f()(xn) / m1));

		this->delta = std::abs(this->get_f()(xn) / m1);
		this->root_value = xn;
    }
};

template<typename FunctionType,
         typename DerivativeType,
         typename IntegerType = int,
         typename RootValueType = double, typename DeltaType = RootValueType>
secant_approximator<FunctionType, DerivativeType, IntegerType, RootValueType, DeltaType>
secant_approximate(FunctionType f, DerivativeType df, DerivativeType ddf) {
    return { f, df, ddf };
}

template<typename FunctionType,
         typename DerivativeType,
         typename IntegerType = int,
         typename RootValueType = double, typename DeltaType = RootValueType>
iteration_approximator<FunctionType, DerivativeType, IntegerType, RootValueType, DeltaType>
iteration_approximate(FunctionType f, DerivativeType df, DerivativeType ddf) {
    return { f, df, ddf };
}

template<class RootApproximatorType, class ErrorRangeType>
class root_table_printer : public RootApproximatorType {
public:
    root_table_printer(std::ostream *stream, RootApproximatorType approximator,
                       ErrorRangeType error_range) :
        RootApproximatorType(approximator),
        error_range(error_range), stream(stream),
		first_bounds_saved(false) {}

	const ErrorRangeType& get_error_range() const {
		return error_range;
	}
	const edges_type& get_first_bounds() const {
		return first_bounds;
	}
	RootApproximatorType& get_approximator() {
		return static_cast<RootApproximatorType&>(*this);
	}
    void operator()(edges_type bounds) {
        bounds.shrink([this](const edges_type& bounds) {
            return bounds.contain_root(this->get_f());
        });
		
		if (!this->first_bounds_saved) {
			this->first_bounds = bounds;
			this->first_bounds_saved = true;
		}

        *stream << "x inside " << bounds << endl;

        int precision = 3;
        print_head(*stream);
        for (auto error : error_range) {
            print_error(*stream, error);
			get_approximator()(bounds, error);
			*stream << std::setprecision(precision);
			print_root_value(*stream, this->get_root_value());
			*stream << std::setprecision(6);
            print_delta(cout, this->get_delta());
			*stream << setw(10) << this->get_nterms();
            *stream << std::endl;
            precision += 3;
        }
        *stream << std::endl;
    }
protected:
    template<typename ErrTy>
    static void print_error(std::ostream& stream, const ErrTy& error) {
        stream << std::setprecision(6) << std::setw(10) << error;
    }
    template<typename RootValTy>
    static void print_root_value(std::ostream& stream, const RootValTy& root_value) {
        stream << std::setw(25) << root_value;
    }
    template<typename DeltaTy>
    static void print_delta(std::ostream& stream, const DeltaTy& delta) {
        stream << std::setw(20) << delta;
    }
    static void print_head(std::ostream& stream) {
        stream << std::left
            << std::setw(10) << "error"
            << std::setw(25) << "root value"
            << std::setw(20) << "delta" 
			<< std::setw(10) << "nterms"
			<< std::endl;
    }

private:
	bool first_bounds_saved;
	edges_type first_bounds;
    ErrorRangeType error_range;
    std::ostream *stream;
    template<class RootApprType, class ErrorRngType>
    friend std::ostream&
    operator<<(std::ostream&,
               root_table_printer<RootApprType, ErrorRngType>&);
};

template<class RootApproximator1Type, class RootApproximator2Type, class ErrorRangeType>
class tables {
	root_table_printer<RootApproximator1Type, ErrorRangeType> tab1;
	root_table_printer<RootApproximator2Type, ErrorRangeType> tab2;
public:
	tables(RootApproximator1Type approximator1, 
		RootApproximator2Type approximator2, 
		ErrorRangeType error_range): 
		tab1(&cout, approximator1, error_range), 
		tab2(&cout, approximator2, error_range) {}

	void print_compared(ostream& stream) {
		stream << setw(10) << "error"
			<< setw(15) << "iteration"
			<< setw(15) << "secant" << std::endl;
		for (auto error : tab1.get_error_range()) {
			stream << setw(10) << error;
			tab1.get_approximator()(tab1.get_first_bounds(), error);
			tab2.get_approximator()(tab2.get_first_bounds(), error);
			stream << setw(15) << tab1.get_nterms();
			stream << setw(15) << tab2.get_nterms();
			stream << std::endl;
		}
	}
	template<class RA1Ty, class RA2Ty, class ERTy>
	friend ostream& operator<<(ostream&, tables<RA1Ty, RA2Ty, ERTy>&);
};

template<class RootApproximator1Type, class RootApproximator2Type, class ErrorRangeType>
tables<RootApproximator1Type, RootApproximator2Type, ErrorRangeType>
make_tables(RootApproximator1Type approximator1,
	RootApproximator2Type approximator2,
	ErrorRangeType error_range) {
	return { approximator1, approximator2, error_range };
}

template<class RA1Ty, class RA2Ty, class ERTy>
ostream& operator<<(ostream& stream, tables<RA1Ty, RA2Ty, ERTy>& tabs) {
	stream << "ITERATION TABLE" << endl;
	stream << tabs.tab1 << endl;
	stream << "SECANT TABLE" << endl;
	stream << tabs.tab2 << endl;
	stream << "COMPARE TABLE" << endl;
	tabs.print_compared(stream);
	return stream;
}

template<class RootApprType, class ErrorRngType>
std::ostream& operator<<(std::ostream& stream,
                         root_table_printer<RootApprType, ErrorRngType>& object) {
    for_each_root<int, double>(object);
    return stream;
}

template<class RootApproximatorType, class ErrorRangeType>
root_table_printer<RootApproximatorType, ErrorRangeType>
root_table_print(std::ostream *stream, RootApproximatorType approximator, ErrorRangeType error_range) {
    return { stream, approximator, error_range };
}



int main()
{
    auto f = [](value_type x) { return 15. * sqrt(1 + cos(x)) + 1.5 * x; };
    function<value_type(value_type)> 
		df = [](value_type x) { return 1.5 - 15. * sin(x) / (2. * sqrt(1 + cos(x))); },
		ddf = [](value_type x) { return -15. / 4. * sqrt(1 + cos(x)); };
    auto error_range = rrange(1e-2, 1e-14, greater<double>(),
                              [](const value_type& x) { return x * 1e-3; });

	auto appr1 = iteration_approximate(f, df, ddf);
	auto appr2 = secant_approximate(f, df, ddf);
	auto tt = make_tables(appr1, appr2, error_range);
	cout << tt;
    return 0;
}
