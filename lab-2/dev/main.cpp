#define _USE_MATH_DEFINES

#include <iostream>
#include <iomanip>
#include <tuple>
#include <cmath>
#include "ranges.h"
#include "isolator.h"
#include "root_approximator.h"
#include "tables.h"

using namespace std;

typedef double value_type;
typedef value_type error_type;
typedef int integer_type;
typedef edges<value_type> edges_type;

/*
template<typename FunctionType,
         typename DerivativeType,
         class ErrorRangeType,
         class RootIsolatorType>
class RootApproximator {
    FunctionType approximated_function;
    DerivativeType derivative;
    ErrorRangeType error_range;
    RootIsolatorType root_isolator;
    mutable value_type delta;
    mutable value_type iters;
public:
    RootApproximator(FunctionType approximated_function,
                     DerivativeType derivative,
                     ErrorRangeType error_range,
                     RootIsolatorType root_isolator):
        approximated_function(approximated_function),
        derivative(derivative),
        error_range(error_range),
        root_isolator(root_isolator) {}

    const FunctionType& get_function() const {
        return approximated_function;
    }


    // Note edges must NOT be at values where does function
    // or its derivative reaches NaN values
    void operator()(edges_type edges) const {
        auto negative = std::bind(std::negate<value_type>(),
                                  std::bind(approximated_function,
                                            std::placeholders::_1));
        //cout << setw(25) << "initial x range: " << edges << endl;
        edges.shrink([this](const edges_type& e) { return e.contain_root(approximated_function); });
        cout << setw(20) << "x = " << fixed << setprecision(4) << edges << endl;
        //cout
        //	<< setw(20) << "y = "
        //	<< '[' << approximated_function(edges.left())
        //	<< ", " << approximated_function(edges.right())
        //	<< ']' << endl;

        cout
            << setw(25) << "eps"
            << setw(25) << "root"
            << setw(25) << "delta"
            << setw(25) << "n"
            << endl;
        int precision = 2;
        for (const auto& error : error_range) {
            value_type result;
            if (derivative(edges.middle()) > 0) {
                result = eval(edges.left(), edges.right(), error, approximated_function);
            } else {
                result = eval(edges.left(), edges.right(), error, negative);
            }
            cout
                << setw(25) << scientific << setprecision(1) << error
                << setw(25) << fixed << setprecision(precision + 1) << result
                << setw(25) << scientific << setprecision(6) << delta
                << setw(25) << fixed << setprecision(0) << dec << iters
                << endl;
            precision += 3;
        }

    }
protected:
    template<typename FTy>
    value_type eval(const value_type& a, const value_type& b, const value_type& error, FTy function) const {
        auto m1 = std::abs(derivative(a)),
             M1 = std::abs(derivative(b));

        if (m1 > M1) std::swap(m1, M1);

        auto xn = (a + b) / 2, x0 = xn,
             lambda = 1 / M1, q = 1 - m1 / M1;
        iters = 0;
        do {
            x0 = xn;
            xn = x0 - lambda * function(x0);
            ++iters;
        } while (std::abs(x0 - xn) > error * (1 - q) / q);

        delta = std::abs(x0-xn) * q / (1-q);
        return xn;
    }
};
template<typename FunctionType, typename DerivativeType, typename ErrorRangeType>
RootApproximator<FunctionType, DerivativeType, ErrorRangeType>
approximate_root(FunctionType approximated_function,
                       DerivativeType derivative,
                       ErrorRangeType error_range) {
    return { approximated_function, derivative, error_range };
}


void print_root_table() {

}
*/
template<typename FunctionType, typename DerivativeType, typename IntegerType, typename RootValueType, typename DeltaType>
class root_approximator {
public:
    root_approximator(FunctionType f, DerivativeType df) : f(f), df(df) {}

    const FunctionType& get_f() const { return f; }

    const DerivativeType& get_df() const { return df; }

    const RootValueType& get_root_value() const { return root_value; }

    const IntegerType& get_nterms() const { return nterms; }

    const DeltaType& get_delta() const { return delta; }

private:
    FunctionType f;
    DerivativeType df;
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
    iteration_approximator(FunctionType f, DerivativeType df) :
        root_approximator<FunctionType, DerivativeType, IntegerType, RootValueType, DeltaType>(f, df) {}

    void operator()(const edges_type& root_range, const error_type& error) {
        const auto& a = root_range.left(), b = root_range.right();
        auto m1 = std::abs(this->get_df()(a)),
            M1 = std::abs(this->get_df()(b));

        if (m1 > M1) std::swap(m1, M1);

        auto xn = (a + b) / 2, x0 = xn,
            lambda = 1 / M1, q = 1 - m1 / M1;
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
    secant_approximator(FunctionType f, DerivativeType df) :
        root_approximator<FunctionType, DerivativeType, IntegerType, RootValueType, DeltaType>(f, df) {}

    void operator()(const edges_type& root_range, const error_type& error) {

    }
};

template<typename FunctionType,
         typename DerivativeType,
         typename IntegerType = int,
         typename RootValueType = double, typename DeltaType = RootValueType>
secant_approximator<FunctionType, DerivativeType, IntegerType, RootValueType, DeltaType>
secant_approximate(FunctionType f, DerivativeType df) {
    return { f, df };
}

template<typename FunctionType,
         typename DerivativeType,
         typename IntegerType = int,
         typename RootValueType = double, typename DeltaType = RootValueType>
iteration_approximator<FunctionType, DerivativeType, IntegerType, RootValueType, DeltaType>
iteration_approximate(FunctionType f, DerivativeType df) {
    return { f, df };
}

template<class RootApproximatorType, class ErrorRangeType>
class root_table_printer : public RootApproximatorType {
public:
    root_table_printer(std::ostream *stream, RootApproximatorType approximator,
                       ErrorRangeType error_range) :
        RootApproximatorType(approximator),
        error_range(error_range), stream(stream) {}

    void operator()(edges_type bounds) {

        bounds.shrink([this](const edges_type& bounds) {
            return bounds.contain_root(this->get_f());
        });

        cout << bounds << endl;

        int precision = 3;
        print_head(*stream);
        for (auto error : error_range) {
            print_error(*stream, error);
            RootApproximatorType::operator()(bounds, error);
            cout << std::setw(30) <<
                       std::setprecision(precision) <<
                    std::right <<
                       this->get_root_value();
            print_delta(cout, this->get_delta());
            *stream << std::endl;
            precision += 3;
        }
        *stream << std::endl;
    }
protected:
    template<typename ErrTy>
    static void print_error(std::ostream& stream, const ErrTy& error) {
        stream << std::setw(20) << error;
    }
    template<typename RootValTy>
    static void print_root_value(std::ostream& stream, const RootValTy& root_value) {
        stream << std::setw(20) << root_value;
    }
    template<typename DeltaTy>
    static void print_delta(std::ostream& stream, const DeltaTy& delta) {
        stream << std::setw(20) << delta;
    }
    static void print_head(std::ostream& stream) {
        stream
            << std::setw(20) << "error"
            << std::setw(20) << "root value"
            << std::setw(20) << "delta" << std::endl;
    }

private:
    ErrorRangeType error_range;
    std::ostream *stream;
    template<class RootApprType, class ErrorRngType>
    friend std::ostream&
    operator<<(std::ostream&,
               root_table_printer<RootApprType, ErrorRngType>&);
};


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
    auto f = [](value_type x) { return 15 * sqrt(1 + cos(x)) + 1.5 * x; };
    auto df = [](value_type x) { return 1.5 - 15 * sin(x) / (2 * sqrt(1 + cos(x))); };
    auto error_range = rrange(1e-2, 1e-14, greater<double>(),
                              [](const value_type& x) { return x * 1e-3; });
    /*//auto ddf = [](value_type x) { return - 15 * sin(x) / (2 * sqrt(1 + cos(x))); };
    //for_each_root<int, double>(approximate_root(f, df,
        //                                        rrange(1e-2, 1e-14, [](const value_type& x) { return x * 1e-3; }, greater<double>())));

    auto iteration = iteration_approximator(f, df);
    auto iteration_table = tables(iteration, error_range);
    for_each_root<int, double>([]() {
        std::cout << iteration_table << std::endl;
    });*/

    auto table = root_table_print(&cout, iteration_approximate(f, df), error_range);
    cout << table;
    return 0;
}
