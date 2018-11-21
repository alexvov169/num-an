#ifndef TRAPEZIUM_H
#define TRAPEZIUM_H

#include <cmath>

template<typename ValTy, typename Integrand>
class trapezium_integrate_impl {
public:
    const int r = 2;
    const int R_denom = 3; // 2^r - 1

    trapezium_integrate_impl(Integrand integrand, ValTy a, ValTy b):
        integrand(integrand), a(a), b(b) {}

    ValTy analytic_step(ValTy max_ddf, ValTy eps) const {
        return std::sqrt((12 * eps) / std::abs((b - a)*max_ddf));
    }

    template<typename IntTy>
    ValTy common_step(IntTy n) const {
        return (b - a) / n;
    }

    ValTy analytical(ValTy epsilon, ValTy max_ddf) {
        h = analytic_step(max_ddf, epsilon);
        int n = int((b-a) / h);
        h = (b-a)/n;
        return value(h, sum(h, int(0), n));
    }

    template<typename IntTy>
    ValTy runge(ValTy epsilon, IntTy n) {
        h = common_step(n);
        ValTy sum_n = sum(h, 1, n),
              sum_2n = 0,
              value_n, value_2n = 0;
        do {
            sum_n = sum(h, 1, n),
            sum_2n = sum(h/2, 1, 2*n);
            value_n = value(h, sum_n);
            value_2n = value(h/2, sum_2n);
            n *= 2;
            h = common_step(n);
        } while (std::abs(value_n-value_2n) / R_denom > epsilon);

        return value_2n;
    }

    ValTy value(ValTy h, ValTy sum) {
        return h * (integrand(a) / 2 + integrand(b) / 2 + sum);
    }

    template<typename IntTy>
    ValTy sum(ValTy h, IntTy begin, IntTy end) {
        ValTy result = 0;
        for (auto i = begin; i != end; ++i) {
            result = integrand(a + h*i) + result;
        }
        return result;
    }

    const ValTy& get_step() const { return h; }

protected:
    Integrand integrand;
    ValTy a;
    ValTy b;
    ValTy h;
};


template<typename ValTy, typename Integrand, typename Primitive>
class trapezium_integrate_with_exact_impl : public trapezium_integrate_impl<ValTy, Integrand> {
public:
    trapezium_integrate_with_exact_impl(Integrand integrand, Primitive primitive, ValTy a, ValTy b):
        trapezium_integrate_impl<ValTy, Integrand>(integrand, a, b), primitive(primitive) {}
    ValTy exact() const {
        return primitive(this->b) - primitive(this->a);
    }
protected:
    Primitive primitive;
};

template<typename ValTy, typename Integrand>
trapezium_integrate_impl<ValTy, Integrand>
trapezium_integrate(Integrand integrand, ValTy a, ValTy b) {
    return { integrand, a, b };
}

template<typename ValTy, typename Integrand, typename Primitive>
trapezium_integrate_with_exact_impl<ValTy, Integrand, Primitive>
trapezium_integrate(Integrand integrand, Primitive primitive, ValTy a, ValTy b) {
    return { integrand, primitive, a, b };
}


#endif // TRAPEZIUM_H
