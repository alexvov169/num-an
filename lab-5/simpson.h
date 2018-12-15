#pragma once

#include <cmath>

template<typename ValTy, typename Integrand>
class simpson_integrate_impl {
public:
	simpson_integrate_impl(ValTy a, ValTy b, Integrand integrand) :
		integrand(integrand), a(a), b(b) {}

	const ValTy& get_step() const { return h; }

	template<typename IntTy>
	ValTy common_step(IntTy n) const {
		return (b - a) / n;
	}

	ValTy relative_runge(ValTy eps) {
		int i;
        ValTy sig1 = 0;
        ValTy sig2 = 0;
        int n = int(1/sqrt(eps));
        ValTy y0 = integrand(a), yn = integrand(b);
        h = common_step(n);
        for (i = 1; i < n; i += 2) {
            sig2 += integrand(a + h * i);
		}
        for (i = 2; i < n; i += 2) {
            sig1 += integrand(a + h * i);
        }
        ValTy curr_int, prev_int = h / 3 * (4 * sig2 + 2 * sig1 + y0 + yn),
                curr_even = sig2 + sig1;
        int curr_n = 2 * n;
		do {
            h = (b - a) / curr_n;
            ValTy curr_odd = 0;
            for (int i = 1; i < curr_n; i += 2) {
                curr_odd += integrand(a + h * i);
			}
            curr_int = h / 3 * (4 * curr_odd + 2 * curr_even + y0 + yn);
            prev_int = curr_int;
            curr_even = curr_odd + curr_even;
            curr_n *= 2;
        } while (!(((abs((curr_int - prev_int) / curr_int)) < eps)));
        return curr_int;
	}

protected:
	Integrand integrand;
	ValTy a;
	ValTy b;
	ValTy h;
};

template<typename ValTy, typename Integrand>
simpson_integrate_impl<ValTy, Integrand>
simpson_integrate(ValTy a, ValTy b, Integrand integrand) {
	return { a, b, integrand };
}
