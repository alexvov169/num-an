#pragma once

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


template<typename ValTy>
class chebyshov_polynomials : public polynomial<ValTy> {
public:
    chebyshov_polynomials(ValTy x) : polynomial<ValTy>(x) {}

    class iterator {
        const chebyshov_polynomials *that;
        ValTy previous;
        ValTy current;
    public:
        iterator(const chebyshov_polynomials *that, ValTy previous, ValTy current):
            that(that), previous(previous), current(current) {}

        const ValTy& operator*() const { return previous; }
        iterator& operator++() {
            ValTy temp = previous;
            previous = current;
            current = 2 * that->get_x() * current - temp;
            return *this;
        }
    };

    iterator begin() const {
        return { this, 1, this->get_x() };
    }
};

template<typename ValTy>
chebyshov_polynomials<ValTy>
make_chebyshov_polynomials(ValTy x) {
    return {x};
}

template<typename ValTy, typename IntTy>
ValTy chebyshov_polynomial(IntTy n, ValTy x) {
    auto t = make_chebyshov_polynomials(x);
    auto i = t.begin();
    while (n--) { ++i; }
    return *i;
}

template<typename ValTy, typename IntTy>
ValTy algebraic_polynomial(IntTy n, ValTy x) {
    auto t = make_algebraic_polynomials(x);
    auto i = t.begin();
    while (n--) { ++i; }
    return *i;
}
