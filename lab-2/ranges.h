#ifndef RANGES_H
#define RANGES_H

#include <iostream>
#include <functional>

/*
 * Numeric iterable with following formula:
 * i = foo(i)
 * */

template<typename ValueType,
         typename NotEqualFunctorType,
         typename EvalNextFunctorType,
         typename EvalPrevFunctorType>
class recursive_iterable {
public:
    recursive_iterable(ValueType initial, ValueType final,
                       NotEqualFunctorType not_equal,
                       EvalPrevFunctorType previous,
                       EvalNextFunctorType next):
        initial(initial), final(final),
        not_equal(not_equal),
        previous(previous), next(next) {}

    class iterator {
    public:
        iterator(const recursive_iterable *that, ValueType initial):
            current(initial), that(that) {}
        inline iterator& operator++() {
            current = that->next(current);
            return *this;
        }
        inline iterator& operator--() {
            current = that->previous(current);
            return *this;
        }
        inline bool operator!=(const iterator& other) const {
            return that->not_equal(this->current, other.current);
        }
        inline const ValueType& operator*() const {
            return current;
        }
    protected:
        ValueType current;
    private:
        const recursive_iterable *that;
    };
    iterator begin() const {
        return iterator(this, initial);
    }
    iterator end() const {
        return iterator(this, final);
    }
private:
    ValueType initial;
    ValueType final;
    NotEqualFunctorType not_equal;
    EvalPrevFunctorType previous;
    EvalNextFunctorType next;
};

template<typename ValueType,
         typename NotEqualFunctorType,
         typename EvalNextFunctorType>
class recursive_iterable<ValueType, NotEqualFunctorType, EvalNextFunctorType, void> {
public:
    recursive_iterable(ValueType initial, ValueType final,
                       NotEqualFunctorType not_equal,
                       EvalNextFunctorType next):
        initial(initial), final(final),
        not_equal(not_equal),
        next(next) {}

    class iterator {
    public:
        iterator(const recursive_iterable *that, ValueType initial):
            current(initial), that(that) {}
        inline iterator& operator++() {
            current = that->next(current);
            return *this;
        }
        inline bool operator!=(const iterator& other) const {
            return that->not_equal(this->current, other.current);
        }
        inline const ValueType& operator*() const {
            return current;
        }
    protected:
        ValueType current;
    private:
        const recursive_iterable *that;
    };
    iterator begin() const {
        return iterator(this, initial);
    }
    iterator end() const {
        return iterator(this, final);
    }
private:
    ValueType initial;
    ValueType final;
    NotEqualFunctorType not_equal;
    EvalNextFunctorType next;
};

/*
 * Numeric iterable with following formula:
 * i = a + j * h
 * where a -- initial constant
 *       j -- passed iterable or if passed number of iterations
 *            then basic_iterable is used
 *       h -- passed step or if passed number of iterations and final value
 *            then step = (final - initial) / (n - 1) is used
 * */
template<typename ValueType, class RecursiveIterableType>
class iterative_iterable : public RecursiveIterableType {
    typedef RecursiveIterableType base_type;
public:
    typedef ValueType value_type;
    iterative_iterable(value_type a, value_type step,
                       base_type iterable):
        base_type(iterable),
        step(step),
        initial(a) {}

    class iterator : public base_type::iterator {
        typedef typename base_type::iterator base_iterator_type;
    public:
        explicit iterator(const iterative_iterable *that, base_iterator_type base):
            base_iterator_type(base),
            that(that), current(that->initial) {}

        const value_type& operator*() const {
            return this->current;
        }
        iterator& operator++() {
            base_iterator_type::operator++();
            this->current = that->initial +
                    that->step * base_iterator_type::operator*();
            return *this;
        }
    protected:
        const iterative_iterable *that;
        value_type current;
    };
    iterator begin() const {
        return iterator(this, base_type::begin());
    }
    iterator end() const {
        return iterator(this, base_type::end());
    }
private:
    value_type step;
    value_type initial;
    friend class iterator;
};

template<typename IntegerType,
         typename NotEqualFunctorType,
         typename EvalNextFunctorType>
recursive_iterable<IntegerType, NotEqualFunctorType, EvalNextFunctorType, void>
rrange(IntegerType a,
       IntegerType b,
       NotEqualFunctorType not_equal,
       EvalNextFunctorType next) {
    return { a, b, not_equal, next };
}

template<typename IntegerType,
         typename NotEqualFunctorType,
         typename EvalNextFunctorType,
         typename EvalPrevFunctorType>
recursive_iterable<IntegerType, NotEqualFunctorType, EvalNextFunctorType, EvalPrevFunctorType>
rrange(IntegerType a,
       IntegerType b,
       NotEqualFunctorType not_equal,
       EvalPrevFunctorType previous,
       EvalNextFunctorType next) {
    return { a, b, not_equal, previous, next };
}
template<typename IntegerType>
auto 
rrange(IntegerType a, IntegerType b, IntegerType h)
-> decltype (rrange(a, b, std::less<IntegerType>(),
                    std::bind(std::plus<IntegerType>(), h, std::placeholders::_1))) {
    return rrange(a, b, std::less<IntegerType>(),
                  std::bind(std::plus<IntegerType>(), h, std::placeholders::_1));
}

template<typename IntegerType, class RecursiveIterableType>
auto irange(IntegerType a, IntegerType b, RecursiveIterableType iterable)
-> decltype (iterative_iterable<IntegerType, RecursiveIterableType>(a, b, iterable)) {
    return iterative_iterable<IntegerType, RecursiveIterableType>(a, b, iterable);
}

template<typename IndexType>
auto range(IndexType n)
-> decltype (rrange(0, n, 1)) {
    return rrange(0, n, 1);
}

template<typename ValueType, typename IndexType>
auto range(ValueType a, ValueType b, IndexType n) 
-> decltype (irange(a, (b - a) / (n - 1), range(n))) {
	return irange(a, (b - a) / (n - 1), range(n));
}

#endif // RANGES_H
