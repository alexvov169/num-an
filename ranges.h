#ifndef RANGES_H
#define RANGES_H

#include <functional>

template<typename IndexType>
class basic_range {
public:
    typedef IndexType index_type;
    basic_range(index_type n) : n(n) {}
    class iterator {
    public:
        iterator(index_type n = 0) : current(n) {}
        inline bool operator!=(const iterator& other) const {
            return this->current < other.current;
        }
        inline const index_type& operator*() const {
            return current;
        }
        inline iterator& operator++() {
            ++current;
            return *this;
        }
    protected:
        index_type current;
    };
    iterator begin() const {
        return iterator(0);
    }
    iterator end() const {
        return iterator(n);
    }
private:
    index_type n;
};

template<typename ValueType, typename ModifierType>
class recursive_iterable {
public:
    typedef ValueType value_type;
    typedef ModifierType modifier_type;

    recursive_iterable(value_type a, value_type b,
                       modifier_type modifier):
        initial(a), final(b), modifier(modifier) {}

    class iterator : public basic_range<value_type>::iterator {
        typedef typename basic_range<value_type>::iterator base_type;
    public:
        iterator(const recursive_iterable *that, value_type initial):
            base_type(initial), that(that) {}
        inline iterator& operator++() {
            base_type::current = that->modifier(base_type::current);
            return *this;
        }
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
    value_type initial;
    value_type final;
    modifier_type modifier;
};

template<typename ValueType, class RecursiveIterableType>
class iterative_iterable : public RecursiveIterableType {
    typedef RecursiveIterableType base_type;
public:
    typedef ValueType value_type;
    iterative_iterable(value_type a, value_type b, base_type iterable):
        base_type(iterable),
        step((b - a) / (n - 1)), initial(a) {}
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
            base_type::operator++();
            this->current = that->initial + that->step * base_type::operator*();
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

template<typename IndexType>
inline basic_range<IndexType> range(IndexType n) {
    return basic_range<IndexType>(n);
}

template<typename IntegerType, typename ModifierType>
inline recursive_iterable<IntegerType> rrange(IntegerType a, IntegerType b,
                                              ModifierType modifier) {
    return recursive_iterable<IntegerType, ModifierType>(a, b, modifier);
}

template<typename IntegerType, class RecursiveIterableType>
inline iterative_iterable<IntegerType> irange(IntegerType a, IntegerType b,
                                              RecursiveIterableType iterable) {
    return iterative_iterable<IntegerType, RecursiveIterableType>(a, b, iterable);
}
#endif // RANGES_H
