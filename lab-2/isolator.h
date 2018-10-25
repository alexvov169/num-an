#ifndef ISOLATOR_H
#define ISOLATOR_H

/* instead of 'obviously'))
 * f(x) = 0
 * f(x) = 15 * sqrt(1 + cos(x)) + 1.5 * x
 * 15sqrt(cos(x) + 1) + 1.5x = 0 <==> x = -10sqrt(cos(x) + 1)
 *                                      = -10sqrt(2 * cos^2(x/2))
 *                                      = -10sqrt(2) * abs(cos(x/2))
 * f'(x) = 1.5 - 15 * sin(x) / (2 * sqrt(1 + cos(x)))
 *       = 1.5 - 15 * sin(x) / (2 * sqrt(2 * cos^2(x/2))
 *       = 1.5 - 15 * sin(x) / (2 * sqrt(2) * abs(cos(x/2)))
 *       = 1.5 - 15 * 2 * sin(x/2) * cos(x/2) / (2 * sqrt(2) * abs(cos(x/2)))
 *       = (x == pi * n) ? NaN
 *                       : (cos(x/2) > 0) ? 1.5 - (15 / sqrt(2)) * sin(x/2)
 *                                        : 1.5 + (15 / sqrt(2)) * sin(x/2)
 * f'(x) = 0
 *      x = 2 * atan(1 / 7) + 2 * pi * n
 * f'(x) > 0 => -pi + 2 * pi * n < x < 2 * atan(1 / 7) + 2 * pi * n
 * f'(x) < 0 => 2 * atan(1 / 7) + 2 * pi * n < x < -pi + 2 * pi * n
 *
 */

#include <iostream>
#include <tuple>
#include <cmath>

template<typename ValueType>
class edges : public std::pair<ValueType, ValueType> {
public:
    typedef ValueType value_type;

private:
    typedef std::pair<value_type, value_type> base_type;

public:
    edges() {}
    edges(value_type first, value_type second) : base_type(first, second) {}
    edges& operator=(const edges& other) {
		this->first = other.first;
		this->second = other.second;
		return *this;
	}

    void set_left(const value_type& other_left) {
        base_type::first = other_left;
    }
    void set_right(const value_type& other_right) {
        base_type::second = other_right;
    }
    const value_type& left() const {
        return base_type::first;
    }
    const value_type& right() const {
        return base_type::second;
    }
    value_type middle() const {
        return 0.5 * (left() + right());
    }

    template<typename UnaryFunction>
    bool contain_root(UnaryFunction f) const {
        return f(left()) * f(right()) < 0;
    }


    // Modifies edges so they don't include their old 'edges'.
	// Needed for approximating methods which don't apply
	// ranges with zero or NaN derivative
    template<typename PredicateType>
    edges& shrink(PredicateType is_appropriate) {
        value_type old_left = left();
        set_left(middle());
        if (is_appropriate(*this)) {
            value_type old_right = right();
            set_right(middle());
            while (!is_appropriate(*this)) {
                //std::cout << *this << std::endl;
                set_left(right());
                set_right(old_right);
                set_right(middle());
            }
        } else {
            while (!is_appropriate(*this)) {
                //std::cout << *this << std::endl;
                set_right(left());
                set_left(old_left);
                set_left(middle());
            }
        }
        return *this;
    }

    // positive and negative stands for derivative's sign
    // in returned edges
    template<typename IntTy>
    static edges get_negative(IntTy n) {
		const value_type atan_1_div_7 = std::atan(1. / 7);
        return { atan_1_div_7 * 2 + 2 * M_PI * n, -M_PI + 2 * M_PI * (n + 1) };
    }
    template<typename IntTy>
    static edges get_positive(IntTy n) {
		const value_type atan_1_div_7 = std::atan(1. / 7);
        return { -M_PI + 2 * M_PI * n, atan_1_div_7 * 2 + 2 * M_PI * n };
    }

    template<typename ValTy>
    friend std::ostream& operator<<(std::ostream&, edges<ValTy>);
};

template<typename ValTy>
std::ostream& operator<<(std::ostream& stream, edges<ValTy> object) {
    stream << '{' << object.left() << ", " << object.right() << '}';
    return stream;
}


// using given function's feature
// we search roots in neighbour(!) edges

// note that this algorithm does not protect
// from endless search when f has no roots
template<typename IntegerType,
         typename ValueType,
         class FunctorType>
void for_each_root(FunctorType& function) {
    IntegerType n{};
    typedef edges<ValueType> edges_type;
    const auto& f = function.get_f();

    // conditions below were intentionally left inefficient
    while (!edges_type::get_positive(n).contain_root(f) &&
           !edges_type::get_negative(n).contain_root(f) &&
           !edges_type::get_positive(-n).contain_root(f) &&
           !edges_type::get_negative(-n).contain_root(f)) {
        ++n;
    }

    edges_type positive_edges, negative_edges;
    bool positive_direction_ended = false,
         negative_direction_ended = false;
    if (n == 0) {
        positive_edges = edges_type::get_positive(n),
        negative_edges = edges_type::get_negative(n);
        if (positive_edges.contain_root(f)) {
            function(positive_edges);
        } else {
            negative_direction_ended = true;
        }
        if (negative_edges.contain_root(f)) {
            function(negative_edges);
        } else {
            positive_direction_ended = true;
        }
        ++n;
    }
    while (!positive_direction_ended || !negative_direction_ended) {
        // positive direction
        if (!positive_direction_ended) {
            positive_edges = edges_type::get_positive(n),
            negative_edges = edges_type::get_negative(n);
            if (positive_edges.contain_root(f)) {
                function(positive_edges);
                if (negative_edges.contain_root(f)) {
                    function(negative_edges);
                } else {
                    positive_direction_ended = true;
                }
            } else {
                positive_direction_ended = true;
            }
        }

        // negative direction
        if (!negative_direction_ended) {
            positive_edges = edges_type::get_positive(-n),
            negative_edges = edges_type::get_negative(-n);
            if (negative_edges.contain_root(f)) {
                function(negative_edges);
                if (positive_edges.contain_root(f)) {
                    function(positive_edges);
                } else {
                    negative_direction_ended = true;
                }
            } else {
                negative_direction_ended = true;
            }
        }
        ++n;
    }
}




#endif // ISOLATOR_H
