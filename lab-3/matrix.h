#ifndef MATRIX_H
#define MATRIX_H

#include <cmath>
#include <iostream>
#include <iomanip>
#include <exception>

// simple pointer-to-vector (matrix) wrapper
template<typename T, int nrows, int ncols>
class matrix {
protected:
    typedef T (*container_ptr)[ncols];
    typedef const T (*const_container_ptr)[ncols];
    typedef T container_type[nrows][ncols];
public:
    explicit matrix(container_ptr data) : m(data) {}

    container_type& c_mtr() {
        return reinterpret_cast<container_type&>(*m);
    }
    const container_type& c_mtr() const {
        return reinterpret_cast<const container_type&>(*m);
    }

private:
    container_ptr m;
};

template<typename T, int nrows, int ncols>
std::ostream& operator<<(std::ostream& stream, const matrix<T, nrows, ncols>& object) {
    for (const auto& row : object.c_mtr()) {
        stream << '(';
        for (const auto& i : row) {
            stream << std::setw(15) << i;
        }
        stream << ")\n";
    }
    return stream;
}

template<typename IndexType>
class zero_leader_exception_impl : public std::exception {
public:
    zero_leader_exception_impl(IndexType i) : i(i) {}
    virtual const char *what() const noexcept override { return "zero leader"; }
    IndexType i;
};
template<typename IndexType>
zero_leader_exception_impl<IndexType>
zero_leader_exception(IndexType i) {
    return zero_leader_exception_impl<IndexType>(i);
}

template<typename T, int nequations>
class sole : public matrix<T, nequations, nequations + 1> {
    typedef matrix<T, nequations, nequations + 1> base_type;
    typedef typename base_type::container_type container_type;

    typedef T roots_container[nequations];
public:
    explicit sole(typename base_type::container_ptr data) : base_type(data) {}

    const roots_container& gaussian_elimination() {
        container_type& m = this->c_mtr();
        const int& n = nequations;

        for (int i = 0; i < n; ++i) {
            if (m[i][i] != 0) {
                const double leader = m[i][i];
                for (int j = 0; j < n + 1; ++j) {
                    m[i][j] /= leader;
                }
                for (int k = i + 1; k < n; ++k) {
                    const double under_leader = m[k][i];
                    for (int j = 0; j < n + 1; j++) {
                        m[k][j] -= m[i][j] * under_leader;
                    }
                }
            } else {
                throw zero_leader_exception(i);
            }
        }

        for (int i = n - 1; i >= 0; --i) {
            for (int ki = i; ki > 0; --ki) {
                const double above_leader = m[ki - 1][i];
                for (int j = 0; j <= n; ++j) {
                    m[ki - 1][j] -= m[i][j] * above_leader;
                }
            }
        }


        for (int i = 0; i < n; ++i) {
            roots[i] = m[i][n];
        }

        return roots;
    }

    const roots_container& seidel_iteration(double eps) {
        container_type& m = this->c_mtr();
        const int n = nequations;

        double previous_roots[n];

        for (int i = 0; i < n; ++i) {
            double leader = m[i][i];
            for (int j = 0; j < n + 1; ++j) {
                m[i][j] = m[i][j] / leader;
            }
            roots[i] = m[i][n];
        }            

        double q = 0;
        for (int i = 0; i < n; ++i){
            double max = 0;
            for (int j = 0; j < i; ++j) {
                max += fabs(m[i][j]);
            }
            for (int j = i + 1; j < n; ++j) {
                max += fabs(m[i][j]);
            }
            if (max > q) {
                q = max;
            }
        }

        int k = 1;
        double norma;
        do {
            for (int i = 0; i < n; i++) {
                previous_roots[i] = roots[i];
            }
            for (int i = 0; i < n; i++) {
                roots[i] = m[i][n];
                for (int j = 0; j < i; ++j) {
                    roots[i] -= m[i][j] * roots[j];
                }
                for (int j = i + 1; j < n; ++j) {
                    roots[i] -= m[i][j] * roots[j];
                }
            }

            norma = fabs(roots[0] - previous_roots[0]);
            for (int i = 1; i < n; i++) {
                double delta = fabs(roots[i] - previous_roots[i]);
                if (delta > norma) norma = delta;
            }
            k++;
        } while (norma > eps * (1 - q) / q);

        iterations = k;

        return roots;
    }
    const int& get_iterations() const { return iterations; }
private:
    roots_container roots;
    int iterations;
};

#endif // MATRIX_H
