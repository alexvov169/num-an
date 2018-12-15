#pragma once

#include <vector>

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

template<typename ValTy, typename IntTy>
class sole : protected std::vector<std::vector<ValTy> > {
public:
    sole(IntTy n) : std::vector<std::vector<ValTy> >(n) {}

    const std::vector<ValTy>& get_roots() const { return roots; }

    void complete_pivoting() {
        const IntTy n = this->size();
        std::vector<IntTy> indices(n);
        for (IntTy m = 0; m < n; ++m) {
            ValTy maximum = (*this)[m][m];
            IntTy imax = m;
            IntTy jmax = m;
            for (IntTy i = m+1; i < n; ++i) {
                for (IntTy j = m+1; j < n; ++j) {
                    if (maximum < (*this)[i][j]) {
                        maximum = (*this)[i][j];
                        imax = i;
                        jmax = j;
                    }
                }
            }
            indices[m] = jmax;

            if (m != imax) {
                std::swap((*this)[m], (*this)[imax]);
            }
            if (m != jmax) {
                for (IntTy i = 0; i < n; ++i) {
                    std::swap((*this)[i][m], (*this)[i][jmax]);
                }
            }
        }
        gaussian_elimination();
        for (IntTy i = 0; i < n; ++i) {
            std::swap(roots[i], roots[indices[i]]);
        }
    }
    void gaussian_elimination() {
        roots.resize(this->size());
        auto& m = *this;
                const int n = this->size();

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
    }
protected:
    std::vector<ValTy> roots;
};
