#ifndef MATRIX_H
#define MATRIX_H

#include <algorithm>
#include <iostream>


// simple pointer-to-vector (matrix) wrapper
template<typename T, size_t nrows, size_t ncols>
class matrix {
protected:
    typedef T (*container_ptr)[ncols];
    typedef const T (*const_container_ptr)[ncols];
    typedef T container_type[nrows][ncols];
public:
    explicit matrix(container_ptr data) : m(data) {}

    // needed for compatibility with for-each loops
    container_type& c_mtr() {
        return reinterpret_cast<container_type&>(*m);
    }
    const container_type& c_mtr() const {
        return reinterpret_cast<const container_type&>(*m);
    }
    operator container_type&() {
        return reinterpret_cast<container_type&>(*m);
    }
    operator const container_type&() const {
        return reinterpret_cast<const container_type&>(*m);
    }
private:
    container_ptr m;
};

template<typename T, size_t nrows, size_t ncols>
std::ostream& operator<<(std::ostream& stream, const matrix<T, nrows, ncols>& object) {
    for (const auto& row : object.c_mtr()) {
    //for (const auto& row : object) {
        stream << '[';
        for (const auto& i : row) {
            stream << i << ", ";
        }
        stream << "],\n";
    }
    return stream;
}

template<typename T, size_t nequations>
class sole : public matrix<T, nequations, nequations + 1> {
    typedef matrix<T, nequations, nequations + 1> base_type;
public:
    sole(typename base_type::container_ptr data) : base_type(data) {}
    void gaussian_elimination() {

    }
    void seidel_iteration() {

    }
};

#endif // MATRIX_H
