#define _USE_MATH_DEFINES

#include <math.h>
#include <iomanip>
#include <iostream>

#include "simpson.h"
#include "matrix.h"

using namespace std;

double f(double x){
    return 1/x - .1*x*x * sin(2*x);
}

template<typename ValTy, typename IntTy>
class legendre_iterator_impl {
public:
    legendre_iterator_impl(ValTy previous, ValTy current, IntTy n):
        previous(previous), current(previous), n(n), x(current) {}

    legendre_iterator_impl& operator++() {
        next = (2*n + 1) / (n+1) * x * current - n / (n+1) * previous;
        previous = current;
        current = next;
        ++n;
        return *this;
    }

    const ValTy& operator*() const { return current; }

private:
    ValTy previous;
    ValTy current;
    ValTy next;
    IntTy n;
    ValTy x;
};

template<typename ValTy, typename IntTy>
legendre_iterator_impl<ValTy, IntTy>
make_legendre(ValTy previous, ValTy current, IntTy n) {
    return { previous, current, n };
}

double L(int n, double x){
    double Ln_next, Ln = x, Ln_first = 1;

    if (n == 0) return Ln_first;
    if (n == 1) return Ln;

    int i = 1;
    while (i < n){
        Ln_next = (1.0 * (2 * i + 1) / (i + 1)) * x * Ln - (1.0 * i / (i + 1)) * Ln_first;
        Ln_first = Ln;
        Ln = Ln_next;
        ++i;
    }
    return Ln;
}

double func(int N, double t, int k1, int k2, double a, double b){
    double x = (2 * t - a - b) / (b - a);
    if (k2==N)
        return  f(t)*L(k1, x);
    else
        return  L(k1, x)*L(k2, x);
}

double simpson(int k1, int k2, double a, double b, int n, int N) {
    int i;
    double h;
    double sig1 = 0;
    double sig2 = 0;
    double y0 = func(N, a, k1, k2, a, b), yn = func(N, b, k1, k2, a, b);

    h = (b - a) / n;
    for (i = 1; i < n; i++) {
        if (i % 2 == 0)
            sig2 += func(N, a + i*h, k1, k2, a, b);
        else
            sig1 += func(N, a + i*h, k1, k2, a, b);
    }

    return h / 3 * (2 * sig2 + 4 * sig1 + y0 + yn);
}

double Integral(int k1, int k2, double a, double b, int N){
    double r, In, I2n;
    double eps = 1e-8;
    int n = (int)ceil((b - a) / sqrt(sqrt(eps)));

    In = simpson(k1, k2, a, b, n, N);
    I2n = simpson(k1, k2, a, b, 2 * n, N);
    r = fabs(In - I2n) / 15;

    while (r > eps) {
        In = I2n;
        n *= 2;
        I2n = simpson(k1, k2, a, b, n, N);
        r = fabs(In - I2n) / 15;
    }
    return I2n;
}

double* scheme_selection_main_element(double **matrix, int N){
    int i, j, line, row, k, line_new = 0;
    double R, max, *M = new double[N];
    double *res = new double[N];
    double ** new_matrix = new double*[N];
    for (int i = 0; i < N + 1; i++)
        new_matrix[i] = new double[N + 1];

    for (k = 0; k < N - 1; k++){

        max = 0;
        for (i = 0; i < N; i++)
            for (j = 0; j < N; j++)
                if (fabs(matrix[i][j]) > max){
            max = fabs(matrix[i][j]);
            line = i;
            row = j;
                }

        for (j = 0; j < N + 1; j++)
            new_matrix[line_new][j] = matrix[line][j];
        line_new++;

        for (i = 0; i < N; i++)
            M[i] = -(matrix[i][row] / max);

        for (i = 0; i < N; i++)
            for (j = 0; j < N + 1; j++)
                if (i != line)
                    matrix[i][j] += matrix[line][j] * M[i];

        for (j = 0; j < N + 1; j++)
            matrix[line][j] = 0;

        for (i = 0; i < N; i++)
            matrix[i][row] = 0;
    }

    for (i = 0; i < N; i++)
        if (matrix[i][N] != 0)
            line = i;

    for (j = 0; j < N + 1; j++)
        new_matrix[line_new][j] = matrix[line][j];

    for (j = 0; j < N; j++)
        res[j] = 0;

    for (i = N - 1; i >= 0; i--){
        R = new_matrix[i][N];
        for (j = 0; j < N; j++)
            if ((new_matrix[i][j] != 0) && (res[j] != 0))
                R -= new_matrix[i][j] * res[j];
        for (j = 0; j < N; j++)
            if ((new_matrix[i][j] != 0) && (res[j] == 0))
                res[j] = R / new_matrix[i][j];
    }

    return res;
}

double * make_A(double a, double b, int N){
    int i, j;
    double *result = new double[N];
    double ** A = new double*[N];
    for (int i = 0; i < N + 1; i++)
        A[i] = new double[N + 1];

    for (i = 0; i < N; i++) {
        for (j = 0; j < N + 1; j++) {
            A[i][j] = A[j][i] = Integral(i, j, a, b, N);
        }
        A[i][N] = Integral(i, N, a, b, N);
    }

    result = scheme_selection_main_element(A, N);
    return result;

}

double Pm(double x, double a, double b, int n){
    double *A = make_A(a, b, n);
    double  t = (2 * x - b - a) / (b - a);
    double y = 0;
    auto L = make_legendre(1., t, 0);
    for (int j = 0; j < n; ++j) {
        y += A[j] * *L;//(j, t);
        ++L;
    }
    return y;
}

int main(){
    double a = 20, b = 30;
    double eps = 0.01;

    int k = (int)(b - a) * 2;//кількість точок, які обчислюються на заданому інтервалі
    double y, x, step = (b - a) / k;
    int n = 2;
    double Dev; // Deviation - середнє квадратичне відхилення

    do{
        n += 1;
        Dev = 0;
        x = a;
        for (int i = 0; i <= k; ++i){
            y = f(x) - Pm(x, a, b, n);
            Dev += y*y;
            x += step;
        }
        Dev = sqrt(Dev / (k+1));
        cout << "Deviation[" << n << "] = " << Dev << endl;
    } while (Dev > eps);

    cout << "N = " << n << endl << endl;

    cout << "   X          P(x)" << endl;
    for (double i = a; i <= b; i += step){
        double h = Pm(i, a, b, n);
        cout << setw(7) << fixed << setprecision(2) << i;
        cout << setw(15) << fixed << setprecision(8) << h << endl;
    }

    return 0;
}
