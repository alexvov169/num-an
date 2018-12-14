#include <math.h>
#include <iostream>
#include <fstream>
#include <functional>

using namespace std;

double simpson(double a, double b, double eps, function<double(double)> y) {
    int i;
    double h;
    double sig1 = 0;
    double sig2 = 0;
    int n = (int)(1.0 / sqrt(eps));
    double y0 = y(a), yn = y(b);
    h = (b - a) / n;
        for (i = 1; i < n; i++) {
        if (i % 2 == 0)
            sig1 += y(a + i*h);
        else
            sig2 += y(a + i*h);
    }
        double cur_int, prev_int = h / 3 * (4 * sig2 + 2 * sig1 + y0 + yn), cur_even = sig2 + sig1;
        int cur_n = 2 * n;
        do{

            h = (b - a) / cur_n;
            double cur_odd = 0, xi;
            for (i = 1, xi = a + h; i < cur_n; i += 2, xi += 2 * h) {
                cur_odd += y(xi);
            }
            cur_int = h / 3 * (4 * cur_odd + 2 * cur_even + y0 + yn);
            prev_int = cur_int;
            cur_even = cur_odd + cur_even;
            cur_n *= 2;
        } while (!(((abs((cur_int - prev_int) / cur_int)) < eps) /*&& ((abs(cur_int - prev_int) / 15) < eps)*/));
    return cur_int;
}

double chebyshev(int n, double x) {
    double Tn1, Tn = x, Tn_1 = 1;

    if (n == 0) return Tn_1;

    int i = 1;
    while (i < n) {
        Tn1 = 2 * x * Tn - Tn_1;
        Tn_1 = Tn;
        Tn = Tn1;
        ++i;
    }

    return Tn;
}


double* gaussian_elimination(double **AM, int N) {
    int i, j, k, l;
    double temp1, temp2, temp3;
    double* result = new double[N];
    for (j = 0; j < N; j++) {
        temp1 = AM[j][j];
        for (i = j; i < N + 1; i++) {
            AM[j][i] /= temp1;
        }
        for (k = j+1; k < N; k++) {
            temp2 = AM[k][j];
            for (l = j; l < N + 1; l++) {
                AM[k][l] -= AM[j][l]*temp2;
            }
        }
    }
    for (j = N - 1; j >= 0; j--) {
        temp3 = 0;
        for (i = N - 1; i > j; i--) {
            temp3 += AM[j][i] * result[i];
        }
        result[j] = AM[j][N] - temp3;
    }
    return result;
}

double * makeA(double a, double b, double eps, int N, function<double(double)> y) {
    double ** A = new double*[N+1];
    int j, i;
    for (j = 0; j < N+1; j++) {
        A[j] = new double[N + 2];
    }
    for (j = 0; j < N+1; j++) {
        for (i = 0; i < N+1; i++) {
            if (j != i)
                A[j][i] = A[i][j] = simpson(a, b, eps, [&](double x) {return chebyshev(j, x) * chebyshev(i, x);});
            else
                A[j][i] = simpson(a, b, eps, [&](double x) {return chebyshev(j, x) * chebyshev(i, x);});
        }
        A[j][N+1] = simpson(a, b, eps, [&](double x) {return chebyshev(j, x) * y(x);});
    }
    auto res = gaussian_elimination(A, N+1);
    return res;
}

class LeastSqAppr {
    double * A;
    int N;
public:
    LeastSqAppr(function<double(double)> y, double a, double b, double eps, int N) {
        A = makeA(a, b, eps, N, y);
        this->N = N;
    }
    double operator()(double x) {
        double Tn1, Tn = x, Tn_1 = 1, sum;

        sum = Tn_1*A[0];
        if (N == 0) return sum;
        sum += Tn*A[1];
        int i = 1;
        while (i < N+1) {
            Tn1 = 2 * x * Tn - Tn_1;
            Tn_1 = Tn;
            Tn = Tn1;
            sum += Tn*A[i+1];
            ++i;
        }

        return sum;
    }
};

int main() {
    auto f = [&](double x) {return 7.5 * log10(x) * sin(x); };
    auto g = [&](double x) {return 1/x - 0.1 * x * x * sin(2*x); };

    double a = 2, b = 11;
    int N = 10;
    double result_n;
    double exact;
    int n = 9;
    double h = (b - a) / n, deviation_sum, delta;
    double x;
    double epsilon = 0.01;
    do {
        deviation_sum = 0;
        for (int i = 1; i < n; ++i) {
            x = a + h*i;
            exact = g(x);
            LeastSqAppr least(g, a, b, 0.00001, N);
            result_n = least(x);
            delta = abs(result_n - exact);
            deviation_sum = delta*delta + deviation_sum;
            printf("%.6f %.6f %.6f\n", x, least(x), g(x));
            //cout << "N = " << N << " exact " << exact << " apprx " << result_n  << endl;
        }
        cout << "N = " << N << " dev = " << (sqrt(deviation_sum / (n + 1))) << endl;
        ++N;
    } while (sqrt(deviation_sum / (n + 1)) > epsilon);
    /*
    for (int ord = 10; ord < 30; ord++) {
        LeastSqAppr least(g, 1, 10, 0.0000001, 8);
        if (sqrt((simpson(1, 10, 0.0001, [&](double x) {return (7.5 * log10(x) * sin(x) - least(x)) * (7.5 * log10(x) * sin(x) - least(x)); })) / 9) < 0.1){
            for (double x = 1; x <= 10; x += 0.5) {
                printf("%.6f %.6f %.6f\n", x, least(x), g(x));
            }
            printf("%f\n", sqrt((simpson(1, 10, 0.0001, [&](double x) {return (7.5 * log10(x) * sin(x) - least(x)) * (7.5 * log10(x) * sin(x) - least(x));})) / 9));
            break;
        }
    }*/
    return 0;
    //printf("%f\n", chebyshev(3, 5));
    //printf("%f", simpson(1, 10, 0.001, [&](double x) {return chebyshev(0, x) * chebyshev(1, x); }));
}
