#include <iostream>
#include <tuple>
#define _USE_MATH_DEFINES
#include <cmath>
#include "ranges.h"
#include "isolator.h"
#include <stdlib.h>
#include <unistd.h>
using namespace std;

typedef double value_type;
typedef int integer_type;
typedef edges<value_type> edges_type;

auto f = [](value_type x) { return 15 * sqrt(1 + cos(x)) + 1.5 * x; };
auto df = [](value_type x) { return 1.5 - 15 * sin(x) / (2 * sqrt(1 + cos(x))); };
auto ddf = [](value_type x) { return - 15 * sin(x) / (2 * sqrt(1 + cos(x))); };
auto phi = [](value_type x) { return 10 * sqrt(1 + cos(x)); };

double increment(double a, double b, double eps, double &delta, int &iters){
    double m1 = abs(df(a)), M1 = abs(df(b));
    if (m1 > M1){
        swap(m1, M1);
    }
    cout << "m1 = " << m1 << " M1 = " << M1 << endl;
    double x0, xk = (a+b) / 2;
    double lambda = 1/M1, q = 1 - m1/M1;
    iters = 0;
    do {
        x0 = xk;
        xk = x0 - lambda * f(x0);
        iters++;
    } while (abs(x0-xk) > eps*(1-q)/q);

    delta = abs(x0-xk)*q/(1-q);
    return xk;
}

void approximate_root(edges_type edges) {
    double delta;
    int iters;
    cout << '[' << edges.left() << ", " << edges.right() << ']' << endl;
    edges.shrink([](const edges_type& e) {
        usleep(100000);
        return e.contain_root(f);
    });
    cout << "result = " << edges << ' '
         << '[' << f(edges.left()) << ", " << f(edges.right()) << ']' << endl;
    cout << increment(edges.left(), edges.right(), .01, delta, iters) << endl;
    cout << "delta = " << delta << " iters = " << iters << endl;
}

void print_root_table() {

}

int main()
{
    isolate_roots<int, double>(f, approximate_root);
    //cout << "Hello World!" << endl;
    edges_type e = {-10, 0};
    cout << e << endl;
    e.shrink([](const edges_type& e) {
        return e.contain_root([](value_type x) {
            usleep(100000);
            return -x*x + 16;
        });
    });
    cout << e << endl;
    e = {0, 10};
        cout << e << endl;
        e.shrink([](const edges_type& e) {
            return e.contain_root([](value_type x) {
                usleep(100000);
                return -x*x + 16;
            });
        });
        cout << e << endl;
    int i;
    double d;
    cout << increment(-11.5, -9.3, .01, d, i);
    return 0;
}
