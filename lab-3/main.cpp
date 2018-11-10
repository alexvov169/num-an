#include <iostream>
#include <iomanip>
#include "matrix.h"

using namespace std;

template<int n>
ostream& operator<<(ostream& stream, const double (&row)[n]) {
    stream << '(';
    for (const auto& i : row) {
        stream << std::setw(2) << i;
    }
    stream << ")\n";
    return stream;
}

int main() {
    double gsole[4][5] = {
        {14,	19,     15,     4,  	180},
        {17,	33,     5,      10, 	208},
        {11,	6,      28,     10, 	230},
        {6,     19,     3,      13, 	149}
    };
    // I've got no other idea than using found roots
    // to modify the matrix to convergable state
    // Even if I gain maximal diagonal item
    // it should be greater than sum of abs other items on row
    // as follows
    /*double ssole[4][5] = {
        //{17-6*2 + 5, 33-19*2 + 5, 5-3*2 + 5, 10-13*2 + 5, 208-149*2 + 5},
        {11,    1,      5,      -10,    208-149*2 + 6},
        {14,	19,     15,     4,  	180},
        {11,	6,      28,     10, 	230},
        {-15,   33,     -43,    46,     4*(149-180)+208}
        //{4*(6-14)+17, 4*(19-19)+33, 4*(3-15)+5, 4*(13-4)+10, 4*(149-180)+208}
    };*/
    // It must be something obvious but I can't see other way
    double ssole[4][5] = {
        {39,	19,     15,     4,  	230},
        {17,	33,     5,      10, 	208},
        {11,	6,      28,     10, 	230},
        {6,     19,     3,      29, 	229}
    };
    sole<double, 4> g(gsole), s(ssole);
    cout << "Gaussian elimination:" << endl;
    cout << "A = \n" << g << endl;
    try {
        cout << "x = " << g.gaussian_elimination() << endl;
    } catch (zero_leader_exception_impl<int> e) {
        cerr << e.what() << endl;
    }

    cout << "Seidel iteration:" << endl;
    cout << "A = \n" << s << endl;
    cout << "x = " << s.seidel_iteration(.00001) << endl;
    cout << "k = " << s.get_iterations() << endl;
    return 0;
}
