#include <iostream>
#include "matrix.h"

using namespace std;

int main()
{
    int mm[3][4] = {
        {1, 2, 3, 4},
        {3, 4, 3, 2},
        {2, 3, 8, 3}
    };
    sole<int, 3> s(mm);
    cout << s << endl;
    cout << "Hello World!" << endl;
    return 0;
}
