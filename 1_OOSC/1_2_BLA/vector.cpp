#include <iostream>
#include <complex>

// include using absolute path or relativ path wrt. current folder: use "..."
// #include "../../SClinalg/SCvector.h"

// include from include directories (target_include_directories in CMakeLists.txt): use <...>
#include <SCvector.h>

using namespace std;
using namespace SC;


int main()
{

    int size = 5 ;
    // allocate vectors of size size
    Vector<> a(size);
    Vector<complex<double>> z(size);


    // set all entries of a to one, of b to i
    a.SetAll(1.);
    z.SetAll(complex<double>(0,1));

    // copy constructor - allocate new vector, set c = a
    Vector <> c(a);
    // assignment operator
    Vector <complex<double>> y;
    y = z;

    // move assignment
    // y = Vector<complex<double>>(4);

    // private member access
    int len_a = a.Size();
    a.SetSize(4);
    a.SetAll(0);
    a(0) = 1;
    a(2) = 2;
    a(1) = 3;

    // range check --------------------
    // in Debug mode: exception
    // in Release mode: undefined behavior
    a(5) = 7;
    cout << "non-existing element a(5) = " << a(5) << endl;

    // compute norm of vector a
    cout << "-----------------------------" << endl;
    double norm = a.Norm();
    cout << "norm of " << a << " is " << norm << endl;

    // ----------------------------
    // template specialization
    cout << "*** set random vectors (template specialization) ***" << endl;
    a.SetRandom(-1,1);
    z.SetRandom(0,100);
    Vector<size_t> v(3);
    v.SetRandom();

    // output of vector
    cout << "a = " << a << endl;
    cout << "z = " << z << endl;

    // -----------------------------------
    // Vector<double> d(2);
    // d(0) = 3; d(1) = 4;


    // cout << "*** basic linear algebra operations ***" << endl;
    // cout << "c = " << c << endl;
    // cout << "y = " << y << endl;
    // c *= 3.;
    // y += z;
    // cout << "3*c = " << c << endl;
    // cout << "y+z = " << y << endl;

    // // does not work, vectors are not of same type - compiler error
    // // y += a;

    // cout << "*** increase size of vector repeatedly ***" << endl;
    // a.AllocateMemory(5);
    // a.SetSize(0);
    // for (int i=1; i<=5; i++)
    // {
    //     a.SetSize(i);
    //     a(i-1) = i;
    //     cout << "iteration " << i << ": a = " << a << endl;
    // }

    // a.SetSize(3);
    // a(0) = 1;
    // a(2) = 2;
    // a(1) = 3;
    // cout << "-----------------------------" << endl;
    // cout << "summing vectors of different size:\nin Debug mode: exception\nin Release mode, there is no check" << endl;
    // cout << "Release mode: adding larger vector does work for this implementation " << endl;
    // cout << a << " += " << c << ":" << endl;
    // a += c;
    // cout << "a = " << a << endl;

    // cout << "-----------------------------" << endl;
    // cout << a << " += " << d << ":" << endl;
    // a += d;
    // cout << "Release mode: adding smaller vector, probably memory corruption (or wrong results in last entry, should be [7, 10, 5]" << endl;
    // cout << "a = " << a << endl;


    return 0;
}