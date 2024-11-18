#include <iostream>
#include <complex>
#include <chrono>

#include <SCvector.h>
#include <SCmatrix.h>

using namespace std;
using namespace SC;


typedef complex<double> COMPLEX;    // write COMPLEX instead of complex<double>



int main()
{

    // call constructor/destructor in debug mode
    {
        Matrix<double> A;
    }

    // a simple example
    int size = 3;
    Vector<COMPLEX> a, b, c;
    a.SetSize(size); a.SetAll(1.);
    cout << "a = " << a << endl;

    Matrix<COMPLEX> M(4,size);
    M.SetAll(2.);
    M(3,2) = COMPLEX(0.,1.);
    M(0,2) = COMPLEX(0.,1.);

    cout << "M = " << M << endl;

    M.Apply(a, b);
    M.ApplyH(b, a);
    
    cout << "b = M.a\nb = " << b << endl;

    Matrix<COMPLEX> N(M);
    N.Transpose();
    cout << "\ntransposed matrix " << N << endl;

    return 0;
}