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
    // LinearOperator<> B;
    // call constructor/destructor in debug mode
    {
        Matrix<double> A(3,2);
        Matrix<double> B;
        B = A;
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
    cout << "b = M.a\nb = " << b << endl;
    M.ApplyH(b, a);
    
    cout << "a = M^H.b\na = " << a << endl;

    // test virtual/non-virtual inheritance -- what happens, if virtual/override keywords are removed in matrix class?
    {
        // operator<< declaration:
        // ostream& operator<<(ostream& os, const LinearOperator<T>& m);
        // M is passed as LinearOperator<T>
        cout << M << endl;
        // M is a Matrix<COMPLEX>
        // at compile-time, refM is treated as LinearOperator<COMPLEX>
        // at run-time, refM will be a Matrix<COMPLEX>, but the compiler does not know
        LinearOperator<COMPLEX>& refM(M);
        refM.Print(cout);
    }

    // transpose a matrix
    Matrix<COMPLEX> N(M);
    N.Transpose();
    cout << "\ntransposed matrix " << N << endl;


    return 0;
}