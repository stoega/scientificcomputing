#include <iostream>

#include <SCmatrix.h>
#include <SCvector.h>
#include <SCLUSolver.h>
#include <SCGaussSolver.h>

using namespace SC;
using namespace std;

template <typename T>
void SetTridiag(Matrix<T> &mat, T diag, T offdiag)
{
    if (!(mat.Height() == mat.Width())) throw "SetTridiag for non-square matrix";

    for (int i=0; i<mat.Height(); i++)
        mat(i,i) = diag;
    for (int i=0; i<mat.Height()-1; i++)
    {
        mat(i,i+1) = offdiag;
        mat(i+1,i) = offdiag;
    }
}


int main()
{
    Matrix<double> a(10);
    SetTridiag(a, 2., -1.);

    cout << "Tridiagonal system matrix "  << endl << a << endl;

    cout << "--- computing the inverse --- \ninverse matrix is dense: " << endl;

    GaussSolver<double> ainv(a);
    cout << "ainv = " << ainv.GetInverse() << endl;

    cin.get();

    cout << "--- computing LU factorization --- \nLU factors are tridiag: " << endl;
    LUSolver<double> a_LU(a);
    cout << a_LU << endl;
    cout << "LU factorization " << endl;
    a_LU.PrintPattern(cout);
    cout << endl;

    cin.get();

    cout << "--- setting additional entries in matrix a ---" << endl;
    a(5,1) = -0.044;
    a(0,7) = -0.011;
    a(8,3) = -0.066;

    // for (int i=1; i<10; i++) a(i,9) = -0.2;
    LUSolver<double> a_LU2(a);



    cout << "a = " << a << endl;
    cin.get();
    cout << "fill-in in LU factorization " <<  endl;
    a_LU2.PrintPattern(cout);
    cout << endl;


    return 0;
}