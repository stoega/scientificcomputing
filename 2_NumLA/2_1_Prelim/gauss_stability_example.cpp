#include <random>
#include <iostream>
#include <iomanip>

#include <SCvector.h>
#include <SCmatrix.h>

using namespace SC;
using namespace std;

template <typename T>
void PrintSystem(const Matrix<T>& mat, const Vector<T>& b)
{
    // cout << std::scientific;
    cout << std::setprecision(16);
    for (int i=0; i<mat.Height(); i++)
    {
        for (int j=0; j<mat.Width(); j++)
            cout << mat(i,j) << " ";
        cout << " | " << b(i) << endl;
    }
    cout << endl;
}


template <typename T>
void GaussSolve(const Matrix<T>& mat, const Vector<T>& b, Vector<T>& x)
{
    // copy matrix a, this copy is destroyed/turned to identity
    Matrix<T> temp(mat);

    x = b;

    T fac, pivinv;

    PrintSystem(temp, x);
    for (int i = 0; i < mat.Height(); i++)
    {
        if (temp(i, i) == 0.0)
            throw("GaussSolve: zero pivot element");

        pivinv = 1.0 / temp(i, i);
        for (int l = 0; l < mat.Height(); l++)
        {
            temp(i, l) *= pivinv;
        }
        x(i) *= pivinv;

        for (int ll = 0; ll < mat.Height(); ll++)
            if (ll != i)
            {
                fac = temp(ll, i);
                for (int l = 0; l < mat.Height(); l++)
                    temp(ll, l) -= temp(i, l) * fac;
                x(ll) -= fac * x(i);
            }
        PrintSystem(temp,x);
    }
}

int main()
{


    Vector<> b(3), x(3), x_exact(3);
    x_exact.SetAll(1.);
    Matrix<> A(3);

    int l0 = 0;
    int l1 = 1;
    A(l0,0) = 1.; 
    // A(l0,0) = 1e-12; 
    A(l0,1) = 8.; A(l0,2) = 1.;
    A(l1,0) = 3.; A(l1,1) = 7.; A(l1,2) = 5.;
    A(2,0) = 8.; A(2,1) = 5.; A(2,2) = 1.;
    
    A.Apply(x_exact, b);
    cout << "System matrix A = " << A << endl;
    cout << "right hand side b = " << b << endl << endl;
    
    GaussSolve(A, b, x);
    x -= x_exact;
    cout << "error vector = " << x << endl;
    cout << "error norm = " << x.Norm() << endl;



}