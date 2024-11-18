#include <SCvector.h>
#include <SCmatrix.h>

#include <SCGMRESSolver.h>

using namespace SC;
using namespace std;

int main()
{
    int size = 10;

    Vector<complex<double>> b(size);
    b.SetAll(3.);

    Matrix<complex<double>> A(size,size);

    A.SetRandom(-1.,1.);
    // diagonal-dominant matrix -- better convergence
    for (int i=0; i<size; i++)
        A(i,i) = complex<double>(5.,5.);

    GMRESSolver<complex<double>> invA(A, 200, 9, 1e-8);

    Vector <complex<double>> c, d;
    invA.Apply(b, c);

    // cout << "invA * b = " << c << endl;

    // A.Apply(c,d);

    // cout << "A * invA * b = " << d << endl;

}