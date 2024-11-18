#include <SCvector.h>
#include <SCmatrix.h>

#include <SCLUSolver.h>
#include <SCQRSolver.h>

using namespace SC;
using namespace std;

typedef complex<double> T;

int main()
{
    int size = 8;

    Vector<T> b(size);
    b.SetAll(T(3.,2.));

    Matrix<T> A(size,size);
    A.SetAll(0.);
    for (int i=0; i<size; i++)
    {
        A(i,i) = T(2.,1.);
        if (i > 0) A(i-1,i) = T(-0.5,-1);
        // hermitian matrix
        // if (i > 0) A(i,i-1) = T(-0.5,1);
    }
    // A.SetRandom();
    LUSolver<T> invA_LU(A);

    Vector <T> c, d;
    invA_LU.Apply(b, c);

    cout << "LU Solver:"  << endl;
    cout << "invA * b = " << c << endl;

    A.Apply(c,d);

    cout << "A * invA * b = " << d << endl;


    QRSolver<T> invA_QR(A);
    Vector <T> cQR, dQR;
    invA_QR.Apply(b, cQR);

    cout << "QR Solver:" << endl;
    cout << "invA * b = " << cQR << endl;

    A.Apply(cQR,dQR);

    cout << "A * invA * b = " << dQR << endl;
    cout << "QR factorization " << invA_QR << endl;


}