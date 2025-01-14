#include <SCvector.h>
#include <SCmatrix.h>

#include <SCLUSolver.h>
#include <SCQRSolver.h>

using namespace SC;
using namespace std;

int main()
{
    int size = 8;

    Vector<> b(size);
    b.SetAll(3.);

    Matrix<> A(size,size);
    A.SetAll(0.);
    for (int i=0; i<size; i++)
    {
        A(i,i) = 2.;
        if (i > 0) A(i-1,i) = -0.5;
        // symmetric A
        if (i > 0) A(i,i-1) = -0.5;
    }
    // A.SetRandom();
    // A(0,0) = 1e-5;

    LUSolver<> invA_LU(A);

    Vector <> c, d;
    invA_LU.Apply(b, c);

    cout << "LU Solver:"  << endl;
    cout << "invA * b = " << c << endl;

    A.Apply(c,d);

    cout << "A * invA * b = " << d << endl;

    d -= b;
    cout << "relative error residual, LU " << d.Norm()/b.Norm() << endl << endl;

    cout << "factorization: " << invA_LU << endl;

    QRSolver<> invA_QR(A);
    Vector <> cQR, dQR;
    invA_QR.Apply(b, cQR);

    cout << "QR Solver:" << endl;
    cout << "invA * b = " << cQR << endl;

    A.Apply(cQR,dQR);

    cout << "A * invA * b = " << dQR << endl;
    // cout << "QR factorization " << invA_QR << endl;
    dQR -= b;
    cout << "relative error residual, QR " << dQR.Norm()/b.Norm() << endl;

    cout << "factorization: " << invA_QR << endl;


}