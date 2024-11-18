#include <SCvector.h>
#include <SCmatrix.h>

#include <SCCGSolver.h>

using namespace SC;
using namespace std;

int main()
{
    int size = 10;

    Vector<> b(size);
    b.SetAll(3.);

    Matrix<> A(size,size);
    A.SetAll(0.);
    for (int i=0; i<size; i++)
    {
        // A = I -> 1 step
        A(i,i) = 1.;
        // // A = diag(i+1) -> Eigenvalues lam_i = i+1, condition n, n steps
        // A(i,i) = i+1;        
        // // A = diag((i+1)^2) -> Eigenvalues lam_i = (i+1)^2, condition n^2, slow convergence
        // A(i,i) = (i+1)*(i+1);
        // // A = diag(i%k+1) -> k different eigenvalues 1..k, k steps needed for exact solution
        // A(i,i) = i%3+1;
        // // "1D finite difference" tridiag matrix
        // A(i,i) = 2.;
        // if (i > 0) A(i-1,i) = -1.;
        // if (i > 0) A(i,i-1) = -1.;
    }

    CGSolver<> invA_CG(A);

    Vector <> c, d;
    invA_CG.Apply(b, c);

    if (size < 20)
    {
        cout << "invA * b = " << c << endl;

        A.Apply(c,d);

        cout << "A * invA * b = " << d << endl;
    }

}