#include <SCvector.h>
#include <SCmatrix.h>

#include <SCGaussSolver.h>

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
        if (i > 0) A(i,i-1) = -0.5;
    }


    SimpleGaussSolver<> invA_simple(A);

    Vector <> c, d;
    invA_simple.Apply(b, c);

    cout << "invA_simple * b = " << c << endl;

    A.Apply(c,d);

    cout << "A * invA_simple * b = " << d << endl;
    

    GaussSolver<> invA(A);

    Vector <> c2, d2;
    invA.Apply(b, c2);

    cout << "invA * b = " << c2 << endl;

    A.Apply(c2,d2);

    cout << "A * invA * b = " << d2 << endl;

    cout << "inverse matrix inva = " << invA.GetInverse() << endl;


}