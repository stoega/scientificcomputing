#include <iostream>
#include <complex>
#include <chrono>

#include <SCvector.h>
#include <SCmatrix.h>

using namespace std;
using namespace SC;

// a template function, takes Matrix<double>, Matrix<complex<double> > etc
// uses operator(), which is defined for Matrix<T> but not LinearOperator<T>
template <typename T>
void SetToIdentity(Matrix<T>& A)
{
    for (int i=0; i<A.Height(); i++)
        for (int j=0; j<A.Width(); j++)
            A(i,j) = (i==j) ? 1 : 0.;
}

// template <typename T>
// void SetToIdentityOP(LinearOperator<T>& A)
// {
//     for (int i=0; i<A.Height(); i++)
//         for (int j=0; j<A.Width(); j++)
//             A(i,j) = (i==j) ? 1 : 0.;
// }

template <class MAT>
void SetToIdentityTemplate(MAT& A)
{
    for (int i=0; i<A.Height(); i++)
        for (int j=0; j<A.Width(); j++)
            A(i,j) = (i==j) ? 1 : 0.;
}

template<typename T>
void SetToIdentityDynamicCast(LinearOperator<T>& A)
{
    // dynamic cast - at runtime, the LinearOperator is cast to Matrix. 
    // if it is no Matrix, an exception is thrown
    Matrix<T> &Amat = dynamic_cast<Matrix<T>&>(A);

    // now, Amat is used instead of A
    for (int i=0; i<Amat.Height(); i++)
        for (int j=0; j<Amat.Width(); j++)
            Amat(i,j) = (i==j) ? 1 : 0.;
}


template<typename T>
void SetToIdentityDynamicCastPointer(LinearOperator<T>& A)
{
    // dynamic cast - at runtime, the LinearOperator is cast to Matrix*. 
    // if it is no Matrix, the null-pointer is returned
    Matrix<T> *Aptr = dynamic_cast<Matrix<T>*>(&A);

    if (!Aptr) 
    {
        cout << "A is not a Matrix, return" << endl;
        return;
    }

    // now, Aptr is used instead of A
    for (int i=0; i<Aptr->Height(); i++)
        for (int j=0; j<Aptr->Width(); j++)
            (*Aptr)(i,j) = (i==j) ? 1 : 0.;
}

int main()
{
    Matrix<> A(3);
    A.SetRandom();

    // A is a Matrix<>. Define a routine that sets A to the unit matrix
    SetToIdentity(A);

    // at some point, A is passed as a LinearOperator
    LinearOperator<>& Aop(A);
    // Aop is still a Matrix<T>, but the compiler doesn't know. virtual inheritance works for Aop, e.g. virtual void Print

    Aop.Print(cout);

    // however, these do not work
    // SetToIdentity(Aop);
    // SetToIdentityOP(Aop);
    // SetToIdentityTemplate(Aop);

    // we can try a dynamic cast
    SetToIdentityDynamicCast(Aop);
    SetToIdentityDynamicCastPointer(Aop);








}