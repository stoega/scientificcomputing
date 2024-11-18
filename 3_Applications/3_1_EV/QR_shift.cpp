#include <iostream>
#include <chrono>

#define _USE_MATH_DEFINES
#include <math.h>
#include <complex>

#include <SCvector.h>
#include <SCmatrix.h>

#include <SCQRSolver.h>

#include <SCHessenbergmatrix.h>
#include <SCQRSolver_Hessenberg.h>
#include <SCArnoldi.h>

#include <Eigen/Eigenvalues>
typedef Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic> EigenMatrix;

using namespace SC;

using namespace std;

int main()
{
    int size = 15;
    Matrix<complex<double>> A(size);
    A.SetRandom(-1,1);

    HessenbergMatrix<complex<double>> H(size);
    Matrix<complex<double>> V;
    Vector<complex<double>> res(size);
    res.SetRandom(0,1);

    auto begin = std::chrono::high_resolution_clock::now();

    // Arnoldi iteration to compute H, V
    Arnoldi(A, res, size, V, H);

    QRSolverHessenberg<complex<double>> QRH(H);
    
    // shifted QR iteration on
    complex<double> mu;
    int index = size-1;
    int step = 0;
    for (step=0; step<1000; step++)
    {
        if (abs(H.Get(index,index-1)) < 1e-10*abs(H.Get(index,index)) || abs(H.Get(index,index-1)) < 1e-12)
        {
            index--;
            if (index<=0) break;
            cout << "step " << step << ": consider submatrix of size .. " << index+1 << endl;
            H.SetSize(index+1);
        }

        // choose shift parameter: use diagonal matrix entry
        mu =  H.Get(index,index);

        // shift matrix
        for (int i=0; i<H.Width(); i++)
            H.AddTo(i,i,-mu);
        // factorize, compute RQ
        QRH.ComputeFactorization();
        QRH.ComputeRQ(H);
        // shift back
        for (int i=0; i<H.Width(); i++)
            H.AddTo(i,i,mu);

        // if (step%1==0)
        // {
        //     cout << "EV(H) at step " << step << ": " << endl;
        //     for (int j=0; j<size; j++) cout << H.Get(j,j) << " ";
        //     cout << endl << "Absolute values: ";
        //     for (int j=0; j<H.Width(); j++) cout << abs(H.Get(j,j)) << " "; 
        //     cout << endl << "Error indicators ";
        //     for (int j=0; j<H.Width()-1; j++) cout << H.Get(j+1,j) << " ";
        //     cout << endl << endl;
        // }
    }
    auto end = std::chrono::high_resolution_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);

    H.SetSize(size);
    cout << "Eigenvalues (H): " << endl;
    for (int j=0; j<H.Width(); j++) cout << H.Get(j,j) << " ";
    cout << endl << "Absolute values: ";
    for (int j=0; j<H.Width(); j++) cout << abs(H.Get(j,j)) << " ";
    cout << endl << "Error indicators ";
    for (int j=0; j<H.Width()-1; j++) cout << H.Get(j+1,j) << " ";
    cout << endl << endl;
    cout << endl << "Timing: " << elapsed.count()/1e9 << endl;
    cout << "Converged after " << step << " steps" << endl;

    // use Eigen lib to compute eigenvalues
    EigenMatrix A_eigen(size,size);
    for (int i=0; i<size; i++)
        for (int j=0; j<size; j++)
            A_eigen(i,j) = A(i,j);
    
    Eigen::ComplexEigenSolver<EigenMatrix> ces;
    ces.compute(A_eigen);

    cout << "-----------------------\n";
    cout << "Eigenvalues computed using Eigen lib" << endl;
    for (int i=0; i<size; i++) cout << ces.eigenvalues()[size-1-i] << " ";
    cout << endl;
    cout << "absolute values: ";
    for (int i=0; i<size; i++) cout << abs(ces.eigenvalues()[size-1-i]) << " ";
    cout << endl;

}