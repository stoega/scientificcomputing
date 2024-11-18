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

    // cout << "A " << A << endl;

    Vector<> *x1;
    x1 = new Vector<>(size);
    x1->SetAll(2.);
    cout << *x1 << endl;
    delete x1;
    
    Matrix<complex<double>> Q(size);
    Matrix<complex<double>> A1(size);
    A1 = A;

    int nsteps = 100000;
    // nsteps = 2;
    auto begin_1 = std::chrono::high_resolution_clock::now();

    // set up QR solver
    QRSolver<complex<double>> QR(A1);

    // QR iteration - QR factorization is performed in same memory space in each step
    // no "new" within loop
    int step=0;
    bool is_converged = false;
    while (step<nsteps && !is_converged)
    {
        step++;

        QR.ComputeFactorization();
        QR.ComputeRQ(A1);

        is_converged = true;
        for (int i=0; i<A1.Width()-1; i++)
            if (abs(A1(i+1,i)) > 1e-10*(abs(A1(i,i)) + abs(A1(i+1,i+1)))) 
            {
                is_converged = false;
                break;
            }

        // cout << "next iterate " << endl << A1 << endl;

        // if (step%100==0)
        // {
        //     cout << "EV at step " << step << ": " << endl;
        //     for (int j=0; j<size; j++) cout << A1(j,j) << " ";
        //     cout << endl << "Absolute values: ";
        //     for (int j=0; j<size; j++) cout << abs(A1(j,j)) << " ";
        //     cout << endl << "Error indicators ";
        //     for (int j=0; j<size-1; j++) cout << A1(j+1,j) << " ";
        //     cout << endl;
        // }

    }
    cout << "QR factorization needed " << step << " steps" << endl;

    auto end_1 = std::chrono::high_resolution_clock::now();
    auto elapsed_1 = std::chrono::duration_cast<std::chrono::nanoseconds>(end_1 - begin_1);


    HessenbergMatrix<complex<double>> H(size);
    Matrix<complex<double>> V;
    Vector<complex<double>> res(size);
    res.SetRandom(0,1);

    auto begin_2 = std::chrono::high_resolution_clock::now();

    // Arnoldi iteration to compute H, V
    Arnoldi(A, res, size, V, H);

    QRSolverHessenberg<complex<double>> QRH(H);
    
    // QR iteration on H

    step = 0;
    is_converged = false;
    while (step<nsteps && !is_converged)
    {
        step++;

        QRH.ComputeFactorization();
        QRH.ComputeRQ(H);

        is_converged = true;
        for (int i=0; i<H.Width()-1; i++)
            if (abs(H.Get(i+1,i)) > 1e-10*(abs(H.Get(i,i)) + abs(H.Get(i+1,i+1)))) 
            {
                is_converged = false;
                break;
            }

        // if (step%100==0)
        // {
        //     cout << "EV(H) at step " << step << ": " << endl;
        //     for (int j=0; j<size; j++) cout << H.Get(j,j) << " ";
        //     cout << endl << "Absolute values: ";
        //     for (int j=0; j<size; j++) cout << abs(H.Get(j,j)) << " ";
        //     cout << endl << "Error indicators ";
        //     for (int j=0; j<size-1; j++) cout << H.Get(j+1,j) << " ";
        // }
    }
    cout << "QR(Hessenberg) factorization needed " << step << " steps" << endl;

    auto end_2 = std::chrono::high_resolution_clock::now();
    auto elapsed_2 = std::chrono::duration_cast<std::chrono::nanoseconds>(end_2 - begin_2);


    // use Eigen lib to compute eigenvalues
    EigenMatrix A_eigen(size,size);
    for (int i=0; i<size; i++)
        for (int j=0; j<size; j++)
            A_eigen(i,j) = A(i,j);
    
    Eigen::ComplexEigenSolver<EigenMatrix> ces;
    ces.compute(A_eigen);





    cout << "-----------------------\n";
    cout << "Eigenvalues computed using QR on full matrix A " << endl;
    for (int j=0; j<size; j++) cout << A1(j,j) << " ";
    cout << endl << "Absolute values: ";
    for (int j=0; j<size; j++) cout << abs(A1(j,j)) << " ";
    cout << endl << "Error indicators ";
    for (int j=0; j<size-1; j++) cout << A1(j+1,j) << " ";
    cout << endl << "Timing: " << elapsed_1.count()/1e9 << endl;
    cout << "-----------------------\n";


    cout << "Eigenvalues computed using QR after Arnoldi reduction to Hessenberg matrix " << endl;
    for (int j=0; j<size; j++) cout << H.Get(j,j) << " ";
    cout << endl << "Absolute values: ";
    for (int j=0; j<size; j++) cout << abs(H.Get(j,j)) << " ";
    cout << endl << "Error indicators ";
    for (int j=0; j<size-1; j++) cout << H.Get(j+1,j) << " ";
    cout << endl << "Timing: " << elapsed_2.count()/1e9 << endl;


    cout << "-----------------------\n";
    cout << "Eigenvalues computed using Eigen lib" << endl;
    for (int i=0; i<size; i++) cout << ces.eigenvalues()[size-1-i] << " ";
    cout << endl;
    cout << "absolute values: ";
    for (int i=0; i<size; i++) cout << abs(ces.eigenvalues()[size-1-i]) << " ";
    cout << endl;


}