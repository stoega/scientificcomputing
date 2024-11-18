
#define _USE_MATH_DEFINES
#include <math.h>
#include <complex>
#include <string>
#include <sstream>
#include <iostream>
#include <fstream>

#include <SCvector.h>
#include <SCmatrix.h>

#include <SCQRSolver.h>

#include <SCHessenbergmatrix.h>
#include <SCQRSolver_Hessenberg.h>
#include <SCArnoldi.h>

using namespace SC;

using namespace std;

void QR(const Matrix<double>& K, Vector<double>& lambda, Vector<double>& error)
{
    int n = K.Width();
    Matrix<double> A(K);
    QRSolver<double> QR(A);
    
    lambda.SetSize(A.Width());
    error.SetSize(A.Width()-1);

    int index = A.Width()-1;
    double mu;

    int step = 0;
    for (step=0; step<10000; step++)
    {
        // factorize, compute RQ
        QR.ComputeFactorization();
        QR.ComputeRQ(A);

        
        for (int i=0; i<A.Width()-1; i++)
            if (abs(A(i+1,i)) < 1e-10*(abs(A(i,i)) + abs(A(i+1,i+1)))) break;

        if (step%1000==0) cout << "step " << step << endl;
    }

    for (int i=0; i<A.Width(); i++)
        lambda(i) = A(i,i);
    for (int i=0; i<A.Width()-1; i++)
        error(i) = A(i+1,i);

}

void ArnoldiQR(const Matrix<double>& K, Vector<double>& lambda, Vector<double>& error)
{
    int n = K.Width();
    HessenbergMatrix<double> H(n);
    Matrix<double> V;
    Vector<double> res(n);
    res.SetRandom(-1,1);

    Arnoldi(K, res, n, V, H,1e-17);
    cout << "H " << H.Width() << endl;
    H.Set(n,n-1,0);

    // shifted QR
    QRSolverHessenberg<double> QRH(H);
    
    lambda.SetSize(H.Width());
    error.SetSize(H.Width()-1);

    // shifted QR iteration on
    double mu = 0;
    int index = H.Width()-1;

    int step = 0;
    for (step = 0; step<4*n; step++)
    {
        if (abs(H.Get(index,index-1)) < 1e-12*abs(H.Get(index,index)))
        {
            index--;
            if (index<=0) break;
            // cout << "step " << i << ": reduce matrix size .. " << index+1 << endl;
            H.SetSize(index+1);
        }

        // choose shift parameter: use diagonal matrix entry
        mu = H.Get(index,index);

        // shift matrix
        for (int i=0; i<H.Width(); i++)
            H.AddTo(i,i,-mu);
        // factorize, compute RQ
        QRH.ComputeFactorization();
        QRH.ComputeRQ(H);
        // shift back
        for (int i=0; i<H.Width(); i++)
            H.AddTo(i,i,mu);

        // if (i%1==0)
        // {
        //     cout << "EV(H) at step " << i << ": " << endl;
        //     for (int j=0; j<n; j++) cout << H.Get(j,j) << " ";//H.Get(j,j) << " ";
        //     cout << endl << "Error indicators ";
        //     for (int j=0; j<H.Width()-1; j++) cout << H.Get(j+1,j) << " "; //H.Get(j+1,j) << " ";
        //     cout << endl << endl;
        // }
        
    }

    H.SetSize(lambda.Size());
    for (int i=0; i<H.Width(); i++)
        lambda(i) = H.Get(i,i);
    for (int i=0; i<H.Width()-1; i++)
        error(i) = H.Get(i+1,i);


}

void PrintEigenfrequenciesToFile(const Vector<double> &lambda, const string& filename)
{
    ofstream outfile(filename.c_str());
    for (int i=0; i<lambda.Size(); i++)
        outfile << lambda(i) << "\t" << sqrt(lambda(i)) << endl;
    outfile.close();
}

int main()
{
    int n_x = 10;
    int n = n_x*n_x;

    double meshsize = 1./(n_x+1);
    double c2_coeff = 1.;

    int n_ex = 8;
    Vector<double> lambda_exact(n_ex*(n_ex+1)/2);
    int ii=0;
    for (int i=1; i<=n_ex; i++)
        for (int j=1; j<=i; j++)
            lambda_exact(ii++) = c2_coeff*M_PI*M_PI*(i*i+j*j);

    double meshsize2 = meshsize*meshsize;

    // stiffness matrix K
    Matrix<> K(n,n);
    K.SetAll(0.);

    for (int i=0; i<n_x; i++)
    {
        for (int j=0; j<n_x; j++)
        {
            K(i*n_x+j, i*n_x+j) = 4.*c2_coeff/meshsize2;
            if (i>0)     K(i*n_x+j, (i-1)*n_x+j) = -c2_coeff/meshsize2;
            if (i<n_x-1) K(i*n_x+j, (i+1)*n_x+j) = -c2_coeff/meshsize2;
            if (j>0)     K(i*n_x+j, (i)*n_x+j-1) = -c2_coeff/meshsize2;
            if (j<n_x-1) K(i*n_x+j, (i)*n_x+j+1) = -c2_coeff/meshsize2;
        }
    }

    // M = I is not generated explicitely

    Vector<double> lambda;
    Vector<double> error;

    // QR(K, lambda, error);
    ArnoldiQR(K, lambda, error);
    
    cout << "EV: " << endl << lambda << endl;

    cout << "---------------" << endl;

    cout << "EV analytic " << endl << lambda_exact << endl;

    stringstream filename, filename_exact;
    filename << "eigenfrequencies_n" << n_x << ".txt"; 
    filename_exact << "eigenfrequencies_exact.txt";
    PrintEigenfrequenciesToFile(lambda, filename.str());
    PrintEigenfrequenciesToFile(lambda_exact, filename_exact.str());
    // cout << endl << "Error indicators " << endl << error <<  endl;
    cout << "Total error: " << error.Norm() << endl;

}