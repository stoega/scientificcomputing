#include <complex>
#include <SCvector.h>
#include <SCmatrix.h>

#include <Eigen/Dense>

#ifdef USE_MKL
#define MKL_Complex16 std::complex<double>
#include <mkl.h>
#else // USE_MKL
extern "C" void dgesvd_(char *jobu, char *jobvt, int *m, int *n,
                double *a, int *lda, double *s, double *u,
                int *ldu, double *vt, int *ldvt, double *work,
                int *lwork, int *info);

extern "C" void  zgesv_(int* N, int* NRHS, std::complex<double>* A, int* LDA, int* IPIV, std::complex<double>* B, int* LDB, int* INFO);
#endif // USE_MKL
// #include <SCSVD.h>

using namespace SC;
using namespace std;

enum { DS = Eigen::Dynamic }; // Eigen-lib dynamic size identifier

void SVDExample()
{
    int m = 4, n = 5;

    Matrix<> A(m,n);
    for (int i=0; i<m; i++)
        for (int j=0; j<n; j++)
            A(i,j) = i+j;

    Eigen::Matrix<double, DS, DS> EigenA(A.Height(),A.Width());
    EigenA.setConstant(0.);
    for (int i=0; i<A.Height(); i++)
        for (int j=0; j<A.Width(); j++)
            EigenA(i,j) = A(i,j);

    Eigen::JacobiSVD<Eigen::Matrix<double, DS, DS>, Eigen::ComputeFullU | Eigen::ComputeFullV> SVD_A(EigenA);

    cout << "Using EIGEN" << endl;
    cout << "singular values: " <<  SVD_A.singularValues() << endl;
    cout << "matrix u" << SVD_A.matrixU() << endl;
    cout << "matrix v" << SVD_A.matrixV() << endl;

    Matrix<> copyA(m,n,TMatrixStorage::COL_MAJOR);
    for (int i=0; i<m; i++)
        for (int j=0; j<n; j++)
            copyA(i,j) = A(i,j);

    Vector<> S(min(m,n));
    Matrix<> U(m, TMatrixStorage::COL_MAJOR), VT(n,TMatrixStorage::COL_MAJOR);
    Vector<double> work(m*n+100);
    int info;
    char jobu = 'A', jobv = 'A';
    int lda = m, ldu = m, ldv = n;
    int lwork = work.Size();

    dgesvd_ ( &jobu, &jobv, &m, &n, &copyA(0,0), &lda,
              &S(0),
              &U(0,0), &ldu, &VT(0,0), &ldv,
              &work(0), &lwork, 
              &info);

    cout << "Using LAPACK" << endl;
    cout << "U " << U << endl;
    cout << "V^T " << VT << endl;
    cout << "sing " << S << endl;
}

void SolveExample()
{
    int n = 5;
    Matrix<complex<double>> B(n, TMatrixStorage::COL_MAJOR);
    B.SetRandom();

    // copy matrix B (is destroyed, contains LU factorization on exit
    Matrix<complex<double>> copyB(B);
    cout << "Matrix B " << copyB << endl;

    Vector<complex<double>> f(n);
    f.SetRandom(0,20,1000);

    // copy right hand side into solution vector, is passed to LAPACK
    Vector<complex<double>> x(f);
    cout << "Right hand side f " << x << endl;

    int NRHS = 1;
    Vector<int> ipiv(n);
    int info;

    zgesv_ 	( &n, &NRHS, &copyB(0,0), &n,
		&ipiv(0),
		&x(0),
		&n, &info );
    if (info>0)
        cout << "Matrix singular" << endl;
    if (info<0)
        cout << "Input " << -info << " illegal" << endl;
        
    cout << "Solution x " << x << endl;

    Vector<complex<double>> f2;
    B.Apply(x,f2);
    f2.AddMultiple(-1, f);
    cout << "error " << f2.Norm() << endl;

}

int main()
{



    SVDExample();


    SolveExample();


}
