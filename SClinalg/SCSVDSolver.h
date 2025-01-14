#pragma once

#include <iostream>
#include <complex>

#include <SCvector.h>
#include <SCmatrix.h>

#include <Eigen/Dense>

#ifdef USE_MKL
#define __HAVE_LAPACK__
#include <mkl.h>
#endif // USE_MKL

#ifdef USE_LAPACK
#define __HAVE_LAPACK__
extern "C" void dgesvd_(char *jobu, char *jobvt, int *m, int *n,
                double *a, int *lda, double *s, double *u,
                int *ldu, double *vt, int *ldvt, double *work,
                int *lwork, int *info);
#endif // USE_LAPACK


namespace SC
{
    class SVDSolver : public LinearOperator<double>
    {
    protected:
        const Matrix<double> & a;
        Matrix<double> U; 
        Matrix<double> VT; 
        Vector<double> Sigma;
    public:

        SVDSolver() = delete;

        SVDSolver(const SVDSolver &) = delete;
        SVDSolver &operator=(const SVDSolver &) = delete;
        
        SVDSolver(const Matrix<double> & A) 
        : LinearOperator<double>(A.Width(), A.Height()), a(A)
        {
            U.SetSize(A.Height(), A.Height());
            VT.SetSize(A.Width(), A.Width());
            // lapack needs col-major storage
            U.SetStorage(TMatrixStorage::COL_MAJOR);
            VT.SetStorage(TMatrixStorage::COL_MAJOR);

            Sigma.SetSize(std::min(A.Height(), A.Width()));

            ComputeFactorization();
        }

        ~SVDSolver() 
        {
            ;
        }

        Matrix<>& GetU() { return U; }
        Matrix<>& GetVT() { return VT; }
        Vector<>& GetSigma() { return Sigma; }
        const Matrix<>& GetU() const { return U; }
        const Matrix<>& GetVT() const { return VT; }
        const Vector<>& GetSigma() const { return Sigma; }

        void ComputeFactorization()
        {
#ifdef __HAVE_LAPACK__
            Vector<double> work(a.Height()*a.Width()+100);
            Matrix<double> copya(a.Height(), a.Width(), TMatrixStorage::COL_MAJOR);
            for (int i=0; i<a.Height(); i++)
                for (int j=0; j<a.Width(); j++)
                    copya(i,j) = a(i,j);

            int info;
            char jobu = 'A', jobv = 'A';
            int m = a.Height(); int n = a.Width();
            int lda = a.Height(), ldu = a.Height(), ldv = a.Width();
            int lwork = work.Size();

            dgesvd_  ( &jobu, &jobv, &m, &n, &copya(0,0), &lda,
                    &Sigma(0),
                    &U(0,0), &ldu, &VT(0,0), &ldv,
                    &work(0), &lwork, 
                    &info);

            // std::cout << "U " << U << std::endl;
            // std::cout << "VT " << VT << std::endl;
            // std::cout << "sing " << Sigma << std::endl;
#else // __HAVE_LAPACK__
            Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> eigenA(a.Height(), a.Width());
            for (int i=0; i<a.Height(); i++)
                for (int j=0; j<a.Width(); j++)
                    eigenA(i,j) = a(i,j);

            Eigen::JacobiSVD<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>, Eigen::ComputeFullU | Eigen::ComputeFullV> SVD_A(eigenA);

            // std::cout << "singular values: " <<  SVD_A.singularValues() << std::endl;
            // std::cout << "matrix u" << SVD_A.matrixU() << std::endl;
            // std::cout << "matrix v" << SVD_A.matrixV() << std::endl;

            for (int i=0; i<a.Height(); i++)
                for (int j=0; j<a.Height(); j++)
                    U(i,j) = SVD_A.matrixU()(i,j);
            for (int i=0; i<a.Width(); i++)
                for (int j=0; j<a.Width(); j++)
                    VT(i,j) = SVD_A.matrixV()(j,i);

            for (int i=0; i<Sigma.Size(); i++)
                Sigma(i) = SVD_A.singularValues()(i);

#endif

        }


        void Apply(const Vector<double> &b, Vector<double> &x, double factor = 1.) const override
        {
            if (b.Size() != a.Height())
                throw "Error in SVDSolver::Apply: matrix dimensions don't fit";
                
            x.SetSize(a.Width());
            x.SetAll(0);

            int max_hw = std::max(a.Height(), a.Width());

            // help-vector, assign maximum needed memory and set to zero
            Vector<> Uinvb(max_hw);
            Uinvb.SetAll(0.);

            // Apply U^T (or U^H for complex systems)
            U.ApplyT(b, Uinvb);

            // divide by singular values, if absolute value is large enough
            for (int i=0; i<Sigma.Size(); i++)
                if (Sigma(i) > 1e-10) Uinvb(i) /= Sigma(i);
                else Uinvb(i) = 0.;

            // set help-vector to size as needed
            // memory is kept as allocated with maximum size above
            // additional entries are set to zero above
            Uinvb.SetSize(a.Width());

            // Apply VT^T (=V, but we have VT from LAPACK)
            VT.ApplyT(Uinvb, x);

            x *= factor;

         }

        void Print(std::ostream& os) const override
        {
            os << "U = " << U << std::endl;
            os << "VT = " << VT << std::endl;
            os << "S = " << Sigma << std::endl;
            
        }
    };



}