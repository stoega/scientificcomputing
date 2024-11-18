#pragma once

#include <iostream>
#include <complex>

#include "SCvector.h"
#include "SCmatrix.h"
#include "SCHessenbergmatrix.h"

namespace SC
{
    template <typename T = double>
    class QRSolverHessenberg : public LinearOperator<T>
    {
    private:

    protected:
        const HessenbergMatrix<T> & a; // must be of Hessenberg form
        Matrix<T> W; // 2xa.height matrix containing householder vectors, only 2 non-zero entries (i,i) and (i+1,i) in each col.
        HessenbergMatrix<T> R; // upper right triangular matrix

    public:

        QRSolverHessenberg() = delete;

        QRSolverHessenberg(const QRSolverHessenberg &) = delete;
        QRSolverHessenberg &operator=(const QRSolverHessenberg &) = delete;
        
        QRSolverHessenberg(const HessenbergMatrix<T> & A) 
        : LinearOperator<T>(A.Width(), A.Height()), a(A), R(A), W(2,A.Width())
        {
            ComputeFactorization();
        }

        ~QRSolverHessenberg() 
        {
            ;
        }

        void ComputeFactorization()
        {
            
            // R contains a copy of A
            // W is empty
            R = a;
            W.SetAll(0.);

            T v[2]; // two non-zero components of each vector v

            for (int i=0;i<a.width;i++) 
            { 
                T sign;
                if(std::abs(R.Get(i,i))==0) sign = T(1);
                else sign = R.Get(i,i)/std::abs(R.Get(i,i));
                double norm = 0;
                for (int j=i; j<i+2; j++) norm += std::norm(R.Get(j,i));
                norm = sqrt(norm);

                v[0] = R.Get(i,i) + sign * norm;
                if (i+1<a.height) v[1] = R.Get(i+1,i);

                double normv = sqrt(std::norm(v[0]) + std::norm(v[1]));
                for (int j=0; j<2; j++)
                    if (i+j<a.height) W(j,i) = v[j]/normv;

                // multiply R by P
                for (int j=i; j<a.width; j++) // column j
                {
                    // inner product col R_j * col W_i
                    T ip = R.Get(i,j)*Conjugate(W(0,i));
                    if (i+1<a.height) ip += R.Get(i+1,j)*Conjugate(W(1,i));

                    // for (int k=0; k<height; k++)
                    //     R(k,j) = R(k,j) - 2. * ip * W(k,i);
                    R.AddTo(i,j, - 2. * ip * W(0,i));
                    if (i+1<a.height) R.AddTo(i+1,j, - 2. * ip * W(1,i));

                }

            }

        }

        const Matrix<T>& GetR() const { return R; }

        void ComputeRQ(HessenbergMatrix<T> &RQ) const
        {
            RQ = R;
            T ip;
            for(int row=0; row<a.width; row++)
            {
                int imin = row>0 ? row-1 : 0;
                for (int i=imin; i<a.width; i++)
                {
                    ip = (W(0,i))*RQ.Get(row,i);
                    if (i+1<a.width) ip += (W(1,i))*RQ.Get(row,i+1);

                    RQ.AddTo(row,i, -2.*ip*Conjugate(W(0,i)));
                    if (i+1<a.width) RQ.AddTo(row,i+1, -2.*ip*Conjugate(W(1,i)));
                }
            }
        }

        void ApplyQinv(Vector<T> &x) const
        {
            T ip;
            for (int i=0; i<a.width; i++)
            {
                ip = Conjugate(W(0,i))*x(i);
                if (i+1<a.height) ip += Conjugate(W(1,i))*x(i+1);

                x(i) -= 2.*ip*W(0,i);
                if (i+1<a.height) x(i+1) -= 2.*ip*W(1,i);
            }
        }


        void ApplyRinv(Vector<T> &x) const
        {
            T sum;
            for (int i=a.width-1; i>=0; i--)
            {   
                sum=x(i);
                for (int j=i+1;j<a.width;j++) sum -= R.Get(i,j)*x(j);
                x(i)=sum/R.Get(i,i); 
            }
        }

        void Apply(const Vector<T> &b, Vector<T> &x, T factor = 1.) const override
        {
            Vector<T> Qinvb(b);
            Qinvb *= factor;

            ApplyQinv(Qinvb);

            x.SetSize(a.width);
            for (int i=0; i<a.width; i++)
                x(i) = Qinvb(i);

            ApplyRinv(x);

         }

        void Print(std::ostream& os) const override
        {
            os << "Householder vectors = [";
            for (int i=0; i<2; i++)
            {
                for (int j=0; j<a.height; j++)
                    os << W(i,j) << " ";
                os << "\n";
            }
            os << "]" << std::endl;
            os << "R = " << R << std::endl;
        }
    };



}
