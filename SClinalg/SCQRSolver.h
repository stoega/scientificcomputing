#pragma once

#include <iostream>
#include <complex>

#include "SCvector.h"
#include "SCmatrix.h"

namespace SC
{
    template <typename T = double>
    class QRSolver : public LinearOperator<T>
    {

    protected:
        const Matrix<T> & A;
        Matrix<T> W; // lower left triangular matrix containing Householder vectors
        Matrix<T> R; // upper right triangular matrix

        Vector<T> v; // help vector for computing factorization

    public:

        QRSolver() = delete;

        QRSolver(const QRSolver &) = delete;
        QRSolver &operator=(const QRSolver &) = delete;
        
        QRSolver(const Matrix<T> & A_) 
        : LinearOperator<T>(A_.Width(), A_.Height()), A(A_), R(A_.Height(),A_.Width()), W(A_.Height()), v(A_.Height())
        {
            ComputeFactorization();
        }

        ~QRSolver() 
        {
            ;
        }

        void ComputeFactorization()
        {
            
            // R contains a copy of A
            // W is empty
            R = A;
            W.SetAll(0.);

            for (int i=0;i<A.width;i++) 
            { 
                T sign;
                sign = (std::abs(R(i, i)) == T(0)) ? T(1) : R(i, i) / std::abs(R(i, i));

                double norm = 0;
                for (int j=i; j<A.height; j++) norm += std::norm(R(j,i)); // std::norm computes absolute value squared

                norm = sqrt(norm);
                for (int j=0; j<i; j++)
                    v(j) = 0;
                v(i) = R(i,i) + sign * norm;
                for (int j=i+1; j<A.height; j++)
                    v(j) = R(j,i);


                double normv = v.Norm();

                for (int j=i; j<A.height; j++)
                    W(j,i) = v(j)/normv;

                // multiply R by P
                for (int j=i; j<A.width; j++) // column j
                {
                    // inner product col R_j * col W_i
                    T ip = 0.;
                    for (int k=i; k<A.height; k++)
                        ip += R(k,j)*Conjugate(W(k,i));

                    for (int k=0; k<A.height; k++)
                        R(k,j) = R(k,j) - 2. * ip * W(k,i);

                }

            }

        }

        // Matrix R
        const Matrix<T>& GetR() const { return R; }

        // compute actual matrix Q -- for debugging reasons
        void ComputeQ(Matrix<T> &Q) const
        {
            Q.SetSize(A.height, A.height);
            Q.SetAll(0.);
            Q.SetDiag(1.);
            T ip;
            for(int col=0; col<A.height; col++)
            {
            for (int i=A.width-1; i>=0; i--)
            {
                ip = 0.;
                for (int j=i; j<A.height; j++)
                    ip += Conjugate(W(j,i))*Q(j,col);
                for (int j=i; j<A.height; j++)
                    Q(j,col) -= 2.*ip*W(j,i);
            }
            }
        }


        void ApplyQinv(Vector<T> &x) const
        {
            T ip;
            for (int i=0; i<A.width; i++)
            {
                ip = 0.;
                for (int j=i; j<A.height; j++)
                    ip += Conjugate(W(j,i))*x(j);
                for (int j=i; j<A.height; j++)
                    x(j) -= 2.*ip*W(j,i);
            }
        }


        void ApplyRinv(Vector<T> &x) const
        {
            T sum;
            for (int i=A.width-1; i>=0; i--)
            {   
                sum=x(i);
                for (int j=i+1;j<A.width;j++) sum -= R(i,j)*x(j);
                x(i)=sum/R(i,i); // Store a component of the solution vector X.
            }
        }

        void Apply(const Vector<T> &b, Vector<T> &x, T factor = 1.) const override
        {
            Vector<T> Qinvb(b);
            Qinvb *= factor;

            ApplyQinv(Qinvb);

            // std::cout << "QinvB = " << Qinvb << std::endl;

            x.SetSize(A.width);
            for (int i=0; i<A.width; i++)
                x(i) = Qinvb(i);

            ApplyRinv(x);

         }

        void Print(std::ostream& os) const override
        {
            os << "Householder vectors = [";
            for (int i=0; i<A.width; i++)
            {
                for (int j=0; j<A.height; j++)
                    os << W(i,j) << " ";
                os << "\n";
            }
            os << "]" << std::endl;
            os << "R = [";
            for (int i=0; i<A.height; i++)
            {
                // for (int j=0; j<i; j++)
                //     os << "0 ";                
                for (int j=0; j<A.width; j++)
                    os << R(i,j) << " ";
                os <<  std::endl;
            }
            os << "]" << std::flush;
        }

        // necessary for QR algorithm for computing eigenvalues
        void ComputeRQ(Matrix<T> &RQ) const
        {
    #ifndef NDEBUG
            if (!A.IsSquare()) std::cout << "Warning - computing RQ for non-square matrix" << std::endl;
    #endif

            RQ = R;
            T ip;
            for(int row=0; row<A.height; row++)
            {
                for (int i=0; i<A.width; i++)
                {
                    ip = 0.;
                    for (int j=i; j<A.height; j++)
                        ip += W(j,i)*RQ(row,j);
                    
                    for (int j=i; j<A.height; j++)
                        RQ(row,j) -= 2.*ip*Conjugate(W(j,i));
                }
            }
        }
    
    };



}
