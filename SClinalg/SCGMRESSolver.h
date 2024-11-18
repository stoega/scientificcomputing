#pragma once

#include <iostream>
#include <complex>

#include "SCvector.h"
#include "SCmatrix.h"

#include "SCArnoldi.h"
#include "SCQRSolver_Hessenberg.h"

namespace SC
{
    template <typename T = double>
    class GMRESSolver : public LinearOperator<T>
    {
        using LinearOperator<T>::height;
        using LinearOperator<T>::width;
    private:


    protected:
        const LinearOperator<T> &a;
        int maxit;
        int k_r;
        double acc;


    public:        
        GMRESSolver() = delete;
        
        GMRESSolver(const LinearOperator<T> &A, int maxiterations=-1, int k_restart = 5, double accuracy=1e-6) 
        : LinearOperator<T>(A.Height()) 
        , a(A)
        , maxit(maxiterations)
        , acc(accuracy)
        , k_r(k_restart)
        {
            if (maxit < 0) maxit = A.Height();
            if (!A.IsSquare())
                throw "Error in GMRESSolver: Matrix A must be square matrix";
        }

        ~GMRESSolver() {}

        GMRESSolver(const GMRESSolver &) = delete;
        GMRESSolver &operator=(const GMRESSolver &) = delete;


        void Apply(const Vector<T> &b, Vector<T> &x, T factor = 1.) const override
        {
            Vector<T> res(b), help(b);
            
            if (x.Size() != height)
            {
                x.SetSize(height);
                x.SetAll(0.);
            }

            // compute the initial residual res = b - a * x
            a.Apply(x,res);
            res *= -1;
            res += b;

            double error = res.Norm();
            double initial_res = error;
            std::cout << "initial res " << initial_res << std::endl;

            Matrix<T> V;
            HessenbergMatrix<T> H;

            int iterations = 0;
            while (error > initial_res * acc && iterations < maxit)
            {
                // Arnoldi iteration -- compute V^(k+1), H^(k+1)
                Arnoldi(a, res, k_r, V, H, acc);

                // comupute minimizer y for
                // || ||r|| e_0 - H^(k+1) y ||
                // solving normal equation via QR decomposition
                // H = QR
                // R y = ||r|| Q^H e_0 
                QRSolverHessenberg<T> QR_H(H);

                Vector<T> y;
                Vector<T> e0(H.Height());
                e0.SetAll(0.);
                e0(0) = res.Norm();
                QR_H.Apply(e0, y);


                // compute new iterate x^(k) = x^(0) + V y
                for (int in=0; in<a.Height(); in++)
                    for (int j=0; j<V.Width()-1; j++)
                        x(in) += V(in,j)*y(j);
                
                iterations += V.Width()-1;
                
                // update the residual
                a.Apply(x,res);
                res *= -1;
                res += b;

                error = res.Norm();
                std::cout << "it " << iterations << ": error " << error << std::endl;
            
            }
            x *= factor;
            
        }

    };
}