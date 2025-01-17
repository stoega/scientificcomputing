#pragma once

#include <iostream>
#include <complex>

#include "SCvector.h"
#include "SCmatrix.h"

namespace SC
{
    template <typename T = double>
    class CGSolver : public LinearOperator<T>
    {
        using LinearOperator<T>::height;
        using LinearOperator<T>::width;
    private:


    protected:
        const LinearOperator<T> &a;
        int maxit;
        double acc;

    public:        
        CGSolver() = delete;
        
        CGSolver(const LinearOperator<T> &A, int maxiterations=-1, double accuracy=1e-8) 
        : LinearOperator<T>(A.Height()) 
        , a(A)
        , maxit(maxiterations)
        , acc(accuracy)
        {
            if (maxit < 0) maxit = A.Height();
            if (!A.IsSquare())
                throw "Error in CGSolver: Matrix A must be square matrix";
        }

        ~CGSolver() {}

        CGSolver(const CGSolver &) = delete;
        CGSolver &operator=(const CGSolver &) = delete;


        void Apply(const Vector<T> &b, Vector<T> &x, T factor = 1.) const override
        {
            Vector<T> res(b);
            Vector<T> p(b);
            Vector<T> Ap(b);
            T r0;
            T r2;
            T r2_prev;

            T alpha, beta;

            if (x.Size() != height)
            {
                x.SetSize(height);
                x.SetAll(0.);
            }

            // compute the initial residual res = b - a * x
            a.Apply(x,res);
            res *= -1;
            res += b;
    
            r0 = InnerProduct(res,res);
            r2 = r0;

#ifndef NDEBUG
            std::cout << "initial res " << sqrt(r0) << std::endl;
#endif



            for (int k=0; k<maxit; k++)
            {
                // Ap = A * p
                a.Apply(p, Ap);
                
                alpha = r2 / InnerProduct(p, Ap);

                // x += alpha * p
                // res -= alpha * Ap
                x.AddMultiple(alpha, p);
                res.AddMultiple(-alpha, Ap);

                r2_prev = r2;
                r2 = InnerProduct(res,res);
                beta = r2/r2_prev;

                // p = r + beta p
                p *= beta;
                p.Add(res);

#ifndef NDEBUG
                std::cout << "step " << k << " res " << sqrt(std::abs(r2)) << std::endl;
#endif

                if (sqrt(std::abs(r2)) < acc*sqrt(std::abs(r0))) break;

            }

            x *= factor;

        }

    };
}
