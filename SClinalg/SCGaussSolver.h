#pragma once

#include <iostream>
#include <complex>

#include "SCvector.h"
#include "SCmatrix.h"

namespace SC
{

    /**
     * Gauss solver, textbook version from semester 1, no pivoting
     * */
    
    template <typename T = double>
    class SimpleGaussSolver : public LinearOperator<T>
    {
    private:

    protected:
        const Matrix<T> & mat;

    public:
        using LinearOperator<T>::height;
        using LinearOperator<T>::width;

        SimpleGaussSolver() = delete;

        SimpleGaussSolver(const SimpleGaussSolver &) = delete;
        SimpleGaussSolver &operator=(const SimpleGaussSolver &) = delete;
        
        SimpleGaussSolver(const Matrix<T> & A) 
        : LinearOperator<T>(A.Height()), mat(A)
        {

            if (!A.IsSquare())
                throw "Error in GaussSolver: Matrix A must be square matrix";

        }

        ~SimpleGaussSolver() 
        { ; }

        void Apply(const Vector<T> &b, Vector<T> &x, T factor = 1.) const override
        {
            // copy matrix a, this copy is destroyed/turned to identity
            Matrix<T> temp(mat);

            x = b;
            x *= factor;

            T fac, pivinv;

            for (int i = 0; i < height; i++)
            {
                if (temp(i, i) == 0.0)
                    throw("SimpleGaussSolver: zero pivot element");
                pivinv = 1.0 / temp(i, i);
                for (int l = 0; l < height; l++)
                {
                    temp(i, l) *= pivinv;
                }
                x(i) *= pivinv;

                for (int ll = 0; ll < height; ll++)
                    if (ll != i)
                    {
                        fac = temp(ll, i);
                        for (int l = 0; l < height; l++)
                            temp(ll, l) -= temp(i, l) * fac;
                        x(ll) -= fac * x(i);
                    }
            }
            // std::cout << "internal temp matrix at end of solution process" << temp << std::endl;
        }

        // missing:
        // void ApplyT(const Vector<T> &b, Vector<T> &x, T factor = 1.) const override
        // void ApplyH(const Vector<T> &b, Vector<T> &x, T factor = 1.) const override

    };

    /*
    * Gauss solver: with pivoting, computes inverse by Gaussian elimination and applies inverse afterwards
    * essentially from Numerical Recipies, Press et al. 2007
    */


    template <typename T = double>
    class GaussSolver : public LinearOperator<T>
    {
    private:

    protected:
        Matrix<T> invmat;

    public:
        using LinearOperator<T>::height;
        using LinearOperator<T>::width;

        GaussSolver() = delete;

        GaussSolver(const GaussSolver &) = delete;
        GaussSolver &operator=(const GaussSolver &) = delete;
        
        GaussSolver(const Matrix<T> & A) 
        : LinearOperator<T>(A.Height()), invmat(A)
        {

            if (!A.IsSquare())
                throw "Error in GaussSolver: Matrix A must be square matrix";

            int i, icol, irow, j, k, l, ll, n = height;
            T dum, pivinv;
            double big;
            Vector<int> indxc(n), indxr(n), ipiv(n); // These integer arrays are used for bookkeeping on
            // the pivoting.
            ipiv.SetAll(0);

            for (i = 0; i < n; i++)
            { // This is the main loop over the columns to be
                // reduced.
                big = 0.0;
                for (j = 0; j < n; j++) // This is the outer loop of the search for a pivot
                    // element.
                    if (ipiv(j) != 1)
                        for (k = 0; k < n; k++)
                        {
                            if (ipiv(k) == 0)
                            {
                                if (std::abs(invmat(j, k)) >= big)
                                {
                                    big = std::abs(invmat(j, k));
                                    irow = j;
                                    icol = k;
                                }
                            }
                        }
                ++(ipiv(icol));

                if (irow != icol)
                {
                    for (l = 0; l < n; l++)
                        std::swap(invmat(irow, l), invmat(icol, l));
                }
                indxr(i) = irow;

                indxc(i) = icol;
                if (invmat(icol, icol) == 0.0)
                    throw("gaussj: Singular Matrix");
                pivinv = 1.0 / invmat(icol, icol);
                invmat(icol, icol) = 1.0;
                for (l = 0; l < n; l++)
                    invmat(icol, l) *= pivinv;
                for (ll = 0; ll < n; ll++) // Next, we reduce the rows...
                    if (ll != icol)
                    { //...except for the pivot one, of course.
                        dum = invmat(ll, icol);
                        invmat(ll, icol) = 0.0;
                        for (l = 0; l < n; l++)
                            invmat(ll, l) -= invmat(icol, l) * dum;
                    }
            }
            for (l = n - 1; l >= 0; l--)
            {
                if (indxr(l) != indxc(l))
                    for (k = 0; k < n; k++)
                        std::swap(invmat(k, indxr(l)), invmat(k, indxc(l)));
            } // And we are done.

        }

        ~GaussSolver() 
        { 
            ;
        }

        const Matrix<T>& GetInverse() const
        {
            return invmat;
        }

        void Apply(const Vector<T> &b, Vector<T> &x, T factor = 1.) const override
        {
            invmat.Apply(b, x, factor);
        }

        void ApplyT(const Vector<T> &b, Vector<T> &x, T factor = 1.) const override
        {
            invmat.ApplyT(b, x, factor);
        }

        void ApplyH(const Vector<T> &b, Vector<T> &x, T factor = 1.) const override
        {
            invmat.ApplyH(b, x, factor);
        }

    };


}