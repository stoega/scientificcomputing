#pragma once

#include <iostream>
#include <complex>

#include "SCvector.h"
#include "SCmatrix.h"

namespace SC
{
    /*
    * LU factorization with pivoting
    * essentially from Numerical Recipies, Press et al. 2007
    */


    template <typename T = double>
    class LUSolver : public LinearOperator<T>
    {
    private:
        using LinearOperator<T>::height;
        using LinearOperator<T>::width;

    protected:
        const Matrix<T> & a;
        Matrix<T> * LU;
        Vector<int> * indx;

    public:

        LUSolver() = delete;

        LUSolver(const LUSolver &) = delete;
        LUSolver &operator=(const LUSolver &) = delete;
        
        LUSolver(const Matrix<T> & A) 
        : LinearOperator<T>(A.Height()), a(A), LU(nullptr), indx(nullptr)
        {

            if (!A.IsSquare())
                throw "Error in LUSolver: Matrix A must be square matrix";
            
            if (A.Height() > 0)
            {
                LU = new Matrix<T>(A);
                indx = new Vector<int>(A.Height());
            }

            ComputeFactorization();
        }

        ~LUSolver() 
        {
            if (LU != nullptr) delete LU;
            if (indx != nullptr) delete indx;
        }

        void ComputeFactorization()
        {
            const double TINY=1.0e-40; // A small number.
            int imax;
            double big;
            double temp;
            T temp_T;
            Vector<double> vv(height); //vv stores the implicit scaling of each row.

            (*LU) = a;

            for (int i=0;i<height;i++) 
            { //Loop over rows to get the implicit scaling information. 
                big=0.0;
                for (int j=0;j<height;j++)
                {
                    big = std::max( std::abs((*LU)(i,j)), big);
                }
                if (big == 0.0) throw("Singular matrix in LUdcmp");
                // No nonzero largest element.
                vv(i)=1.0/big; //Save the scaling.
            }

            for (int k = 0; k < height; k++)
            {              // This is the outermost kij loop.
                big = 0.0; // Initialize for the search for largest pivot element.
                for (int i = k; i < height; i++)
                {
                    temp = vv(i) * std::abs((*LU)(i, k));
                    if (temp > big)
                    { // Is the figure of merit for the pivot better than the best so far?
                        big = temp;
                        imax = i;
                    }
                }
                if (k != imax)
                { // Do we need to interchange rows?
                    for (int j = 0; j < height; j++)
                    { // Yes, do so...
                        std::swap((*LU)(imax, j), (*LU)(k, j));
                    }
                    // d   = -d; ...and change the parity of d.
                    vv(imax) = vv(k); // Also interchange the scale factor.
                }
                (*indx)(k) = imax;
                if ((*LU)(k, k) == 0.0)
                    throw("LUSolver: matrix exactly singular");
                // If the pivot element is zero, the matrix is singular and 1/piv cannot be computed
                for (int i = k + 1; i < height; i++)
                {
                    (*LU)(i, k) /= (*LU)(k, k); // Divide by the pivot element.
                    temp_T = (*LU)(i, k);
                    for (int j = k + 1; j < height; j++) // Innermost loop: reduce remaining submatrix.
                        (*LU)(i, j) -= temp_T * (*LU)(k, j);
                }
            }
        }


        void Apply(const Vector<T> &b, Vector<T> &x, T factor = 1.) const override
        {
            int ii=0;
            int ip;
            T sum;
            x.SetSize(height);
            if (b.Size() != height || x.Size() != height)
                throw("LUdcmp::solve bad sizes");
            
            x = b;
            x *= factor;

            for (int i=0;i<height;i++) 
            { // When ii is set to a positive value, it will become the
                //index of the first nonvanishing element of b. We now
                //do the forward substitution, equation (2.60). The
                //only new wrinkle is to unscramble the permutation
                //as we go.
                ip=(*indx)(i);
                sum=x(ip);
                x(ip)=x(i);
                if (ii != 0)
                {
                    for (int j=ii-1;j<i;j++) sum -= (*LU)(i,j)*x(j);
                }
                else if (sum != 0.0) // A nonzero element was encountered, so from now on we
                    //will have to do the sums in the loop above. 
                    ii=i+1;
                x(i)=sum;
            }
            for (int i=height-1;i>=0;i--) 
            {// Now we do the backsubstitution, equation (2.61).
                sum=x(i);
                for (int j=i+1;j<height;j++) sum -= (*LU)(i,j)*x(j);
                x(i)=sum/(*LU)(i,i); // Store a component of the solution vector X.
            }
        }

        void Print(std::ostream& os) const override
        {
            os << "L = [";
            for (int i=0; i<height; i++)
            {
                for (int j=0; j<i; j++)
                    os << (*LU)(i,j) << " ";
                os << "1.0" << std::endl;
            }
            os << "]" << std::endl;
            os << "U = [";
            for (int i=0; i<height; i++)
            {
                for (int j=0; j<i; j++)
                    os << "0 ";                
                for (int j=i; j<height; j++)
                    os << (*LU)(i,j) << " ";
                os <<  std::endl;
            }
            os << "]" << std::flush;
        }

        void PrintPattern(std::ostream& os) const
        {
            os << "L = [\n";
            for (int i=0; i<height; i++)
            {
                for (int j=0; j<i; j++)
                    if (std::abs((*LU)(i,j)) > 1e-16 ) os << "* ";
                    else os << "0 ";
                os << "1 " << std::endl;
            }
            os << "]" << std::endl;
            os << "U = [\n";
            for (int i=0; i<height; i++)
            {
                for (int j=0; j<i; j++)
                    os << "0 ";                
                for (int j=i; j<height; j++)
                    if (std::abs((*LU)(i,j)) > 1e-16 ) os << "* ";
                    else os << "0 ";
                os <<  std::endl;
            }
            os << "]" << std::flush;
        }
    };


}
