#include <SCvector.h>
#include <SCmatrix.h>
#include "SCHessenbergmatrix.h"

namespace SC
{
    template<typename T>
    void Arnoldi(const LinearOperator<T>& A, const Vector<T>& r, int k, Matrix<T> &V, HessenbergMatrix<T> &H, double acc=1e-8)
    {
        int n = A.Height();
        V.SetSize(n,k+1);
        V.SetStorage(TMatrixStorage::COL_MAJOR);
        H.SetSize(k);
        V.SetAll(0.);
        H.SetAll(0.);

        Vector<T> q(r), w(r);


        double normw = w.Norm();
        w *= 1./normw;
        for (int i=0; i<n; i++)
            V(i,0) = w(i);
        
        for (int j=1; j<=k; j++)
        {
            A.Apply(w,q);
            // (modified) Gram Schmidt orthogonalization, reverse order
            for (int i=j-1; i>=0; i--)
            {
                // (q , v_i) = H_i,j-1
                for (int in=0; in<n; in++)
                {
                    H.AddTo(i,j-1, Conjugate(V(in,i))*q(in));
                }
                // q -= (q , v_i) v_i
                for (int in=0; in<n; in++)
                    q(in) -= H.Get(i,j-1)*V(in,i);
            }
            H.Set(j,j-1, q.Norm());

            if (j<k)
            {
                w = q;
                w *= 1./H.Get(j,j-1);
                
                for (int in=0; in<n; in++)
                    V(in,j) = w(in);
            }
            
            if (std::abs(H.Get(j,j-1)) < acc*std::abs(H.Get(1,0))) 
            {
                H.SetSize(j);
                V.SetSize(n,j+1);
                return;
            }

        }



    }


}