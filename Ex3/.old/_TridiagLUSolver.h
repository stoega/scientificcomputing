// #include <iostream>
// #include <complex>
// #include <type_traits>
// #include <iomanip>
// #include "SCmatrix.h"

namespace SC
{
    template <typename T = double>
    class TridiagLUSolver : public LinearOperator<T>
    {
    private:
        const TridiagSparseMatrix<T> &matrix;
        SC::SparseMatrix<T> L;
        SC::SparseMatrix<T> U;

    public:
        // Default constructor
        // TridiagLUSolver() = delete;

        // Constructor for given matrix
        TridiagLUSolver(const TridiagSparseMatrix<T> &mat) : matrix(mat){
            Factorize();
        }

        // Copy
        TridiagLUSolver(const TridiagLUSolver &other) = delete;

        // Destructor
        ~TridiagLUSolver() = default;

        // Assignment operator
        TridiagLUSolver &operator=(const TridiagLUSolver &other) = delete;

        void Factorize()
        {
            int n = matrix.Size();
            L = SC::SparseMatrix(n, n);
            // L.SetAll(0);
            U = SC::SparseMatrix(n, n);
            // U.SetAll(0);

            for (int j = 0; j < n; j++)
            {
                for (int i = 0; i <= j; i++)
                {
                    U(i, j) = A(i, j);
                    for (int k = 0; k < i; k++)
                    {
                        U(i, j) -= L(i, k) * U(k, j)
                    }
                }
                for (int i = j + 1; i <= n; i++)
                {
                    U(i, j) = A(i, j) / U(j, j);
                    for (int k = 0; k < j; k++)
                    {
                        L(i, j) -= L(i, k) * U(k, j) / U(j, j);
                    }
                }
            }
        }

        void Apply(const Vector<T> &x, Vector<T> &r, T factor = 1.) const override
        {
            int n = x.size();
            SC::Vector<T> y(n);

            // forward and backward substitution (2.60)
            // Forward Ly = b with y_0 = b_0
            y(0) = x(0);
            for (int i = 1; i < n; i++)
            {
                y(i) = x(i) - L(i - 1) * y(i - 1); // (2.61)
            }

            // Backward Ux = y with x_n-1 = y_n-1 / u_n-1,n-1
            r[n - 1] = y[n - 1] / U[n - 1, n - 1];
            for (int i = n - 2; i >= 0; i--)
            {
                r[i] = (y[i] - U[i] * r[i + 1]) / U[i];
            }

            // Scaling
            for(int i = 0; i < n; i++){
                r[i] *= factor;
            }
        }
    };
}