#pragma once

// #include "SCmatrix.h"
// #include "SCTridiagSparseMatrix.h"
// #include "SCSparseMatrix.h"

namespace SC
{
    template <typename T = double>
    class TridiagLUSolver : public LinearOperator<T>
    {
    private:
        const TridiagSparseMatrix<T> &matrix;
        // std::vector<std::vector<double>> LU;    // since L has 1 as diagonal, we can store L and U in one matrix
        Vector<T> L;
        Vector<T> D;
        Vector<T> U;

    public:
        TridiagLUSolver() = delete;

        TridiagLUSolver(const TridiagSparseMatrix<T> &mat)
            : LinearOperator<T>(mat.Height(), mat.Width()), 
              matrix(mat), 
              L(mat.Height() - 1), 
              D(mat.Height()), 
              U(mat.Height() - 1)
        {
            // LU.resize(mat.Height(), std::vector<double>(mat.Width(), 0.0));
            Factorize();
        }

        // Copy constructor and assignment operator deleted
        TridiagLUSolver(const TridiagLUSolver &other) = delete;
        TridiagLUSolver &operator=(const TridiagLUSolver &other) = delete;

        ~TridiagLUSolver() = default;

        void Factorize()
        {
            int n = matrix.Height();

            // Perform LU decomposition
            // for (int i = 0; i < n; i++)
            // {
            //     // Calculate upper triangular matrix (U)
            //     for (int k = i; k < n; k++)
            //     {
            //         T sum;
            //         for (int j = 0; j < i; j++)
            //             sum += LU[i][j] * LU[j][k];
            //         LU[i][k] = matrix.Get(i, k) - sum;
            //     }

            //     // Calculate lower triangular matrix (L)
            //     for (int k = i; k < n; k++)
            //     {
            //         T sum = 0.0;
            //         for (int j = 0; j < i; j++)
            //             sum += LU[k][j] * LU[j][i];
            //         LU[k][i] = (matrix.Get(k, i) - sum) / LU[i][i];
            //     }
            // }

            // First row
            D(0) = matrix.Get(0, 0);
            U(0) = matrix.Get(0, 1);

            // Middle rows
            for (int i = 1; i < n - 1; ++i)
            {
                L(i - 1) = matrix.Get(i, i - 1) / D(i - 1); // l(i, j) = a(i, j) / u(j, j)
                D(i) = matrix.Get(i, i) - L(i - 1) * U(i - 1);  // u(i, j) -= l(i, k) * u(k, j)
                U(i) = matrix.Get(i, i + 1);    // u(i, j) = a(i, j)
            }

            // Last row
            if (n > 1)
            {
                L(n - 2) = matrix.Get(n - 1, n - 2) / D(n - 2);
                D(n - 1) = matrix.Get(n - 1, n - 1) - L(n - 2) * U(n - 2);
            }
        }

        void Apply(const Vector<T> &x, Vector<T> &r, T factor = 1.) const override
        {
            int n = matrix.Height();
            r.SetSize(n);
            Vector<T> y(n); // Temporary vector for forward substitution

            // Psuedo-Code: sums can be reduced to single statement since tridiagform
            // Forward substitution (Ly = x)
            y(0) = x(0);
            for (int i = 1; i < n; i++)
            {
                // T sum;
                // for (int j = 0; j <= i - 1; j++)
                // {
                //     // sum += LU.Get(i, j) * y(j);
                //     sum += LU[i][j] * y(j);
                // }
                // y(i) = x(i) - sum;
                y(i) = x(i) - L(i - 1) * y(i - 1);
            }

            // Backward substitution (Ur = y)
            // r(n - 1) = y(n - 1) / LU[n - 1][n - 1];
            r(n - 1) = y(n - 1) / D(n - 1);

            for (int i = n - 2; i >= 0; i--)
            {
                // T sum;
                // for (int j = i + 1; j < n; j++)
                // {
                //     sum += LU[i][j] * x(j);
                // }
                // r(i) = (y(i) - sum) / LU[i][i];
                r(i) = (y(i) - U(i) * r(i + 1)) / D(i);
            }

            // Apply scaling factor
            if (factor != T(1))
            {
                r *= factor;
            }
        }
    };
}
