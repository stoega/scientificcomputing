namespace SC
{
    template <typename T = double>
    class TridiagGaussSeidelSolver : public LinearOperator<T>
    {
    private:
        const TridiagSparseMatrix<T> &matrix;
        double acc;
        int maxit;

    public:
        // Constructors
        TridiagGaussSeidelSolver() = delete;
        TridiagGaussSeidelSolver(const TridiagSparseMatrix<T> &a, size_t maxIterations = 1000, double accuracy = 1e-6)
            : LinearOperator<T>(a.Height(), a.Width()),
              matrix(a),
              acc(accuracy),
              maxit(maxIterations) {}

        // Delete constructor
        ~TridiagGaussSeidelSolver() = default;

        // Copy constructor and assignment operator deleted
        TridiagGaussSeidelSolver(const TridiagGaussSeidelSolver &other) = delete;
        TridiagGaussSeidelSolver &operator=(const TridiagGaussSeidelSolver &other) = delete;

        void Apply(const Vector<T> &b, Vector<T> &x, T factor = 1.) const override
        {
            int n = matrix.Height();
            Vector<T> x_prev(n);

            // Formula 2.116 (https://en.wikipedia.org/wiki/Gauss%E2%80%93Seidel_method - Element-based formula)
            for (int iter = 0; iter < maxit; ++iter)
            {
                x_prev = x; // Store current values

                for (int i = 0; i < n; ++i)
                {
                    T sum = 0.0;

                    // Lower triangular part
                    if (i > 0)
                    {
                        sum += matrix.Get(i, i - 1) * x(i - 1); // Sub-diagonal
                    }

                    // Upper triangular part
                    if (i < n - 1)
                    {
                        sum += matrix.Get(i, i + 1) * x_prev(i + 1); // Super-diagonal
                    }

                    // Update x[i]
                    x(i) = (factor * b(i) - sum) / matrix.Get(i, i);
                }

                // Check convergence
                T error = 0.0;
                for (int i = 0; i < n; ++i)
                {
                    error += std::abs(x(i) - x_prev(i));
                }
#ifndef NDEBUG
                std::cout << "Iteration: " << iter << "\tResidual = " << error << std::endl;
#endif

                if (error < acc)
                {
#ifndef NDEBUG
                    std::cout << "Error < Acc" << std::endl;
#endif
                    return; // Convergence achieved
                }
            }
        }
    };
}