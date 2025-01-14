#include <iostream>
#include <vector>
#include <cmath>
#include "SCTridiagSparseMatrix.h" 
#include "SCTridiagLUSolver.h"       
#include "SCvector.h"

using namespace SC;

int main()
{
    // Problem parameters
    const double k = 10.0;          // Diffusion coefficient
    const double b = 50.0;          // Convection coefficient
    const int n = 10;              // Number of unknowns (discretization points)
    const double h = 1.0 / (n + 1); // Grid spacing

    // Create tridiagonal matrix
    TridiagSparseMatrix<double> A(n);

    // Fill the tridiagonal matrix according to equations
    for (int i = 0; i < n; ++i)
    {
        // Diagonal elements
        A.Set(i, i, 2.0 * k / (h * h) + b / h);

        // Subdiagonal elements (left)
        if (i > 0)
        {
            A.Set(i, i - 1, -k / (h * h) - b / h);
        }

        // Superdiagonal elements (right)
        if (i < n - 1)
        {
            A.Set(i, i + 1, -k / (h * h));
        }
    }

    // Create right-hand side vector b
    Vector<double> rhs(n);
    rhs.SetAll(0.0);
    for (int i = n / 2; i < n; ++i)
    { // f(x) = 1 for x >= 0.5
        rhs(i) = 1.0;
    }

    // Solve the system Ax = b using TridiagLUSolver
    TridiagLUSolver<double> solver(A);
    Vector<double> solution(n);
    solution.SetAll(0.0);

    solver.Apply(rhs, solution);

    // Output the solution
    std::cout << "Solution u(x):" << std::endl;
    for (int i = 0; i < n; ++i)
    {
        double xi = (i + 1) * h;
        std::cout << "x = " << xi << ", u(x) = " << solution(i) << std::endl;
    }

    return 0;
}
