#include <iostream>
#include <vector>
#include <cmath>
#include "SCTridiagSparseMatrix.h"
#include "SCGaussSeidelSolver.h"
#include "SCvector.h"
#include <iostream>
#include <fstream>

using namespace SC;

/// @brief Generates a .csv file for given u and x vector
/// @param filename Name of output file with extension embedded
/// @param u First output vector
/// @param x Second output vector
void write_to_csv(const std::string &filename, const SC::Vector<double> &u, const SC::Vector<double> &x)
{
    // Check dimensions of vectors
    if (u.Size() != x.Size())
    {
#ifndef NDEBUG
        throw std::out_of_range("Error: Vectors u and x must have the same size.");
#endif
        std::cerr << "Error: Vectors u and x must have the same size." << std::endl;
        return;
    }

    std::ofstream file(filename);

    if (!file.is_open())
    {
#ifndef NDEBUG
        throw std::out_of_range("Error opening file for writing: ");
#endif
        std::cerr << "Error opening file for writing: " << filename << std::endl;
        return;
    }

    // Write headers
    file << "x,u\n";

    // Add values for -1 index as stated in task for visualisation
    file << 0 << "," << 0 << "\n";

    for (size_t i = 0; i < u.Size(); ++i)
    {
        file << x(i) << "," << u(i) << "\n";
    }

    // Add values for n index as stated in task for visualisation
    file << 1 << "," << 0 << "\n";

    file.close();
    std::cout << "Data written to " << filename << std::endl;
}

int main()
{
    // Problem parameters
    const double k = 10.0;             // Diffusion coefficient
    const double bArr[] = {0.0, 50.0}; // Convection coefficient
    const int n = 100;                 // Number of unknowns (discretization points)
    const double h = 1.0 / (n + 1);    // Grid spacing

    // Simulation for b=0 and b = 50 as well as saving data
    for (int idx = 0; idx < 2; idx++)
    {
        const double b = bArr[idx];
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
        TridiagGaussSeidelSolver<double> solver(A, 1e6, 1e-9);
        Vector<double> solution(n);
        solution.SetAll(0.0);

        solver.Apply(rhs, solution);

        // Output the solution
        Vector<double> x_i(n);  // linspace for x-axis
        std::cout << "Solution u(x):" << std::endl;
        for (int i = 0; i < n; ++i)
        {
            x_i(i) = (i + 1) * h;
            std::cout << "x = " << x_i(i) << ", u(x) = " << solution(i) << std::endl;
        }

        // Constructiong filename and writing to csv
        std::ostringstream filename_stream;
        filename_stream << "AE_b_" << (int)b << ".csv";
        std::string filename = filename_stream.str();
        write_to_csv(filename, solution, x_i);
    }

    return 0;
}
