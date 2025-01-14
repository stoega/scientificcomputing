#include <iostream>
#include <vector>
#include <cmath>
#include "SCTridiagSparseMatrix.h"
#include "SCGaussSeidelSolver.h"
#include "SCTridiagLUSolver.h"
#include "SCvector.h"
#include <iostream>
#include <fstream>
#include <chrono>

using namespace SC;

double timeLUSolver(TridiagSparseMatrix<double> &A, Vector<double> &rhs, Vector<double> &solution)
{
    TridiagLUSolver<double> luSolver(A);

    // https://en.cppreference.com/w/cpp/chrono/high_resolution_clock/now
    auto start = std::chrono::high_resolution_clock::now();
    luSolver.Apply(rhs, solution);
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> diff = end - start;
    return std::chrono::duration<double>(end - start).count();
}

double timeGSSolver(TridiagSparseMatrix<double> &A, Vector<double> &rhs, Vector<double> &solution)
{
    TridiagLUSolver<double> gsSolver(A);

    // https://en.cppreference.com/w/cpp/chrono/high_resolution_clock/now
    auto start = std::chrono::high_resolution_clock::now();
    gsSolver.Apply(rhs, solution);
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> diff = end - start;
    return std::chrono::duration<double>(end - start).count();
}

int main()
{
    // Problem parameters
    const double k = 10.0; // Diffusion coefficient
    const double b = 50.0; // Convection coefficient

    // Starting size
    int n = 100;
    const int maxSize = 838860799;  // Adjust based on your system
    const double timeLimit = 180.0; // 3 minutes in seconds

    // Timing variables
    double luTime = 0.0;
    double gsTime = 0.0;

    // Open output file for plotting data
    std::ofstream outFile("Timing_Data.csv");
    outFile << "N,LU_Time,GS_Time\n";

    // Reusable vectors declared outside the loop
    Vector<double> *solution = nullptr;
    Vector<double> *rhs = nullptr;

    while (n <= maxSize)
    {
        delete solution;
        delete rhs;

        const double h = 1.0 / (n + 1); // Grid spacing
        const double diagEntry = 2.0 * k / (h * h) + b / h;
        const double lowerEntry = -k / (h * h) - b / h;
        const double upperEntry = -k / (h * h);

        // Resize vectors instead of recreating them
        solution = new Vector<double>(n);
        solution->SetAll(0.0);
        rhs = new Vector<double>(n);
        rhs->SetAll(0.0);

        // Create tridiagonal matrix
        TridiagSparseMatrix<double> A(n);
        // Vector<double> U(n);
        // Vector<double> D(n);
        // Vector<double> L(n);

        // U.SetAll(upperEntry);
        // D.SetAll(diagEntry);
        // L.SetAll(lowerEntry);

        // Fill the tridiagonal matrix according to equations
#pragma omp parallel for
        for (int i = 0; i < n; ++i)
        {
            // Diagonal elements
            A.Set(i, i, diagEntry);
            // D(i) = 2.0 * k / (h * h) + b / h;

            // Subdiagonal elements (left)
            if (i > 0)
            {
                A.Set(i, i - 1, lowerEntry);
                // L(i) = -k / (h * h) - b / h;
            }

            // Superdiagonal elements (right)
            if (i < n - 1)
            {
                A.Set(i, i + 1, upperEntry);
                // U(i) = -k / (h * h);
            }
        }
        // TridiagSparseMatrix<double> A(n, D, L, U);

        // Create right-hand side vector b
        // Vector<double> *rhs = new Vector<double>(n);
        rhs->SetAll(0.0);
        for (int i = n / 2; i < n; ++i)
        { // f(x) = 1 for x >= 0.5
            (*rhs)(i) = 1.0;
        }

        // Definition of solvers
        // TridiagGaussSeidelSolver<double> gsSolver(A);

        // Creating solution vector
        // Vector<double> *solution = new Vector<double>(n);
        solution->SetAll(0.0);

        luTime = timeLUSolver(A, *rhs, *solution);

        solution->SetAll(0.0);
        // start = std::chrono::high_resolution_clock::now();
        // gsSolver.Apply(rhs, solution);
        // end = std::chrono::high_resolution_clock::now();
        // gsTime = std::chrono::duration<double>(end - start).count();
        gsTime = timeGSSolver(A, *rhs, *solution);

        // Write to file
        outFile << n << " " << luTime << " " << gsTime << "\n";
        outFile.flush(); // Ensure data is written in case of crash

        // Print progress
        std::cout << "N = " << n << ":\n";
        std::cout << "  LU time: " << luTime << "s\n";
        std::cout << "  GS time: " << gsTime << "s\n\n";

        // Check time limit
        if (luTime > timeLimit || gsTime > timeLimit)
            break;

        // Increase size (LU solver: multiply by 10, GS solver: multiply by 2)
        // n *= (luTime > gsTime) ? 2 : 10;
        n *= 2;
    }
    delete rhs;
    delete solution;
    outFile.close();
    return 0;
}
