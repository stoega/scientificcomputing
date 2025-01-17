#include <iostream>
#include <vector>
#include <cmath>
#include "SCTridiagSparseMatrix.h"
#include "SCGaussSeidelSolver.h"
#include "SCTridiagLUSolver.h"
#include "SCCGSolver.h"
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
    TridiagGaussSeidelSolver<double> gsSolver(A);

    // https://en.cppreference.com/w/cpp/chrono/high_resolution_clock/now
    auto start = std::chrono::high_resolution_clock::now();
    gsSolver.Apply(rhs, solution);
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> diff = end - start;
    return std::chrono::duration<double>(end - start).count();
}

double timeCGSolver(TridiagSparseMatrix<double> &A, Vector<double> &rhs, Vector<double> &solution)
{
    CGSolver<double> cgSolver(A);

    // https://en.cppreference.com/w/cpp/chrono/high_resolution_clock/now
    auto start = std::chrono::high_resolution_clock::now();
    cgSolver.Apply(rhs, solution);
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> diff = end - start;
    return std::chrono::duration<double>(end - start).count();
}

int main()
{
    // Problem parameters
    const double k = 10.0; // Diffusion coefficient
    const double b = 0; // Convection coefficient

    // Starting size
    int n = 100;
    const int maxSize = 838860799; // Max Value before int overflow
    const double timeLimit = 120.0; // 2 minutes as time limit

    // Timing variables
    double luTime = 0.0;
    double gsTime = 0.0;
    double cgTime = 0.0;

    // Open output file for plotting data
    std::ofstream outFile("Timing_Data.csv");
    outFile << "N,LU_Time,GS_Time,CG_Time\n";

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

        // Fill the tridiagonal matrix according to equations
#pragma omp parallel for
        for (int i = 0; i < n; ++i)
        {
            // Diagonal elements
            A.Set(i, i, diagEntry);

            // Subdiagonal elements (left)
            if (i > 0)
            {
                A.Set(i, i - 1, lowerEntry);
            }

            // Superdiagonal elements (right)
            if (i < n - 1)
            {
                A.Set(i, i + 1, upperEntry);
            }
        }

        // Fill rhs vector
        rhs->SetAll(0.0);
        for (int i = n / 2; i < n; ++i)
        { // f(x) = 1 for x >= 0.5
            (*rhs)(i) = 1.0;
        }

        // Time solvers and write to out file
        std::cout << "N = " << n << ":\n";
        outFile << n << ',';
        if (luTime < timeLimit)
        {
            solution->SetAll(0.0);
            luTime = timeLUSolver(A, *rhs, *solution);
            outFile << luTime;
            std::cout << "\tLU time: " << luTime << "s\n";

        }
        outFile << ',';

        if (gsTime < timeLimit)
        {
            solution->SetAll(0.0);
            gsTime = timeGSSolver(A, *rhs, *solution);
            outFile << gsTime;
            std::cout << "\tGS time: " << gsTime << "s\n";
        }
        outFile << ',';

        if (cgTime < timeLimit)
        {
            solution->SetAll(0.0);
            cgTime = timeCGSolver(A, *rhs, *solution);
            outFile << cgTime;
            std::cout << "\tCG time: " << cgTime << "s\n\n";

        }
        outFile << "\n";
        outFile.flush();

        // If every solver is above time limit -> break out of while
        if (luTime > timeLimit && gsTime > timeLimit && cgTime > timeLimit)
        {
            break;
        }

        n *= 2;
    }

    delete rhs;
    delete solution;
    outFile.close();
    return 0;
}
