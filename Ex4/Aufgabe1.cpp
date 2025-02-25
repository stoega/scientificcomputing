#include <iostream>
#include <vector>
#include <cmath>
#include "SCTridiagSparseMatrix.h"
#include "SCTridiagLUSolver.h"
#include "SCvector.h"
#include <iostream>
#include <fstream>

using namespace SC;

// Normalize a vector (divide by its Euclidean norm)
template <typename T>
void normalize(SC::Vector<T> &v)
{
    T norm = 0.0;
    // for (const auto& val : v) {
    //     norm += val * val;
    // }
    norm = v.Norm();

    v.Mult(1 / norm);
}

void ComputeInverseIteration(const TridiagSparseMatrix<double> &K, const TridiagSparseMatrix<double> &M, const Vector<double> &u, double mu, Vector<double> &eigenVec, double &eigenValue, double tol = 1e-9, int maxIter = 1000)
{
    int n = K.Height();
    TridiagSparseMatrix<double> Kshifted(n);

    // Shift K with K - mu*M
    for (int i = 0; i < n; i++)
    {
        double shiftedValue = K.Get(i, i) - mu * M.Get(i, i);
        Kshifted.Set(i, i, shiftedValue);

        // Subdiagonal elements (left)
        if (i > 0)
        {
            shiftedValue = K.Get(i, i - 1) - mu * M.Get(i, i - 1);
            Kshifted.Set(i, i - 1, shiftedValue);
        }

        // Superdiagonal elements (right)
        if (i < n - 1)
        {
            shiftedValue = K.Get(i, i + 1) - mu * M.Get(i + 1, i);
            Kshifted.Set(i, i + 1, shiftedValue);
        }
    }

    // Define solver, also computes LU factorization
    TridiagLUSolver<double> solver(Kshifted);

    // Create solution vector
    Vector<double> z(n);
    z.SetAll(0.0);

    Vector<double> v(u);
    normalize(v);
    double lambda = 0.0;

    // Residual vector
    Vector<double> r(n);
    r.SetAll(0.0);

    // Temporary vectors used for calculations
    Vector<double> Kv(n);
    Vector<double> Mv(n);

    std::cout << "Starting inverse iteration with shift = " << mu << std::endl;
    for (int iter = 0; iter < maxIter; iter++)
    {
        // (3.32) left
        M.Apply(v, Mv);

        // Solve Kshifted * z = M*v
        solver.Apply(Mv, z);

        // Rayleigh quotient (3.33)
        Kshifted.Apply(z, Kv);
        M.Apply(z, Mv);
        double numerator = SC::InnerProduct(z, Kv);
        double denominator = SC::InnerProduct(z, Mv);
        lambda = numerator / denominator;

        for (int j = 0; j < n; j++)
        {
            r(j) = Kv(j) - lambda * Mv(j);
        }

        normalize(z);

        if (r.Norm() <= tol)
        {
            std::cout << "Converged after " << iter + 1 << " iterations." << std::endl;
            break;
        }

        // (3.32) right
        v = z;

        // Print every 10 iterations
        if (iter % 10 == 0)
        {
            std::cout << "Iteration " << iter
                      << ", eigenvalue = " << lambda
                      << ", residual = " << r.Norm() << std::endl;
        }
    }
    eigenValue = lambda;
    eigenVec = z;
}

int main()
{
    // Problem parameters
    const float L = .5;   // length of flute
    const int c = 343;    // speed of sound
    const int n = 100;    // number of unknowns
    const int ne = n + 1; // number of elements

    const double tol = 1e-9;  // error tolerance
    const int maxIter = 1000; // max iterations

    const double h = L / ne; // mesh fineness

    // Create tridiagonal matrix
    TridiagSparseMatrix<double> K(n);
    TridiagSparseMatrix<double> M(n);

    // Fill the tridiagonal matrices according to equation (5) and (6)
    for (int i = 0; i < n; ++i)
    {
        // Diagonal elements
        K.Set(i, i, 2.0 * pow(c, 2) / h);
        M.Set(i, i, 4.0 * h / 6.0);

        // Subdiagonal elements (left)
        if (i > 0)
        {
            K.Set(i, i - 1, -pow(c, 2) / h);
            M.Set(i, i - 1, h / 6.0);
        }

        // Superdiagonal elements (right)
        if (i < n - 1)
        {
            K.Set(i, i + 1, -pow(c, 2) / h);
            M.Set(i, i + 1, h / 6.0);
        }
    }

    // Create solution vector u for given mesh
    Vector<double> u(n);
    u.SetAll(0.0);
    for (int i = 0; i < n; ++i)
    {
        double xi = (i + 1) * L / ne;
        u(i) = sin(n * M_PI * xi / L);
    }

    /*
     *   Aufgabe 1
     */
    Vector<double> eigenVector1(n);
    double eigenValue1;
    ComputeInverseIteration(K, M, u, 0, eigenVector1, eigenValue1);

    double omega1 = std::sqrt(eigenValue1);
    double analyticalOmega1 = c * M_PI / L;
    double error1 = std::abs(omega1 - analyticalOmega1) / analyticalOmega1 * 100.0;

    std::cout << "First eigenfrequency: " << omega1 << " rad/s" << std::endl;
    std::cout << "Analytical value: " << analyticalOmega1 << " rad/s" << std::endl;
    std::cout << "Relative error: " << error1 << "%" << std::endl;
    std::cout << std::endl;

    // Save eigenmode to file
    std::ofstream out1("eigenmode1.csv");
    out1 << "x,u" << std::endl;
    for (int i = 0; i < n; ++i)
    {
        double x = (i + 1) * h;
        out1 << x << "," << eigenVector1(i) << std::endl;
    }
    out1.close();

    /*
     *   Aufgabe 2 - omega10
     *   Fazit for Alex: approximation is pretty scheiße due to n = 100 (too low like alexlow)
     */
    Vector<double> eigenVector10(n);
    double eigenValue10;

    ComputeInverseIteration(K, M, u, 5e9, eigenVector10, eigenValue10);

    eigenValue10 = std::abs(eigenValue10);
    double omega10 = std::sqrt(eigenValue10);
    double analyticalOmega10 = c * 10 * M_PI / L;
    double error10 = std::abs(omega10 - analyticalOmega10) / analyticalOmega10 * 100.0;

    std::cout << "Tenth eigenfrequency: " << omega10 << " rad/s" << std::endl;
    std::cout << "Analytical value: " << analyticalOmega10 << " rad/s" << std::endl;
    std::cout << "Relative error: " << error10 << "%" << std::endl;
    std::cout << std::endl;

    // Save eigenmode to file
    std::ofstream out2("eigenmode10.csv");
    out2 << "x,u" << std::endl;
    for (int i = 0; i < n; ++i)
    {
        double x = (i + 1) * h;
        out2 << x << "," << eigenVector10(i) << std::endl;
    }
    out2.close();

    /*
     *   Aufgabe 2 - omega5
     *   Here we ask ourself what the true meaning of life is
     */
    const int iterMax = 1e9;
    double targetError = 0.1;
    double errorIter = 100.0;
    int nIter = 100;
    while (errorIter > targetError && nIter <= iterMax)
    {
        std::cout << "Trying for n = " << nIter << "\t Please send help." << std::endl;
        int neIter = nIter + 1;    // number of elements
        double hIter = L / neIter; // mesh fineness

        // Create tridiagonal matrix
        TridiagSparseMatrix<double> KIter(nIter);
        TridiagSparseMatrix<double> MIter(nIter);

        // Fill the tridiagonal matrices according to equation (5) and (6)
        for (int i = 0; i < nIter; ++i)
        {
            KIter.Set(i, i, 2.0 * pow(c, 2) / hIter);
            MIter.Set(i, i, 4.0 * hIter / 6.0);

            if (i > 0)
            {
                KIter.Set(i, i - 1, -pow(c, 2) / hIter);
                MIter.Set(i, i - 1, hIter / 6.0);
            }

            if (i < nIter - 1)
            {
                KIter.Set(i, i + 1, -pow(c, 2) / hIter);
                MIter.Set(i, i + 1, hIter / 6.0);
            }
        }

        Vector<double> uIter(nIter);
        for (int i = 0; i < nIter; ++i)
        {
            double xi = (i + 1) * L / neIter;
            uIter(i) = sin(nIter * M_PI * xi / L);
        }

        Vector<double> eigenVectorIter(nIter);
        double eigenValueIter;
        double analyticalOmegaIter = c * 5 * M_PI / L;
        ComputeInverseIteration(KIter, MIter, uIter, pow(analyticalOmegaIter, 2) * 0.95, eigenVectorIter, eigenValueIter);

        eigenValueIter = std::abs(eigenValueIter);
        double omegaIter = std::sqrt(eigenValueIter);
        errorIter = std::abs(omegaIter - analyticalOmegaIter) / analyticalOmegaIter * 100.0;
        nIter *= 10;
    }
    std::cout << "Number of unknowns needed for 0.1\% accuracy for 5th eigenmode: " << nIter << std::endl;

    return 0;
}