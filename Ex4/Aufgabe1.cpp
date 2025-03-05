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
/**
 * @param v Vector to be normalized in place
 */
void normalize(SC::Vector<T> &v)
{
    T norm = v.Norm();
    v.Mult(1 / norm);
}

/**
 * Compute the eigenvector and value using the inverse iteration with shift.
 *
 * @param[in] K Tridiagonal sparse matrix
 * @param[in] M Tridiagonal sparse matrix
 * @param[in] u Starting value of v
 * @param[in] mu Shift of eigenvalue
 * @param[out] eigenVec Approximated eigenvector
 * @param[out] eigenValue Approximated eigenvalue
 * @param[in] tol (Optional) Specifies the tolerance for the residual criteria
 * @param[in] maxIter (Optional) Specifies the maximum allowed iterations
 */
void ComputeInverseIteration(const TridiagSparseMatrix<double> &K,
                             const TridiagSparseMatrix<double> &M,
                             const Vector<double> &u,
                             double mu,
                             Vector<double> &eigenVec,
                             double &eigenValue,
                             double tol = 1e-9,
                             int maxIter = 1000)
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

#ifndef NDEBUG
    std::cout << "Starting inverse iteration with shift = " << mu << std::endl;
#endif
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
#ifndef NDEBUG
            std::cout << "Converged after " << iter + 1 << " iterations." << std::endl;
#endif
            break;
        }

        // (3.32) right
        v = z;

#ifndef NDEBUG
        // Print every 10 iterations
        if (iter % 10 == 0)
        {
            std::cout << "Iteration " << iter
                      << ", eigenvalue = " << lambda + mu
                      << ", residual = " << r.Norm() << std::endl;
        }
#endif
    }
    // lamda = eigenvalue - shift
    eigenValue = lambda + mu;
    eigenVec = z;
}

template <typename T>
/**
 * Fills a matrix of type SC::TridiagSparseMatrix with proviede diagonal and offdiagonal values
 * @param[out] A Matrix to be filled
 * @param[in] diagonalValue Values on the diagonal
 * @param[in] offdiagonalValue Values of the offdiagonal
 */
void FillTridiagSparseMatrix(TridiagSparseMatrix<T> &A, T diagonalValue, T offdiagonalValue)
{
    size_t n = A.Height();
    for (int i = 0; i < n; ++i)
    {
        // Diagonal elements
        A.Set(i, i, diagonalValue);

        // Subdiagonal elements (left)
        if (i > 0)
        {
            A.Set(i, i - 1, offdiagonalValue);
        }

        // Superdiagonal elements (right)
        if (i < n - 1)
        {
            A.Set(i, i + 1, offdiagonalValue);
        }
    }
}

/**
* Writes a vector to a .csv file
*
* @param modeNr Number of eigenmode
* @param vec The vector to write
* @param h Mesh fineness of vector
*/
void WriteModeToCSV(int modeNr, Vector<double> &vec, double h)
{
    std::string fullFileName = "Ex4_A1_w" + std::to_string(modeNr) + ".csv";
    size_t n = vec.Size();
    std::ofstream out(fullFileName);
    out << "x,u" << std::endl;
    for (int i = 0; i < n; ++i)
    {
        double x = (i + 1) * h;
        out << x << "," << vec(i) << std::endl;
    }
    out.close();
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

    // Values of eigenproblem
    Vector<double> eigenVector(n);
    double eigenValue = 0.0;
    double omega = 0.0;
    double analyticalOmega = 0.0;
    double error = 0.0;

    // Create tridiagonal matrix
    // Create and fill the tridiagonal matrices according to equation (5) and (6)
    TridiagSparseMatrix<double> K(n);
    K.Print(std::cout);
    FillTridiagSparseMatrix(K, 2.0 * pow(c, 2) / h, -pow(c, 2) / h);
    K.Print(std::cout);

    TridiagSparseMatrix<double> M(n);
    FillTridiagSparseMatrix(M, 4.0 * h / 6.0, h / 6.0);

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
    std::cout << "----- Aufgabe 1: omega_1, n = 100 -----" << std::endl;

    ComputeInverseIteration(K, M, u, 0, eigenVector, eigenValue);

    omega = std::sqrt(eigenValue);
    analyticalOmega = c * M_PI / L;
    error = std::abs(omega - analyticalOmega) / analyticalOmega * 100.0;

    std::cout << "Target eigenfrequency: " << omega << " rad/s" << std::endl;
    std::cout << "Analytical value: " << analyticalOmega << " rad/s" << std::endl;
    std::cout << "Relative error: " << error << "%" << std::endl;
    std::cout << std::endl;

    WriteModeToCSV(1, eigenVector, h);

    /*
     *   Aufgabe 2 - omega10
     *   Fazit for Alex: approximation is pretty scheiße due to n = 100 (too low like alexlow)
     */
    std::cout << "----- Aufgabe 2: omega_10, n = 100 -----" << std::endl;

    analyticalOmega = c * 10 * M_PI / L;

    ComputeInverseIteration(K, M, u, analyticalOmega * analyticalOmega * 0.99, eigenVector, eigenValue);

    omega = std::sqrt(eigenValue);
    error = std::abs(omega - analyticalOmega) / analyticalOmega * 100.0;

    std::cout << "Target eigenfrequency: " << omega << " rad/s" << std::endl;
    std::cout << "Analytical value: " << analyticalOmega << " rad/s" << std::endl;
    std::cout << "Relative error: " << error << "%" << std::endl;
    std::cout << std::endl;

    WriteModeToCSV(10, eigenVector, h);

    /*
     *   Aufgabe 2 - omega5
     *   Here we ask ourself what the true meaning of life is
     */
    // Iteration params
    const int iterMax = 1e6;
    double targetError = 0.1;
    int nIter = 10;
    int iterStep = 10;

    // reset error
    error = 100.0;

    std::cout << "----- Aufgabe 2: n = ? -----" << std::endl;
    while (nIter <= iterMax)
    {
        int neIter = nIter + 1;    // number of elements
        double hIter = L / neIter; // mesh fineness

        // Create and fill the tridiagonal matrix
        TridiagSparseMatrix<double> KIter(nIter);
        FillTridiagSparseMatrix(KIter, 2.0 * pow(c, 2) / hIter, -pow(c, 2) / hIter);
        TridiagSparseMatrix<double> MIter(nIter);
        FillTridiagSparseMatrix(MIter, 4.0 * hIter / 6.0, hIter / 6.0);

        Vector<double> uIter(nIter);
        for (int i = 0; i < nIter; ++i)
        {
            double xi = (i + 1) * L / neIter;
            uIter(i) = sin(nIter * M_PI * xi / L);
        }

        // Vector<double> eigenVectorIter(nIter);
        // double eigenValueIter;
        eigenVector = Vector<double>(nIter);

        analyticalOmega = c * 5 * M_PI / L;

        ComputeInverseIteration(KIter, MIter, uIter, analyticalOmega * analyticalOmega * 0.99, eigenVector, eigenValue);
        omega = std::sqrt(eigenValue);
        error = std::abs(omega - analyticalOmega) / analyticalOmega * 100.0;

        std::cout << "n = " << nIter << ": \n\tanalytical omega = " << analyticalOmega << "\n\tapproximated omega = " << omega << "\n\terror = " << error << "%\n"
                << std::endl;

        
        if(error <= targetError){
            break;
        }

        nIter += iterStep;
    }
    std::cout << "Number of unknowns needed for < 0.1\% rel. error for 5th eigenmode: " << nIter << std::endl;
    WriteModeToCSV(5, eigenVector, L / (nIter + 1));

    return 0;
}