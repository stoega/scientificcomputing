#include <iostream>
#include <fstream>

#include "SCTridiagSparseMatrix.h"
#include "SCTridiagLUSolver.h"

#include <Eigen/Eigenvalues>
#include <Eigen/Sparse>
#include <Eigen/SparseLU>

using namespace SC;
using namespace std;

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
    std::string fullFileName = "Ex4_A2_eigen_w" + std::to_string(modeNr) + ".csv";
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

/**
 * Computes the first m Eigenvalues and Eigenvectors using Accelerated Inverse Iteration using Rayleigh quotient
 * K * v = lambda * M * v
 * @param K Tridiagonal sparse matrix of the lhs
 * @param M Tridiagonal sparse matrix of the rhs
 * @param lambda Vector which stores/returns the m-Eigenvalues
 * @param V Vector of arrays which stores/returns the Eigenvector of the corresponding Eigenvalue
 */
void modifiedInverseIterationNRayleigh(TridiagSparseMatrix<double> &K,
                                       TridiagSparseMatrix<double> &M, Vector<double> &lambda, Vector<double> *V)
{
    // This function is based on "InverseIterationNRayleigh" from waveequation_Eigen.cpp

    // Convergence parameters
    const int MAX_ITERATIONS = 100;
    const double TOLERANCE = 1e-8;

    // dimensions
    int n = K.Height();
    int m = lambda.Size();

    // LU decomposition for K
    TridiagLUSolver<double> solver(K);

    // init V with random values and normalize
    for (int i = 0; i < m; i++)
    {
        V[i].SetSize(n);
        V[i].SetRandom(0, 1, i);
        V[i] *= 1. / V[i].Norm();
    }

    // init matrices & vectors
    // alternative Vector<Vector<double>> could be used instead of plain Arrays
    Vector<double> *MV = new Vector<double>[m];
    Vector<double> *KV = new Vector<double>[m];
    Vector<double> *W = new Vector<double>[m];
    Vector<double> *KW = new Vector<double>[m];
    Vector<double> *MW = new Vector<double>[m];
    Vector<double> *Vnew = new Vector<double>[m];

    for (int i = 0; i < m; i++)
    {
        MV[i].SetSize(n);
        KV[i].SetSize(n);
        W[i].SetSize(n);
        MW[i].SetSize(n);
        KW[i].SetSize(n);
        Vnew[i].SetSize(n);
        Vnew[i].SetAll(0.);
    }

    // Error tracking
    double error_init = m;
    Vector<double> error(m);
    error.SetAll(1.);

    int step = 0;

    // Reduced 2m x 2m Matrices
    Eigen::MatrixXd Ms(2 * m, 2 * m), Ks(2 * m, 2 * m);
    Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXd> es;
    while (step < MAX_ITERATIONS && error.Norm() > TOLERANCE * error_init)
    {
        step++;

        // Compute MVi and KVi for all basis vectors
        for (int i = 0; i < m; i++)
        {
            // MV = M*V, KV = K*V
            M.Apply(V[i], MV[i]);
            K.Apply(V[i], KV[i]);

            // next update direction : w_i = v_i - lambda_i K ^ -1 M v_i
            solver.Apply(MV[i], W[i]);
            W[i] *= -lambda(i);
            W[i] += V[i];

            // use residual w as error indicator
            error(i) = W[i].Norm();

            // KW = K*W, MW = M*W for subspace projection
            M.Apply(W[i], MW[i]);
            K.Apply(W[i], KW[i]);
        }

        // Initialize error norm on first iteration
        if (step == 1)
            error_init = error.Norm();

        // Setup reduced 2m x 2m eigenvalue Problem
        for (int i = 0; i < m; i++)
            for (int j = 0; j < m; j++)
            {
                Ks(i, j) = InnerProduct(V[i], KV[j]);
                Ks(i, m + j) = InnerProduct(W[j], KV[i]);
                Ks(m + j, i) = Ks(i, m + j);
                Ks(m + i, m + j) = InnerProduct(W[i], KW[j]);

                Ms(i, j) = InnerProduct(V[i], MV[j]);
                Ms(i, m + j) = InnerProduct(W[j], MV[i]);
                Ms(m + j, i) = Ms(i, m + j);
                Ms(m + i, m + j) = InnerProduct(W[i], MW[j]);
            }

        // Solve the reduced eigenvalue problem using EIGEN
        es.compute(Ks, Ms);
        const auto &lam_small = es.eigenvalues();
        const auto &y_small = es.eigenvectors();

        // compute Vnew based on y_small and V & W
        for (int s = 0; s < m; s++) // s = spalte
        {
            // Reset new vector
            Vnew[s].SetAll(0.0);

            for (int z = 0; z < n; z++) // z = zeile
            {
                for (int i = 0; i < m; i++)
                {
                    Vnew[s](z) += V[i](z) * y_small(i, s) + W[i](z) * y_small(i + m, s);
                }
            }
        }

        // normalize Vnew and update V, store new Eigenvalue
        for (int i = 0; i < m; i++)
        {
            Vnew[i] *= (1. / Vnew[i].Norm());
            V[i] = Vnew[i];
            lambda(i) = lam_small(i);
        }
    }
}

int main()
{
    // Problem parameters
    const int m = 5;
    const float L = .5;   // length of flute
    const int c = 343;    // speed of sound
    const int n = 100;    // number of unknowns
    const int ne = n + 1; // number of elements

    const double tol = 1e-9;  // error tolerance
    const int maxIter = 1000; // max iterations

    const double h = L / ne; // mesh fineness

    TridiagSparseMatrix<double> K(n);
    FillTridiagSparseMatrix(K, 2.0 * pow(c, 2) / h, -pow(c, 2) / h);

    TridiagSparseMatrix<double> M(n);
    FillTridiagSparseMatrix(M, 4.0 * h / 6.0, h / 6.0);

    SC::Vector<double> V[m];
    SC::Vector<double> lambda(m);
    lambda.SetAll(1.);

    modifiedInverseIterationNRayleigh(K, M, lambda, V);

    for (int idx = 0; idx < m; idx++)
    {
        std::cout << "Eigenwert lambda_" << idx + 1 << " = " << lambda(idx) << std::endl;
        std::cout << "Eigenfrequenz w_" << idx + 1 << " = " << sqrt(lambda(idx)) << std::endl;
        WriteModeToCSV(idx + 1, V[idx], h);
        // V[idx].Print(std::cout);
        // std::cout << std::endl;
    }

    return 0;
}