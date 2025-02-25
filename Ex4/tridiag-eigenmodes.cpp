#include <iostream>
#include <cmath>
#include <vector>
#include <fstream>
#include "SCTridiagSparseMatrix.h"

using namespace SC;

// Normalize a vector (divide by its Euclidean norm)
template <typename T>
void normalize(std::vector<T>& v) {
    T norm = 0.0;
    for (const auto& val : v) {
        norm += val * val;
    }
    norm = std::sqrt(norm);
    
    for (auto& val : v) {
        val /= norm;
    }
}

// Inverse iteration with shift to find eigenvalue and eigenvector
template <typename T>
void inverseIteration(const TridiagSparseMatrix<T>& K, const TridiagSparseMatrix<T>& M, 
                     T shift, std::vector<T>& eigenVector, T& eigenValue,
                     int maxIter = 100, T tol = 1e-10) {
    
    int n = K.getSize();
    std::vector<T> v(n, 0.0);
    std::vector<T> w(n, 0.0);
    
    // Initial guess for eigenvector (sine function)
    for (int i = 0; i < n; ++i) {
        v[i] = std::sin(M_PI * (i + 1) / (n + 1));
    }
    normalize(v);
    
    // Create K - shift*M matrix
    TridiagSparseMatrix<T> A = K;
    for (int i = 0; i < n; ++i) {
        A(i, i) -= shift * M(i, i);
        if (i > 0) {
            A(i, i-1) -= shift * M(i, i-1);
        }
        if (i < n - 1) {
            A(i, i+1) -= shift * M(i, i+1);
        }
    }
    
    // LU decomposition of A
    TridiagSparseMatrix<T> L(n), U(n);
    A.LUdecomposition(L, U);
    
    T lambdaOld = 0.0;
    T lambdaNew = 0.0;
    
    std::cout << "Starting inverse iteration with shift = " << shift << std::endl;
    
    // Inverse iteration
    for (int iter = 0; iter < maxIter; ++iter) {
        // Solve Aw = v
        A.LUSolve(L, U, v, w);
        
        // Rayleigh quotient
        T numerator = 0.0;
        T denominator = 0.0;
        
        for (int i = 0; i < n; ++i) {
            T Kw_i = 0.0;
            if (i > 0) Kw_i += K(i, i-1) * w[i-1];
            Kw_i += K(i, i) * w[i];
            if (i < n - 1) Kw_i += K(i, i+1) * w[i+1];
            
            T Mw_i = 0.0;
            if (i > 0) Mw_i += M(i, i-1) * w[i-1];
            Mw_i += M(i, i) * w[i];
            if (i < n - 1) Mw_i += M(i, i+1) * w[i+1];
            
            numerator += w[i] * Kw_i;
            denominator += w[i] * Mw_i;
        }
        
        lambdaNew = numerator / denominator;
        
        // Normalize w
        normalize(w);
        
        // Check convergence
        T diff = std::abs(lambdaNew - lambdaOld);
        if (diff < tol) {
            std::cout << "Converged after " << iter + 1 << " iterations." << std::endl;
            break;
        }
        
        lambdaOld = lambdaNew;
        v = w;
        
        if (iter % 10 == 0) {
            std::cout << "Iteration " << iter << ", eigenvalue = " << lambdaNew 
                      << ", diff = " << diff << std::endl;
        }
    }
    
    eigenVector = w;
    eigenValue = lambdaNew;
}

int main() {
    // Parameters
    const double L = 0.5;  // Pipe length in meters
    const double c = 343.0;  // Speed of sound in m/s
    const int n = 100;  // Number of unknowns (grid points)
    const int ne = n + 1;  // Number of elements
    const double h = L / ne;  // Grid spacing
    
    // Create tridiagonal matrices K and M
    TridiagSparseMatrix<double> K(n);
    TridiagSparseMatrix<double> M(n);
    
    // Fill K and M according to the problem description
    for (int i = 0; i < n; ++i) {
        K(i, i) = 2.0 * c * c / h;
        M(i, i) = 4.0 * h / 6.0;
        
        if (i > 0) {
            K(i, i-1) = -c * c / h;
            M(i, i-1) = h / 6.0;
        }
        
        if (i < n - 1) {
            K(i, i+1) = -c * c / h;
            M(i, i+1) = h / 6.0;
        }
    }
    
    // Part 1: Compute the first eigenmode without shift
    std::vector<double> eigenVector1(n);
    double eigenValue1;
    
    inverseIteration(K, M, 0.0, eigenVector1, eigenValue1, 1000, 1e-10);
    
    // Convert eigenvalue to frequency
    double omega1 = std::sqrt(eigenValue1);
    double analyticalOmega1 = c * M_PI / L;
    double error1 = std::abs(omega1 - analyticalOmega1) / analyticalOmega1 * 100.0;
    
    std::cout << "First eigenfrequency: " << omega1 << " rad/s" << std::endl;
    std::cout << "Analytical value: " << analyticalOmega1 << " rad/s" << std::endl;
    std::cout << "Relative error: " << error1 << "%" << std::endl;
    
    // Save eigenmode to file
    std::ofstream out1("eigenmode1.csv");
    out1 << "x,u" << std::endl;
    for (int i = 0; i < n; ++i) {
        double x = (i + 1) * h;
        out1 << x << "," << eigenVector1[i] << std::endl;
    }
    out1.close();
    
    // Part 2: Compute the tenth eigenmode with shift
    std::vector<double> eigenVector10(n);
    double eigenValue10;
    
    // Choose a shift close to the expected eigenvalue
    double expectedOmega10 = c * 10 * M_PI / L;
    double shift = expectedOmega10 * expectedOmega10 * 0.99; // Slightly less than expected value
    
    inverseIteration(K, M, shift, eigenVector10, eigenValue10, 1000, 1e-10);
    
    // Convert eigenvalue to frequency
    double omega10 = std::sqrt(eigenValue10);
    double analyticalOmega10 = c * 10 * M_PI / L;
    double error10 = std::abs(omega10 - analyticalOmega10) / analyticalOmega10 * 100.0;
    
    std::cout << "Tenth eigenfrequency: " << omega10 << " rad/s" << std::endl;
    std::cout << "Analytical value: " << analyticalOmega10 << " rad/s" << std::endl;
    std::cout << "Relative error: " << error10 << "%" << std::endl;
    
    // Save eigenmode to file
    std::ofstream out10("eigenmode10.csv");
    out10 << "x,u" << std::endl;
    for (int i = 0; i < n; ++i) {
        double x = (i + 1) * h;
        out10 << x << "," << eigenVector10[i] << std::endl;
    }
    out10.close();
    
    // Find number of unknowns needed for 0.1% accuracy for the 5th eigenmode
    int testN = 10;
    double targetError = 0.1; // 0.1%
    double error5 = 100.0;
    
    while (error5 > targetError && testN <= 1000) {
        int testNe = testN + 1;
        double testH = L / testNe;
        
        TridiagSparseMatrix<double> testK(testN);
        TridiagSparseMatrix<double> testM(testN);
        
        for (int i = 0; i < testN; ++i) {
            testK(i, i) = 2.0 * c * c / testH;
            testM(i, i) = 4.0 * testH / 6.0;
            
            if (i > 0) {
                testK(i, i-1) = -c * c / testH;
                testM(i, i-1) = testH / 6.0;
            }
            
            if (i < testN - 1) {
                testK(i, i+1) = -c * c / testH;
                testM(i, i+1) = testH / 6.0;
            }
        }
        
        std::vector<double> testEigenVector(testN);
        double testEigenValue;
        
        double expectedOmega5 = c * 5 * M_PI / L;
        double testShift = expectedOmega5 * expectedOmega5 * 0.99;
        
        try {
            inverseIteration(testK, testM, testShift, testEigenVector, testEigenValue, 1000, 1e-8);
            
            double testOmega5 = std::sqrt(testEigenValue);
            double analyticalOmega5 = c * 5 * M_PI / L;
            error5 = std::abs(testOmega5 - analyticalOmega5) / analyticalOmega5 * 100.0;
            
            std::cout << "With n = " << testN << ", 5th eigenfrequency error: " << error5 << "%" << std::endl;
        }
        catch (...) {
            std::cout << "Calculation failed for n = " << testN << std::endl;
        }
        
        testN *= 2; // Increase grid size
    }
    
    std::cout << "Number of unknowns needed for 0.1% accuracy for 5th eigenmode: " << testN / 2 << std::endl;
    
    return 0;
}
