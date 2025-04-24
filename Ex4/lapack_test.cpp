#include <iostream>

// LAPACK function prototypes
extern "C" {
    void dgesv_(int* n, int* nrhs, double* a, int* lda, int* ipiv, double* b, int* ldb, int* info);
}

int main() {
    int n = 3;
    int nrhs = 1;
    int lda = n;
    int ldb = nrhs;
    int ipiv[n];
    int info;

    // Coefficient matrix A
    double A[9] = {
        2, 1, -1,
        -3, -1, 2,
        -2, 1, 2
    };

    // Right-hand side vector b
    double b[3] = {8, -11, -3};

    // Solve the linear system A * x = b
    dgesv_(&n, &nrhs, A, &lda, ipiv, b, &ldb, &info);

    // Check for convergence
    if (info == 0) {
        std::cout << "Solution:" << std::endl;
        for (int i = 0; i < n; ++i) {
            std::cout << b[i] << std::endl;
        }
    } else {
        std::cerr << "LAPACK DGESV failed with info = " << info << std::endl;
    }

    return 0;
}
