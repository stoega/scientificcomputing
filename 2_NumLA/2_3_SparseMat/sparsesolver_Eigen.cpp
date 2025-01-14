#define _USE_MATH_DEFINES
#include <math.h>

#include <Eigen/Sparse>
#include <Eigen/SparseLU>
#include <SCvector.h>
#include <SCmatrix.h>
#include <iostream>
#include <vector>
#include <chrono>
#include <fstream>

using namespace std;

int main()
{
    int n_x = 100;
    int n_vec = n_x*(n_x);
    double meshsize = 1./(n_x+1);
    double k_coeff = 1.;
    double q_coeff = 1.;

    double meshsize2 = meshsize*meshsize;

    Eigen::SparseMatrix<double, Eigen::RowMajor> Amat(n_vec,n_vec);

    std::vector<Eigen::Triplet<double> > tripletList;
    tripletList.reserve(5*n_vec);
    for(int i=0; i<n_x; i++)
    {
        for (int j=0; j<n_x; j++)
        {
            tripletList.push_back(Eigen::Triplet<double>(i+j*n_x, i+j*n_x, 4.*k_coeff/meshsize2));
            if (j>0) tripletList.push_back(Eigen::Triplet<double>(i+(j-1)*n_x, i+j*n_x, -k_coeff/meshsize2));
            if (j<n_x-1) tripletList.push_back(Eigen::Triplet<double>(i+(j+1)*n_x, i+j*n_x, -k_coeff/meshsize2));
            if (i<n_x-1) tripletList.push_back(Eigen::Triplet<double>(i+1+(j)*n_x, i+j*n_x, -k_coeff/meshsize2));
            if (i>0) tripletList.push_back(Eigen::Triplet<double>(i-1+(j)*n_x, i+j*n_x, -k_coeff/meshsize2));
        }
    }

    Amat.setFromTriplets(tripletList.begin(), tripletList.end());

    Eigen::Matrix<double, Eigen::Dynamic, 1> b(n_vec), u(n_vec);
    b.setConstant(q_coeff);
    cout << "system size " << n_vec << endl;
    cout << "setup completed, factorization..." << endl;
    auto begin_sparse = std::chrono::high_resolution_clock::now();

    Eigen::SparseLU<Eigen::SparseMatrix<double, Eigen::RowMajor>> Ainv;
    Ainv.analyzePattern(Amat);
    Ainv.factorize(Amat);

    u = Ainv.solve(b);
    auto end_sparse = std::chrono::high_resolution_clock::now();
    auto elapsed_sparse = std::chrono::duration_cast<std::chrono::nanoseconds>(end_sparse - begin_sparse);

    // std::cout << "solution computed via Eigen Sparse LU" << std::endl;
    // for (int i=0; i<n_x; i++)
    // {
    //     for (int j=0; j<n_x; j++)
    //         std::cout<< u[i*n_x+j] << " ";
    //     std::cout << "\n";
    // }

    cout << "Timing: using Eigen SparseLU: " << elapsed_sparse.count()/1e9 << "s" << endl;

    std::cout << "export solution to file" << std::endl;
    ofstream outfile("poisson_solution_eigen.txt");

    for (int j=0; j<n_x+2; j++)
        outfile << "0 ";
    outfile << "\n";
    for (int i=0; i<n_x; i++)
    {
        outfile << "0 ";
        for (int j=0; j<n_x; j++)
            outfile << u[i*n_x+j] << " ";
        outfile << "0\n";
    }
    for (int j=0; j<n_x+2; j++)
        outfile << "0 ";
    outfile << "\n";
}
