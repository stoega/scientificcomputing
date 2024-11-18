#include <vector>
#include <iostream>
#include <fstream>

#include <SCvector.h>
#include <SCmatrix.h>
#include <SCSparseMatrix.h>

#include <SCCGSolver.h>
#include <SCGMRESSolver.h>

using namespace SC;
using namespace std;

int main()
{
    // solve
    // -div( k nabla u) + (bx,0).nabla u = q

    // use CG and GMRES
    // CG works only for symmetric positive definite -> bx = 0!

    int n_x = 100;
    int n = n_x*n_x;

    double meshsize = 1./(n_x+1);
    double k_coeff = 1.;
    double q_coeff = 1.;
    double bx_coeff = 20;

    double meshsize2 = meshsize*meshsize;

    std::vector<int> row, col;
    std::vector<double> val;


    int k = 0;
    for (int i=0; i<n_x; i++)
    {
        for (int j=0; j<n_x; j++)
        {
            // equation i*n_x + j

            row.push_back(i*n_x+j);
            col.push_back(i*n_x+j);
            val.push_back(4.*k_coeff/meshsize2 + bx_coeff/meshsize);

            if (i>0)
            {
                row.push_back(i*n_x+j);
                col.push_back((i-1)*n_x+j);
                val.push_back(-1.*k_coeff/meshsize2);
            }
            if (i<n_x-1)
            {
                row.push_back(i*n_x+j);
                col.push_back((i+1)*n_x+j);
                val.push_back(-1.*k_coeff/meshsize2);
            }
            if (j>0)
            {
                row.push_back(i*n_x+j);
                col.push_back((i)*n_x+j-1);
                val.push_back(-1.*k_coeff/meshsize2 - bx_coeff/meshsize);
            }
            if (j<n_x-1)
            {
                row.push_back(i*n_x+j);
                col.push_back((i)*n_x+j+1);
                val.push_back(-1.*k_coeff/meshsize2);
            }

        }
    }

    
    SparseMatrix<double> A(n, n, row, col, val);

    Vector<double> b(n);
    b.SetAll(q_coeff);

    auto begin = std::chrono::high_resolution_clock::now();

    // CGSolver<> invA(A);
    GMRESSolver<> invA(A, 1000,20);

    Vector<double> u;
    invA.Apply(b, u);
    auto end = std::chrono::high_resolution_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);

    cout << "Timing: " << elapsed.count()/1e9 << "s" << endl;


    std::cout << "export solution to file" << std::endl;
    ofstream outfile("poisson_solution.txt");

    for (int j=0; j<n_x+2; j++)
        outfile << "0 ";
    outfile << "\n";
    for (int i=0; i<n_x; i++)
    {
        outfile << "0 ";
        for (int j=0; j<n_x; j++)
            outfile << u(i*n_x+j) << " ";
        outfile << "0\n";
    }
    for (int j=0; j<n_x+2; j++)
        outfile << "0 ";
    outfile << "\n";
}