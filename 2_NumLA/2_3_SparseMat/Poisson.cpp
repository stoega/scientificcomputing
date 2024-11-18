#include <vector>
#include <chrono>
#include <iostream>
#include <fstream>

#include <SCvector.h>
#include <SCmatrix.h>
#include <SCSparseMatrix.h>
#include <SCLUSolver.h>


using namespace SC;
using namespace std;

// implements method of steepest descent - slow but simple
template <typename T>
int GradientSolver(const LinearOperator<T>& K, const Vector<T>& f, Vector<T>& u, double damping, int maxsteps=1000)
{
    u.SetAll(0.);

    // compute the residual res = f - K u
    Vector<T> res;
    K.Apply(u,res,-1);
    res += f;
    double normres_init = res.Norm();
    double normres;

    int step;
    for (step = 0; step < maxsteps; step++)
    {
        u.AddMultiple(damping, res);

        K.Apply(u,res,-1);
        res += f;

        normres = res.Norm();
        if (normres < normres_init*1e-6)
            break;

        if ((step)%100==0) cout << "step " << step << ": error " << normres << endl;

    }
    cout << "Gradient solver needed " << step << " steps, residual error " << normres << endl;
    return normres;

}

int main()
{
    int n_x = 50;
    int n = n_x*n_x;

    double meshsize = 1./(n_x+1);
    double k_coeff = 1.;
    double q_coeff = 1.;

    double meshsize2 = meshsize*meshsize;

    std::vector<int> row, col;
    std::vector<double> val;

    Matrix<> Adense(n,n);
    Adense.SetAll(0.);


    int k = 0;
    for (int i=0; i<n_x; i++)
    {
        for (int j=0; j<n_x; j++)
        {
            // equation i*n_x + j

            row.push_back(i*n_x+j);
            col.push_back(i*n_x+j);
            val.push_back(4.*k_coeff/meshsize2);
            Adense(i*n_x+j, i*n_x+j) = 4.*k_coeff/meshsize2;

            if (i>0)
            {
                row.push_back(i*n_x+j);
                col.push_back((i-1)*n_x+j);
                val.push_back(-1.*k_coeff/meshsize2);
                Adense(i*n_x+j, (i-1)*n_x+j) = -1.*k_coeff/meshsize2;
            }
            if (i<n_x-1)
            {
                row.push_back(i*n_x+j);
                col.push_back((i+1)*n_x+j);
                val.push_back(-1.*k_coeff/meshsize2);
                Adense(i*n_x+j, (i+1)*n_x+j) = -1.*k_coeff/meshsize2;
            }
            if (j>0)
            {
                row.push_back(i*n_x+j);
                col.push_back((i)*n_x+j-1);
                val.push_back(-1.*k_coeff/meshsize2);
                Adense(i*n_x+j, (i)*n_x+j-1) = -1.*k_coeff/meshsize2;
            }
            if (j<n_x-1)
            {
                row.push_back(i*n_x+j);
                col.push_back((i)*n_x+j+1);
                val.push_back(-1.*k_coeff/meshsize2);
                Adense(i*n_x+j, (i)*n_x+j+1) = -1.*k_coeff/meshsize2;
            }

        }
    }
    
    SparseMatrix<double> A(n, n, row, col, val);

    // std::cout << "dense matrix " << Adense << std::endl;
    // std::cout << "sparse matrix" << A << std::endl;

    Vector<double> b(n);
    b.SetAll(q_coeff);

    // solve dense problem via LU factorization
    auto begin_dense = std::chrono::high_resolution_clock::now();
    
    LUSolver<> invA(Adense);

    Vector<double> u;
    invA.Apply(b, u);

    auto end_dense = std::chrono::high_resolution_clock::now();
    auto elapsed_dense = std::chrono::duration_cast<std::chrono::nanoseconds>(end_dense - begin_dense);

    // std::cout << "solution " << std::endl;
    // for (int i=0; i<n_x; i++)
    // {
    //     for (int j=0; j<n_x; j++)
    //         std::cout<< u[i*n_x+j] << " ";
    //     std::cout << "\n";
    // }

    // cout << "pattern of LU factorization " << endl;
    // invA.PrintPattern(cout);
    // cout << endl << endl;

    // solve sparse problem via iterative GradientSolver
    auto begin_sparse = std::chrono::high_resolution_clock::now();

    Vector<double> u2(n);
    GradientSolver(A, b, u2, 1./A.Get(0,0), 10000);
    auto end_sparse = std::chrono::high_resolution_clock::now();
    auto elapsed_sparse = std::chrono::duration_cast<std::chrono::nanoseconds>(end_sparse - begin_sparse);

    // std::cout << "solution computed iteratively" << std::endl;
    // for (int i=0; i<n_x; i++)
    // {
    //     for (int j=0; j<n_x; j++)
    //         std::cout<< u2[i*n_x+j] << " ";
    //     std::cout << "\n";
    // }

    cout << "Timing: using dense matrix: " << elapsed_dense.count()/1e9 << "s" << endl;
    cout << "Timing: using sparse matrix: " << elapsed_sparse.count()/1e9 << "s" << endl;


    std::cout << "export solution to file" << std::endl;
    ofstream outfile("poisson_solution.txt");

    for (int j=0; j<n_x+2; j++)
        outfile << "0 ";
    outfile << "\n";
    for (int i=0; i<n_x; i++)
    {
        outfile << "0 ";
        for (int j=0; j<n_x; j++)
            outfile << u2(i*n_x+j) << " ";
        outfile << "0\n";
    }
    for (int j=0; j<n_x+2; j++)
        outfile << "0 ";
    outfile << "\n";
}
