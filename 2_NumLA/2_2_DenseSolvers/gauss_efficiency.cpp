#include <chrono>

#include <SCvector.h>
#include <SCmatrix.h>

#include <SCGaussSolver.h>

using namespace SC;
using namespace std;

int main()
{
    int size = 1000;

    int n_solve = 1;
    cout << "system size " << size << endl;
    cout << "number of repeated solves " << n_solve << endl;

    Vector<> b(size);
    b.SetAll(3.);

    Matrix<> A(size,size);
    A.SetAll(0.);
    for (int i=0; i<size; i++)
    {
        A(i,i) = 2.;
        if (i > 0) A(i-1,i) = -0.5;
        if (i > 0) A(i,i-1) = -0.5;
    }

    Vector <> c, d, c2, d2;
    
    // Test 1 - simple Gauss solver, repeated solve
    auto begin_1 = std::chrono::high_resolution_clock::now();

    SimpleGaussSolver<> invA_simple(A);

    for (int i=0; i<n_solve; i++)
        invA_simple.Apply(b, c);

    auto end_1 = std::chrono::high_resolution_clock::now();
    auto elapsed_1 = std::chrono::duration_cast<std::chrono::nanoseconds>(end_1 - begin_1);
    cout << "Using simple Gauss solver, sec                 " << elapsed_1.count()/1e9 << endl;


    // Test 1 - simple Gauss solver, repeated solve
    auto begin_2 = std::chrono::high_resolution_clock::now();

    GaussSolver<> invA(A);

    for (int i=0; i<n_solve; i++)
        invA.Apply(b, c2);

    auto end_2 = std::chrono::high_resolution_clock::now();
    auto elapsed_2 = std::chrono::duration_cast<std::chrono::nanoseconds>(end_2 - begin_2);
    cout << "Using Gauss inverse, full pivoting solver, sec " << elapsed_2.count()/1e9 << endl;



}