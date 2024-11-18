#include <random>
#include <iostream>
#include <iomanip>

#include <SCvector.h>
#include <SCmatrix.h>

#include <SCGaussSolver.h>

using namespace SC;
using namespace std;

int main()
{
    
    int size = 10;

    Vector<> b(size);
    b.SetAll(3.);


    Matrix<> A(size,size);
    A.SetAll(0.);
    for (int i=0; i<size; i++)
    {
        A(i,i) = 2.;
        if (i > 0) A(i-1,i) = -1;
        if (i > 0) A(i,i-1) = -1;
    }

    Matrix<> A_random(A);
    std::uniform_real_distribution<double> unif(0.0, 100);
    std::default_random_engine re;
    re.seed(0);
    double factor = 1.;
    A_random.SetAll(0);
    for (int i=0; i<size; i++)
    {
        for (int j=0; j<size; j++)
        {
            A_random(i,j) += unif(re);
        }
    }
    // A_random(0,0) = 1e-5;

    

    SimpleGaussSolver<> invA_simple(A);
    SimpleGaussSolver<> invA_random_simple(A_random);
    GaussSolver<> invA(A);
    GaussSolver<> invA_random(A_random);

    Vector<> Ab;
    A.Apply(b,Ab);
    Vector<> A_random_b;
    A_random.Apply(b,A_random_b);

    Vector <> res, res_piv;
    Vector <> res_random, res_random_piv;
    invA_simple.Apply(Ab, res);
    invA_random_simple.Apply(A_random_b, res_random);

    invA.Apply(Ab, res_piv);
    invA_random.Apply(A_random_b, res_random_piv);

    if (size <= 10)
    {
        cout << "---- SPD matrix A ----" << endl;
        cout << "-- Gauss solver (no pivoting) --" << endl;
        cout << setprecision(15) << "res = " << res << endl;
        cout << "-- Gauss inverse with pivoting --" << endl;
        cout << "res = " << res_piv << endl;


        cout << "---- random matrix A ----" << endl;
        cout << "-- Gauss solver (no pivoting) --" << endl;
        cout << "res = " << res_random << endl;
        cout << "-- Gauss inverse with pivoting --" << endl;
        cout << "res = " << res_random_piv << endl;
    }
    Vector<> error(res), error_piv(res_piv), error_random_piv(res_random_piv), error_random(res_random);
    error -= b;
    error_random -= b;
    error_random_piv -= b;
    error_piv -= b;

    cout << "problem size " << size << endl;
    cout << "relative error: SPD matrix, simple Gauss:    " << error.Norm()/b.Norm() << endl;
    cout << "relative error: SPD matrix, pivoting:        " << error_piv.Norm()/b.Norm() << endl;
    cout << "relative error: random matrix, simple Gauss: " << error_random.Norm()/b.Norm() << endl;
    cout << "relative error: random matrix, pivoting:     " << error_random_piv.Norm()/b.Norm() << endl;


}