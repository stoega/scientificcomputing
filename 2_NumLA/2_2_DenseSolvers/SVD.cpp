#include <random>
#include <fstream>
#include <iostream>

#include <SCvector.h>
#include <SCmatrix.h>

#define _USE_MATH_DEFINES
#include <math.h>

#include <Eigen/Dense>

#include <SCSVDSolver.h>

using namespace SC;
using namespace std;

typedef complex<double> COMPLEX;    // write COMPLEX instead of complex<double>
enum { DS = Eigen::Dynamic }; // Eigen-lib dynamic size identifier

void squareSVD(){
        // Problem 1 - solve problem with square, nonsymm matrix
    int size = 8;

    Vector<> x1(size);
    x1.SetAll(3.);

    Matrix<> A(size, size, TMatrixStorage::COL_MAJOR);
    A.SetRandom(-1,1);

    Vector<> b;
    A.Apply(x1,b);


    SVDSolver invA(A);
    Vector <> x;
    invA.Apply(b, x);

    cout << "SVD Solver, square matrix:" << endl;
    cout << "invA * b = " << x << endl;


}


void squareSVD_singular(){
    // Problem 2 - solve problem with square, nonsymm, singular matrix
    int size = 8;

    Matrix<> A(size, size, TMatrixStorage::COL_MAJOR);
    A.SetRandom(-1,1);

    // two identical lines - one singular value
    for (int i=0; i<size; i++)
        A(0,i) = A(1,i);
        
    // right hand side is in range of singular matrix
    Vector<> b(size);
    b.SetRandom();
    b(1) = b(0);

    SVDSolver invA(A);
    Vector <> x;
    invA.Apply(b, x);

    Vector<> test;
    A.Apply(x,test);
    cout << "SVD Solver, square singular matrix:" << endl;
    cout << "original right hand side " << b << endl;
    cout << "invA * b = " << x << endl;
    cout << "A * invA * b = " << test << endl;
    cout << "singular values " << invA.GetSigma() << endl;


}


void SVD_overdetermined(){
    // Problem 3 - solve overdetermined system
    // polynomial regression, number of points = 5, polynomial order 2

    int n_points = 5;
    int order_poly = 2;
    Matrix<> A(n_points, order_poly+1, TMatrixStorage::COL_MAJOR);
    Vector<double> x_vals(n_points);
    Vector<double> y_vals(n_points);
    Vector<double> c(order_poly+1);

    x_vals(0) = 0.; y_vals(0) = 2;
    x_vals(1) = 1.; y_vals(1) = 2;
    x_vals(2) = 2.; y_vals(2) = 6;
    x_vals(3) = 3.; y_vals(3) = 7;
    x_vals(4) = 5.1; y_vals(4) = 7;

    for (int i=0; i<n_points; i++)
    {
        double x = 1.;
        for(int j=0; j<=order_poly; j++)
        {
            A(i,j) = x;
            x *= x_vals(i);
        }
    }

    SVDSolver invA(A);
    invA.Apply(y_vals, c);
    
    cout << "SVD Solver, overdetermined system:" << endl;
    cout << "coefficient vector: invA * y = " << c << endl;
    cout << "singular values " << invA.GetSigma() << endl;

    ofstream outfile("plot_poly_regression.py");
    outfile << "import matplotlib.pyplot as plt\nimport numpy as np" << endl << endl;
    outfile << "xvals = np.zeros(" << n_points << ")\n";
    outfile << "yvals = np.zeros(" << n_points << ")\n";
    for (int i=0; i<n_points; i++) outfile << "xvals[" << i << "] = " << x_vals(i) << endl; 
    for (int i=0; i<n_points; i++) outfile << "yvals[" << i << "] = " << y_vals(i) << endl; 
    outfile << "xvals_fine = np.linspace(np.min(xvals), np.max(xvals), 100)\n\n";
    outfile << "def f_reg(x):\n\treturn ";
    for (int i=0; i<=order_poly; i++)
    {
        outfile << c(i) << "* x**" << i;
        if (i < order_poly ) outfile << " + ";
    }
    outfile << endl << endl;;

    outfile << "plt.plot(xvals, yvals, 'x')\n";
    outfile << "plt.plot(xvals_fine, f_reg(xvals_fine))\n";
    outfile << "plt.show()" << endl;;

    outfile.close();


}

void SVD_underdetermined(){
    // Problem 4 - solve underdetermined
    int height = 6;
    int width = 8;

    Matrix<> A(height, width, TMatrixStorage::COL_MAJOR);
    A.SetRandom(-1,1);


    // right hand side vector
    Vector<> b(height);
    b.SetAll(3.);

    SVDSolver invA(A);
    Vector <> x;
    invA.Apply(b, x);

    Vector<> test;
    A.Apply(x,test);
    cout << "SVD Solver, underdetermined system:" << endl;
    cout << "(pseudo)invA * b = " << x << endl;
    cout << "A * invA * b = " << test << endl;
    cout << "singular values " << invA.GetSigma() << endl;


}

int main()
{
    squareSVD();
    squareSVD_singular();
    SVD_overdetermined();
    SVD_underdetermined();
}