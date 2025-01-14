#include <iostream>
#include <fstream>

#include <SCvector.h>
#include <SCmatrix.h>

#include <SCQRSolver.h>

using namespace SC;
using namespace std;

int main()
{
    // polynomial regression

    int n_points = 5;
    int order_poly = 1;
    Matrix<double> B(n_points, order_poly+1);
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
            B(i,j) = x;
            x *= x_vals(i);
        }
    }

    QRSolver<double> invB_QR(B);
    invB_QR.Apply(y_vals, c);
    
    cout << "coefficient vector of polynomial = " << c << endl;

    cout << "QR factorization " << invB_QR << endl;

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