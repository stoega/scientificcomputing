#include <iostream>

#include <SCmatrix.h>
#include <SCvector.h>
#include <SCQRSolver.h>
#include <SCLUSolver.h>

using namespace std;
using namespace SC;

void f(const double& x, double& r)
{
    r = x*x*x - 27.;
}

void fvec(const Vector<double>& x, Vector<double>& r)
{
    r.SetSize(3);
    r(0) = x(0)*x(0)*x(0) + x(0)*x(1) - 27.;
    r(1) = x(1)*x(1) - x(2)*x(2) - 4.;
    r(2) = sin(x(2));
}

void J(const double& x, double &jac)
{
    jac = 3*x*x;
}

// not efficient, as dynamic memory allocation -- to be seen as proof of concept
void Jmat_numdiff(const Vector<double>& x,Matrix<double> &jac)
{
    double eps = 1e-6;
    jac.SetSize(3,3);
    Vector<double> xeps(x), f(3), feps(3);
    fvec(x, f);

    for (int j=0; j<3; j++)
    {
        xeps(j) += eps;
        fvec(xeps, feps);
        for (int i=0; i<3; i++)
            jac(i,j) = (feps(i) - f(i))/eps;
        xeps(j) -= eps;
    }

}


void Jmat(const Vector<double>& x,Matrix<double> &jac)
{
    jac.SetSize(3,3);
    jac(0,0) = 3*x(0)*x(0) + x(1);
    jac(0,1) = x(0);
    jac(0,2) = 0.;

    jac(1,0) = 0.;
    jac(1,1) = 2*x(1);
    jac(1,2) = -2*x(2);

    jac(2,0) = 0.;
    jac(2,1) = 0;
    jac(2,2) = cos(x(2));

}


template<typename T>
void ScalarFixedPointSolver(T& x, void (*function)(const T&, T&), const T& Pinv, int maxit = 200, double acc = 1e-8) 
{
    T res;
    T deltax;

    double res_init;
    function(x, res);
    res_init = std::abs(res);

    for (int i=0; i<maxit; i++)
    {
        deltax = Pinv*res;
        x -= deltax;
        
        function(x, res);
        cout << "step " << i << ": res " << res << endl;
        if (std::abs(res) < acc*res_init)
            break;
    }
}


template<typename T>
void ScalarNewtonSolver(T& x, void (*function)(const T&, T&), void (*jacobian)(const T&, T&), int maxit = 20, double acc = 1e-8) 
{
    T res;
    T deltax;
    T jac;

    double res_init;
    function(x, res);
    res_init = fabs(res);

    for (int i=0; i<maxit; i++)
    {
        jacobian(x, jac);
        deltax = res/jac;
        x -= deltax;
        
        function(x, res);
        cout << "step " << i << ": res " << res << endl;
        if (fabs(res) < acc*res_init)
            break;
    }
}

template<typename T>
void VectorNewtonSolver(Vector<T> &x, void (*function)(const Vector<T>&, Vector<T>&) , void (*jacobian)(const Vector<T>&, Matrix<T>&), int maxit = 20, double acc = 1e-8)
{
    Vector<T> res;
    Matrix<T> jac;
    Vector<T> deltax(x);

    function(x, res);
    jacobian(x, jac); 
    QRSolver<double> invJ(jac);
    double res_init = res.Norm();
    double res_norm;
    cout << "initial residual " << res_init << endl;

    for (int i=0; i<maxit; i++)
    {
        jacobian(x, jac); 
        invJ.ComputeFactorization();
        invJ.Apply(res, deltax); //deltax = f(x)/J(x);
        x -= deltax;

        function(x, res);
        res_norm = res.Norm();
        cout << "step " << i << ": res " << res_norm << endl;
        if (res_norm < acc * res_init)
            break;
    }
}


int main()
{
    double x=4;
    ScalarFixedPointSolver(x, f, 1e-2);

    /*
    x=4;
    ScalarNewtonSolver(x, f, J);
    */

    double r;
    f(x, r);

    cout << "x = " << x << "  f(x) = " << r << endl;

    /*
    // solve 3x3 system
    Vector<double> xvec(3);
    xvec.SetAll(0.5);
    VectorNewtonSolver(xvec, fvec, Jmat);

    // xvec.SetAll(0.5);
    // VectorNewtonSolver(xvec, fvec, Jmat_numdiff);

    Vector<double> rvec;
    fvec(xvec, rvec);

    cout << "x = " << xvec << "  f(x) = " << rvec << endl;
    */
}