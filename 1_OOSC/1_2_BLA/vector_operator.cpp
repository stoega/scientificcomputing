#include <iostream>
#include <complex>

#include <chrono>

#include <SCvector.h>

#include <Eigen/Dense>

using namespace std;
using namespace SC;


int main()
{
    int size = 3;
    Vector<> a(size);
    Vector<> b(size);
    Vector<> c(size);

    a.SetAll(0);
    a(0) = 1;
    a(2) = 2;
    a(1) = 3;
    b.SetAll(5);


    // example from lecture notes
    cout << "using operators.." << endl;
    cout << "-----------------" << endl;
    c = 2.*a + b;
    cout << "-----------------" << endl;

    cout << "using member functions.." << endl;
    cout << "-----------------" << endl;
    c = b;
    c.AddMultiple(2., a);
    cout << "-----------------" << endl;

    cout << "*** now check efficiency ***" << endl;
    // compute c = 2*a + b, a = 1/|c|*c repeatedly
    // using * straightforward operators
    //       * member functions
    //       * Eigen lib - expression templates
    size = 4;
    int count = 1000*1000*10;

    cout << "vector size " << size << endl;
    cout << "number of repeats " << count << endl;

    a.SetSize(size); b.SetSize(size); c.SetSize(size);
    a.SetAll(0);
    a(0) = 1; a(2) = 2; a(1) = 3;
    b.SetAll(5);

    auto begin_op = std::chrono::high_resolution_clock::now();

    for (int i=0; i<count; i++)
    {
        c = 2.*a + b;
        a = 1./c.Norm()* c;
    }
    auto end_op = std::chrono::high_resolution_clock::now();
    auto elapsed_op = std::chrono::duration_cast<std::chrono::nanoseconds>(end_op - begin_op);
    cout << "a = " << a << endl;

    a.SetAll(0);
    a(0) = 1; a(2) = 2; a(1) = 3;
    b.SetAll(5);
    auto begin_fun = std::chrono::high_resolution_clock::now();
    for (int i=0; i<count; i++)
    {
        c = b;
        c.AddMultiple(2, a);
        a = c;
        a *= 1./c.Norm();
    }
    auto end_fun = std::chrono::high_resolution_clock::now();
    auto elapsed_fun = std::chrono::duration_cast<std::chrono::nanoseconds>(end_fun - begin_fun);
    cout << "a = " << a << endl;


    // using the eigen lib
    Eigen::Matrix<double, Eigen::Dynamic, 1> aa(size), bb(size);
    Eigen::Matrix<double, Eigen::Dynamic, 1> cc(size);
    bb.setConstant(5.);
    aa.setConstant(0.);
    aa[0] = 1;
    aa[1] = 3.;
    aa[2] = 2.;

    auto begin_expr = std::chrono::high_resolution_clock::now();

    for (int i=0; i<count; i++)
    {
        cc = 2.*aa + bb;
        aa = cc.normalized();
    }
    auto end_expr = std::chrono::high_resolution_clock::now();
    auto elapsed_expr = std::chrono::duration_cast<std::chrono::nanoseconds>(end_expr - begin_expr);

    cout << "aa = " << aa << endl;
    
    cout << "Using straightforward operators, s   " << elapsed_op.count()/1e9 << endl;
    cout << "Using member functions, s            " << elapsed_fun.count()/1e9 << endl;
    cout << "Eigen: using expression templates, s " << elapsed_expr.count()/1e9 << endl;

    return 0;
}