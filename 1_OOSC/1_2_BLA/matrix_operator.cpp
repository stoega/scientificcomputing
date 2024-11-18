#include <iostream>
#include <complex>
#include <chrono>

#include <Eigen/Dense>

#include "SCvector.h"
#include "SCmatrix.h"

using namespace std;
using namespace SC;

typedef complex<double> COMPLEX;    // write COMPLEX instead of complex<double>
enum { DS = Eigen::Dynamic }; // Eigen-lib dynamic size identifier

template <typename T>
inline SC::Vector<T> operator*(const SC::Matrix<T>& mat, const SC::Vector<T>& vec)
{
    SC::Vector<T> multvec(vec);
    mat.Apply(vec, multvec);
    return multvec;
}

int main()
{
    // compute b = M.a, a += b repeatedly
    // using * straightforward operators
    //       * member functions
    //       * Eigen lib - expression templates
    int size = 4;
    Vector<COMPLEX> a, b, c;
    a.SetSize(size); a.SetAll(1.);

    Matrix<COMPLEX> M(size);
    M.SetAll(0.1);
    M(2,2) = COMPLEX(0.,1.);
    M(0,2) = COMPLEX(0.,1.);

    cout << "Timings: " << endl;
    int count = 1000*1000;
    cout << "matrix/vector size " << size << endl;
    cout << "number of repeats " << count << endl;

    // Test 1 - count matrix-vector multiplications, using Apply
    auto begin_app = std::chrono::high_resolution_clock::now();

    for (int i=0; i<count; i++)
    {
        M.Apply(a, b);
        a += b;
    }
    auto end_app = std::chrono::high_resolution_clock::now();
    auto elapsed_app = std::chrono::duration_cast<std::chrono::nanoseconds>(end_app - begin_app);
    cout << "Using apply routine, s                " << elapsed_app.count()/1e9 << endl;

    // Test 2 - count matrix-vector multiplications, using operator* and operator=
    auto begin_op = std::chrono::high_resolution_clock::now();
    for (int i=0; i<count; i++)
    {
        b = M * a;
        a += b;
    }
    auto end_op = std::chrono::high_resolution_clock::now();
    auto elapsed_op = std::chrono::duration_cast<std::chrono::nanoseconds>(end_op - begin_op);
    cout << "Using operators, s                    " << elapsed_op.count()/1e9 ;
    cout << "\tslower by factor " << (1.*elapsed_op.count())/elapsed_app.count() << endl;

    // Test 3 - count matrix-vector multiplications, expression templates from Eigen lib
    Eigen::Matrix<COMPLEX, DS, 1> aa(size), bb;
    Eigen::Matrix<COMPLEX, DS, DS> MM(size,size);
    MM.setConstant(0.1);
    MM(3,2) = COMPLEX(0.,1.);
    MM(0,2) = COMPLEX(0.,1.);
    aa.setConstant(1.);

    auto begin_expr = std::chrono::high_resolution_clock::now();
    for (int i=0; i<count; i++)
    {
        bb = MM * aa;
        aa += bb;
    }
    auto end_expr = std::chrono::high_resolution_clock::now();
    auto elapsed_expr = std::chrono::duration_cast<std::chrono::nanoseconds>(end_expr - begin_expr);


    cout << "Using expression templates (Eigen), s " << elapsed_expr.count()/1e9 ;
    cout << "\tcompares by factor " << (1.*elapsed_expr.count())/elapsed_app.count() << endl;


    return 0;
}