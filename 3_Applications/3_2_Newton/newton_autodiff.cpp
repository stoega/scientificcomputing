#include <iostream>
#include <cmath>

#include <SCmatrix.h>
#include <SCvector.h>
#include <SCLUSolver.h>

namespace SC
{
template <int DIM, typename T>
class Vec
{
private:
    T data[DIM];

public:
    Vec(T val) 
    {
        for (int i=0; i<DIM; i++) data[i] = val;
    }

    ~Vec() = default;

    Vec() 
    {
        for (int i=0; i<DIM; i++) data[i] = T(0);
    }
    Vec(const Vec& vec)
    {
        for (int i=0; i<DIM; i++) data[i] = vec.data[i];
    } 

    Vec& operator=(const Vec& vec)
    {
        for (int i=0; i<DIM; i++) data[i] = vec.data[i];
        return *this;
    } 
    
    int Size() const { return DIM;}

    T& operator()(int i) 
    { 
        #ifndef NDEBUG
        if(i < 0 || i >= DIM)
            throw "Error in Vector::operator(): out of range";
        #endif
        return data[i]; 
    }    
    
    const T& operator()(int i) const 
    { 
        #ifndef NDEBUG
        if(i < 0 || i >= DIM)
            throw "Error in Vector::operator(): out of range";
        #endif
        return data[i]; 
    }    

    // 
    void SetAll(T val) { for (int i=0; i<DIM; i++) data[i] = val; }

};


template <int DIM, typename T>
class AutoDiff
{
private:
    T value;
    Vec<DIM,T> diffvalue;

public:
    // inialize constant function f = val, f' = 0
    AutoDiff(const T& val) : value(val)
    {   
        for (int i=0; i<DIM; i++) diffvalue(i) = 0;
    }

    // initialize coordinate function  f = x_dir, f' = e_dir, and set value x_dir = val
    AutoDiff(T val, int dir) : value(val)
    {   
        for (int i=0; i<DIM; i++) diffvalue(i) = 0;
        diffvalue(dir) = 1;

    }

    // default initialize zero function
    AutoDiff(): value(T(0)) { diffvalue.SetAll(0); }

    AutoDiff(const AutoDiff& ad): value(ad.value), diffvalue(ad.diffvalue) { ; }

    AutoDiff& operator=(const AutoDiff& ad)
    {
        value = ad.value;
        diffvalue = ad.diffvalue;
        return *this;
    }

    T& Value() { return value; }
    const T& Value() const { return value; }
    Vec<DIM,T>& DiffValue() { return diffvalue; }
    const Vec<DIM,T>& DiffValue() const { return diffvalue; }

};

// copy values from Vector<AutoDiff> to dynamic Vector<T> fvec, Matrix<T> jmat
template<int DIM, typename T>
void CopyValues(const Vector<AutoDiff<DIM, T> >& f, Vector<T>& fvec, Matrix<T>& jmat)
{
    fvec.SetSize(f.Size());
    for (int i=0; i<f.Size(); i++)
        fvec(i) = f(i).Value();

    jmat.SetSize(f.Size(), DIM);
    for (int i=0; i<f.Size(); i++)
        for (int j=0; j<DIM; j++)
            jmat(i,j) = f(i).DiffValue()(j);
}

// set Vector<AutoDiff> values according to fvec 
template<int DIM, typename T>
void CopyValues(const Vector<T>& fvec, Vector<AutoDiff<DIM, T>>& f )
{
    f.SetSize(fvec.Size());
    for (int i=0; i<f.Size(); i++)
        f(i).Value() = fvec(i);
}

} // end namespace SC


// operators, functions defined outside namespace SC(!)
template<int DIM, typename T>
SC::AutoDiff<DIM, T> operator+(const SC::AutoDiff<DIM, T>& a, const SC::AutoDiff<DIM, T>& b)
{
    SC::AutoDiff<DIM, T> difference;
    difference.Value() = a.Value() + b.Value();
    for (int i=0; i<DIM; i++) difference.DiffValue()(i) = a.DiffValue()(i) + b.DiffValue()(i);
    return difference;
}

template<int DIM, typename T>
SC::AutoDiff<DIM, T> operator-(const SC::AutoDiff<DIM, T>& a, const SC::AutoDiff<DIM, T>& b)
{
    SC::AutoDiff<DIM, T> difference;
    difference.Value() = a.Value() - b.Value();
    for (int i=0; i<DIM; i++) difference.DiffValue()(i) = a.DiffValue()(i) - b.DiffValue()(i);
    return difference;
}

// product rule
template<int DIM, typename T>
SC::AutoDiff<DIM, T> operator*(const SC::AutoDiff<DIM, T>& a, const SC::AutoDiff<DIM, T>& b)
{
    SC::AutoDiff<DIM, T> product;
    product.Value() = a.Value()*b.Value();
    for (int i=0; i<DIM; i++) product.DiffValue()(i) = a.Value()*b.DiffValue()(i) + a.DiffValue()(i)*b.Value();
    return product;
}

// derive sin
template<int DIM, typename T>
SC::AutoDiff<DIM, T> sin(const SC::AutoDiff<DIM, T>& x)
{
    SC::AutoDiff<DIM, T> sinx;
    sinx.Value() = sin(x.Value());
    for (int i=0; i<DIM; i++) sinx.DiffValue()(i) = cos(x.Value())*x.DiffValue()(i);
    return sinx;
}


using namespace SC;
using namespace std;



template<typename T>
void ScalarNewtonSolver(T& x_init, void (*function)(const AutoDiff<1,T>&, AutoDiff<1,T>&), int maxit = 20, double acc = 1e-8) 
{
    AutoDiff<1,T> x(x_init,0);
    AutoDiff<1,T> f;

    function(x, f);
    double res_init = fabs(f.Value());
    cout << "res_init = " << res_init << endl;

    for (int i=0; i<maxit; i++)
    {
        x.Value() -= 1./(f.DiffValue())(0)*f.Value();


        function(x, f);
        cout << "step " << i << ": res " << fabs(f.Value()) << endl;
        if (fabs(f.Value()) < acc*res_init)
            break;    
    }
    x_init = x.Value();
}


template<int DIM, typename T>
void VectorNewtonSolver(Vector<T>& xvec, void (*function)(const Vector<AutoDiff<DIM,T>>&, Vector<AutoDiff<DIM,T>>&), int maxit = 20, double acc = 1e-8) 
{
    if (DIM != xvec.Size())
        throw "Error in VectorNewtonSolver - dimension of vector x != number of degrees of freedom";

    Vector<AutoDiff<DIM,T>> x(DIM);
    Vector<AutoDiff<DIM,T>> f;
    for(int i=0; i<DIM; i++)
        x(i) = AutoDiff<DIM,T>(xvec(i), i);

    Matrix<T> Jmat(DIM);
    Vector<T> deltaxvec(DIM);
    Vector<T> fvec(DIM);

    function(x, f);
    // copy AutoDiff f to fvec
    CopyValues(f, fvec, Jmat);

    double res_init = fvec.Norm();


    cout << "res_init = " << res_init << endl;

    for (int i=0; i<maxit; i++)
    {
        LUSolver<T> invJ(Jmat);
        invJ.Apply(fvec, deltaxvec);
        xvec -= deltaxvec;

        CopyValues(xvec, x);
        function(x, f);
        CopyValues(f, fvec, Jmat);

        cout << "step " << i << ": res " << fvec.Norm() << endl;
        if (fvec.Norm() < acc*res_init)
            break;    
    }
}




template <class T>
void ComputeF(const T& x, T& f)
{
    f = x*x*x - T(27.);
}

template <class T>
void ComputeFVec(const Vector<T>& x, Vector<T>& f)
{
    f.SetSize(3);
    f(0) = x(0)*x(0)*x(0) + x(0)*x(1) - T(27.);
    f(1) = x(1)*x(1) - x(2)*x(2) - T(4.);
    f(2) = sin(x(2));
}



int main()
{
    double x=0.4;
    ScalarNewtonSolver(x, ComputeF);

    double r;
    ComputeF(x, r);

    cout << "x = " << x << "  f(x) = " << r << endl;


    Vector<double> xvec(3);
    xvec.SetAll(1.);
    VectorNewtonSolver<3>(xvec, ComputeFVec);

    Vector<double> rvec;
    ComputeFVec(xvec, rvec);

    cout << "x = " << xvec << "  f(x) = " << rvec << endl;

}