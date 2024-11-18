#include <iostream>

using namespace std;

// declarations of functions BEFORE usage
// if not defined, compiler error
// inline declaration and definition of Sum function
double Sum(double a, double b)
{
    return a + b;
}

// declaration of Product function
double Product(double a, double b);

int main()
{
    double a, b, c, d;
    a = 3.;
    b = 4.;

    // using previously declared functions
    c = Sum(a, b);
    d = Product(a, b);

    cout << "Sum(a,b) = " << c << endl;
    cout << "Product(a,b) = " << d << endl;
}

// definitions of function somewhere else in the code
// if not defined (comment out below), linker error
double Product(double a, double b)
{
    return a*b;
}
