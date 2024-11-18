#include <iostream>

#ifndef __COMPLEXCLASS_INCLUDEGUARD__
#define __COMPLEXCLASS_INCLUDEGUARD__

// implement a complex number
class Complex
{
    public:

    // two members: real and imaginary part
    double real;
    double imag;

    // the default constructor
    // Complex a;
    Complex()
    { 
        // this compares to python's self, is a pointer to the object itself
        // can be omitted most of the time
        this->real = 0.; 
        this->imag = 0.;
    }

    // constructor for given real and imaginary part, e.g.
    // Complex a(1., 0.);
    Complex(double r, double i)
    {
        real = r; imag = i;
    }

    // copy constructor
    // Complex b(a);
    Complex(const Complex& c)
    {
        real = c.real;
        imag = c.imag;
    }

    // destructor -- is called when complex number is deleted at end of scope/at delete
    ~Complex() 
    { 
        // cout << "calling ~Complex()" << endl;
    }

    // assignment operator
    // c = b;
    // is equivalent to
    // c.operator=(b);
    Complex& operator=(const Complex& c)
    {
        real = c.real;
        imag = c.imag;
        return *this;
    }

    // member functions add functionality to class

    // get the absolute value via c.abs()
    double abs() const
    {   
        return sqrt(real*real + imag*imag);
    }

    // get the argument using atan2
    double arg() const
    {
        return atan2(imag, real);
    }

    // print the complex number -- inline definition
    // void Print(std::ostream& out)
    // {
    //     out << real << "+" << imag << "j";
    // }

    // print the complex number -- declaration, definition below
    void Print(std::ostream& out);

};

// definition of class member Print
void Complex::Print(std::ostream& out)
{
    out << real << "+" << imag << "j";
}

// operator "+" implements c1 + c2;
Complex operator+(const Complex& c1, const Complex& c2)
{   
    return Complex(c1.real+c2.real, c1.imag+c2.imag);
}

// operator "*" implements c1 * c2;
Complex operator*(const Complex& c1, const Complex& c2)
{
    return Complex(c1.real*c2.real - c1.imag*c2.imag, c1.real*c2.imag + c1.imag*c2.real);
}

// operator "<<" for streaming, e.g. cout << c1;
std::ostream& operator<<(std::ostream& out, Complex& c)
{
    c.Print(out);
    return out;
}


#endif // __COMPLEXCLASS_INCLUDEGUARD__
