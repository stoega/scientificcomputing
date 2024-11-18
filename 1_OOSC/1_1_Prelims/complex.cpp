#include <iostream>
#include <math.h>

// #include "..." -> local path
#include "complexclass.h"

using namespace std;


int main()
{
    Complex c(1.,2.);
    Complex d, e;

    // access to public members, similar to c struct
    d.real = 1.5;
    d.imag = 2*c.imag;

    // access to public member functions, such as absolute value
    double abs_c = c.abs();
    cout << "abs(c) = " << abs_c << endl;
    cout << "arg(d) = " << d.arg() << endl;

    // output via Print function
    c.Print(cout);
    cout << endl;

    // output via cout and operator<< (see streams section)
    cout << "c = " << c << endl;
    cout << "d = " << d << endl;

    // pointer to Complex
    Complex* c_ptr = &c;
    cout << "address of c" << c_ptr << endl;

    cout << "object at address of c" << *c_ptr << endl;

    // use -> for access to members of class pointer, same as (* ).
    cout << "real part of object at address c_ptr (c_ptr->real) " << c_ptr->real << endl;
    cout << "real part of object at address c_ptr (*c_ptr).real " <<(*c_ptr).real << endl;

    // (de)allocate complex number via new/delete
    Complex* z_ptr;
    z_ptr = new Complex();

    z_ptr->real = 0.;
    z_ptr->imag = 1.;

    e = (*z_ptr) * d;

    delete z_ptr;

    // dynamic allocation of Complex array
    cout << "a dynamic array of complex numbers" << endl;
    Complex* c_array_dyn;
    c_array_dyn = new Complex[5];
    for (int i=0; i<5; i++)
        cout << c_array_dyn[i] << endl;
    delete [] c_array_dyn;


    // --------------------------------
    // operator overloading
    e.operator=(operator+(d, c));
    cout << "d + c = e = " << e << endl;
    e = d + c;
    cout << "d + c = e = " << e << endl;



}