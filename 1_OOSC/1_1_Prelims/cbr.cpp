#include <iostream>

using namespace std;

// call by value, does not work
void Increase_cbv(double d)
{
    // print memory address of d inside function
    cout << "CBV: internal address " << &d << endl;
    d++;
    cout << "CBV: internal value of d " << d << endl;
}

// call by reference
void Increase_CBR(double &d)
{
    // print memory address of d inside function
    cout << "CBR: internal address " << &d << endl;
    d++;
    cout << "CBR: internal value of d " << d << endl;
}

// call by reference, pointer argument
void Increase_CBR(double *p)
{
    // print memory address of pointer inside function
    cout << "CBR - ptr: internal address p " << p << endl;
    // increase value
    (*p)++;
    cout << "CBR - ptr: internal value (*p) " << (*p) << endl;
}


int main()
{
    double d = 5.;
    // print memory address of "outside" d
    cout << "CBV: external address " << &d << endl;
    Increase_cbv(d);
    // value of d
    cout << "CBV: after function call: value of d " << d << endl;

    cout << "---------------------------------------------" << endl;
    cout << "Call by reference" << endl;

    // print memory address of "outside" d
    cout << "CBR: external address " << &d << endl;
    Increase_CBR(d);
    // value of d
    cout << "CBR: after function call: value of d " << d << endl;


    cout << "---------------------------------------------" << endl;
    cout << "Call by reference, using pointer argument" << endl;

    // print memory address of "outside" d
    cout << "CBR: external address " << &d << endl;
    Increase_CBR(&d);
    // value of d
    cout << "CBR: after function call: value of d " << d << endl;

    // default return value of int main()
    return 0;
}

