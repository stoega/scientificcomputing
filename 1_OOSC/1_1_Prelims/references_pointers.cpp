#include <iostream>

using namespace std;

int main()
{
    double d = 5.;

    
    cout << "value of d:   " << d << endl;   // print value of d
    cout << "address of d: " << &d << endl;  // print address of d

    cout << "value at address of d: " <<  *(&d) << endl;  // print value of double stored at address of d



    cout << "---------------------------------------------" << endl;
    cout << "Pointers" << endl;
    double* ptr = &d;
    cout << "value of pointer " << ptr << endl;
    cout << "de-referencing the pointer gives the value " << *ptr << endl;

    cout << "---------------------------------------------" << endl;
    cout << "References" << endl;
    double& refd = d;
    cout << "address of d    " << &d << endl;
    cout << "address of refd " << &refd << endl;

    // compiler error if reference is not referring to existing valid object
    // double& refb;

    cout << "when changing d, refd is affected" << endl;
    cout << "value of refd before increasing d: " << refd << endl;
    // increase d
    d++;
    // print refd
    cout << "value of refd after increasing d: " << refd << endl;
    // default return value of int main()
    return 0;
}

