#include <iostream>

using namespace std;


int main()
{
    // declare different objects a, d, e, b
    int a;
    double d;
    char e;
    bool b;
    size_t len;

    // set values accoring to type
    d = 4.5;
    a = d;
    e = 'e';
    b = true;
    len = 1000000000000;
    

    printf("as char: %c\n", e);
    printf("as int: %i\n", e);
    printf("as hex: %x\n", e);

    // cout is used to write to console, endl means newline
    cout << "int a  = " << a << endl;
    cout << "double d  = " << d << endl;
    cout << "char e  = " << e << endl;
    cout << "bool b  = " << b << endl;
    cout << "size_t len = " << len << endl;

    // double d has already been declared, cannot be redeclared within the same scope.
    // however, it is possible to declare int d within its own scope (use {})
    // outside this scope, d is still double and not changed
    {
        int d;
        d = 2;
        cout << "int d = " << d << endl;
    }

    cout << "double d = " << d << endl;

    cout << endl << "range of integer types:" << endl;
    // scope of integer types: 
    // when setting outside range, over/underflow occurs
    cout << "int ranges from -2147483648 to 2147483647: " << endl;
    a = 3000000000;
    cout << "int a  = " << a << " != 3000000000" <<  endl;

    cout<< "size_t ranges from 0 to the maximum size of a theoretically possible object of any type" << endl;
    len =3000000000;
    cout << "size_t len  = " << len << " == 3000000000" <<  endl;
    len = -1;
    cout << "size_t len  = " << len << " != -1" <<  endl;

    cout << "sizeof " << sizeof(b) << endl;

    cout << "standard mistake 1" << endl;
    // int/int = int
    d = 1/3;
    // implicit conversion to double -> works
    // d = 1./3;
    // explicit conversion to double -> works
    // d = double(1)/3;
    cout << "1/3 = " << d << endl;

    cout << "standard mistake 2" << endl;
    // int*int = int
    len = 1000*1000*1000*1000;
    cout << "1000*1000*1000*1000 = " << len << endl;

    cout << "standard mistake 3" << endl;
    // signed and unsigned
    size_t var1 = 10, var2 = 20;
    double difference = var1 - var2;
    cout << "10 - 20 = " << difference << endl;
}
