#include <iostream>

using namespace std;

int main ()
{
    // static allocation
    int a;
    double b[29];

    a = 0;
    for (int i=0; i<29; i++)
    {
        b[i] = i;
        a += b[i];
    }

    // dynamic allocation

    // dynamic size n, not known at compile time (user input)
    int n;
    cout << "n = ";
    cin >> n;

    // dynamically allocated double
    double *c = new double;
    // dynamically allocated double array
    double *d = new double[n];

    *c = 0;
    for (int i=0; i<n; i++)
    {
        d[i] = i;
        *c += d[i];
    }

    cout << "sum of entries of d is " << *c << endl; 

    // free memory of new double c
    delete c;
    // free memory of new double[n] d
    delete [] d;

    // c, d are still double pointers, also the address is still stored however:
    // now memory of c, d can be overwritten, or newly allocated
    c = new double[n];
    //...
    c[0] = 10;
    c[1] = 2;
    c[2] = 4;


    cout << "d has been deleted ... contained d = [0,1,...n]\nexperimental: see what is stored in d[0] now: " << d[0] << endl;

    delete [] c;


}