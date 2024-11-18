#include <iostream>

// std:: is visible below
using namespace std;

namespace A
{
    void myfunction() { cout << "A::myfunction" << endl; }
}

namespace B
{
    void myfunction() { cout << "B::myfunction" << endl; }
}

// A:: is visible below, i.e. A::myfunction and myfunction are equivalent
using namespace A;
int main()
{
    // B:: is visible within the scope of main
    // using namespace B;

    A::myfunction();
    B::myfunction();

    // compiler error if A:: and B:: are both used (ambiguity) or none is used (not defined)
    myfunction();
    // try renaming myfunction -> function: ambiguity with std::function pops up

}
