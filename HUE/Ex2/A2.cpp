#include <stdio.h>
#include "HermiteBandMatrix.h"
#include <complex>
#include "SCvector.h"

int main()
{
    int n = 10;
    int b = 4;

    // Unit-Test
    n = 6;
    b = 3;
    SC::HermiteBandMatrix<std::complex<double>> A(n, b);

    // Set diagonal
    for (int i = 0; i < n; i++)
    {
        A.Set(i, i, std::complex<double>(i + 1, 0));
    }
    A.Set(0, 1, std::complex<double>(2, 3));
    A.Set(0, 2, std::complex<double>(1, 1));
    A.Set(1, 2, std::complex<double>(1.4, 0));
    A.Set(1, 3, std::complex<double>(0, 3));
    A.Set(2, 3, std::complex<double>(1, -1));
    A.Set(2, 4, std::complex<double>(0, 2));
    A.Set(3, 4, std::complex<double>(1, 0));
    A.Set(3, 5, std::complex<double>(1, 0));
    A.Set(4, 5, std::complex<double>(1, -2));

    A.Print(std::cout);

    SC::Vector<std::complex<double>> x(n);
    SC::Vector<std::complex<double>> r(n);

    x(0) = std::complex<double>(1, 0);
    x(1) = std::complex<double>(0, 1);
    x(2) = std::complex<double>(-1, 0);
    x(3) = std::complex<double>(0, -1);
    x(4) = std::complex<double>(1, 0);
    x(5) = std::complex<double>(0, 1);

    A.Apply(x, r);
    r.Print(std::cout);
    std::cout << "\n";
    A.ApplyT(x, r);
    r.Print(std::cout);
    std::cout << "\n";
    A.ApplyH(x, r);
    r.Print(std::cout);
    std::cout << "\n";
}