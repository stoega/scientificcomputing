#include <iostream>
#include <vector>
#include <cmath>
#include "SCTridiagSparseMatrix.h" // Die Tridiagonalmatrix-Klasse
#include "TridiagLUSolver.h"       // Die Solver-Klasse
#include "SCvector.h"

using namespace SC;

int main()
{
    // Problemparameter
    const double k = 10.0;          // Diffusionskoeffizient
    const double b = 50.0;          // Konvektionskoeffizient
    const int n = 100;              // Anzahl der Unbekannten (Diskretisierungspunkte)
    const double h = 1.0 / (n + 1); // Gitterabstand

    // Erstelle Tridiagonalmatrix
    TridiagSparseMatrix<double> A(n);

    // Fülle die Tridiagonalmatrix gemäß den Gleichungen
    for (int i = 0; i < n; ++i)
    {
        double xi = (i + 1) * h;

        // Diagonalelemente
        A.Set(i, i, 2.0 * k / (h * h) + b / h);

        // Subdiagonalelemente
        if (i > 0)
        {
            A.Set(i, i - 1, -k / (h * h) - b / h);
        }

        // Superdiagonalelemente
        if (i < n - 1)
        {
            A.Set(i, i + 1, -k / (h * h));
        }
    }

    // Erstelle rechte Seite b
    Vector<double> rhs(n);
    rhs.SetAll(0.0);
    for (int i = n / 2; i < n; ++i)
    { // f(x) = 1 für x >= 0.5
        rhs(i) = 1.0;
    }

    // Löse das System Ax = b mit TridiagLUSolver
    TridiagLUSolver<double> solver(A);
    Vector<double> solution(n);
    solution.SetAll(0.0);

    solver.Apply(rhs, solution);

    // Ausgabe der Lösung
    std::cout << "Lösung u(x):" << std::endl;
    for (int i = 0; i < n; ++i)
    {
        double xi = (i + 1) * h;
        std::cout << "x = " << xi << ", u(x) = " << solution(i) << std::endl;
    }

    return 0;
}
