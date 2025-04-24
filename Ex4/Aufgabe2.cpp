#include "SCTridiagSparseMatrix.h"
#include "SCTridiagLUSolver.h"


using namespace SC;
using namespace std;

// LAPACK Funktion für Eigenwertproblem (spezifische eigenwerte)
extern "C"
{
    void dsygvx_(int *itype, char *jobz, char *range, char *uplo, int *n,
                 double *a, int *lda, double *b, int *ldb, double *vl, double *vu, int *il,
                 int *iu, double *abstol, int *m, double *w, double *z, int *ldz, double *work,
                 int *lwork, int *iwork, int *ifail, int *info);

    void dsyev_(char *jobz, char *uplo, int *n, double * a, int *lda, double *w, double *work, int *lwork, int *info);

    void dggev_(const char* JOBVL, const char* JOBVR, const int* N,
        const double* A, const int* LDA, const double* B, const int* LDB,
        double* ALPHAR, double* ALPHAI, double* BETA,
        double* VL, const int* LDVL, double* VR, const int* LDVR,
        double* WORK, const int* LWORK, int* INFO);
}

template <typename T>
/**
 * Fills a matrix of type SC::TridiagSparseMatrix with proviede diagonal and offdiagonal values
 * @param[out] A Matrix to be filled
 * @param[in] diagonalValue Values on the diagonal
 * @param[in] offdiagonalValue Values of the offdiagonal
 */
void FillTridiagSparseMatrix(TridiagSparseMatrix<T> &A, T diagonalValue, T offdiagonalValue)
{
    size_t n = A.Height();
    for (int i = 0; i < n; ++i)
    {
        // Diagonal elements
        A.Set(i, i, diagonalValue);

        // Subdiagonal elements (left)
        if (i > 0)
        {
            A.Set(i, i - 1, offdiagonalValue);
        }

        // Superdiagonal elements (right)
        if (i < n - 1)
        {
            A.Set(i, i + 1, offdiagonalValue);
        }
    }
}



void modifiedInvIterNRayleigh(TridiagSparseMatrix<double> &K, TridiagSparseMatrix<double> &M, Vector<double> &lambda, Vector<double> *V)
{
    int n = K.Height();
    int m = lambda.Size();

    TridiagLUSolver<double> solver(K);

    // Initialisierung der Eigenvektor-Näherungen U (n x m)
    for (int i = 0; i < m; i++)
    {
        V[i].SetSize(n);
        V[i].SetRandom(0, 1, i);
        V[i] *= 1. / V[i].Norm();
    }

    Vector<double> *MV = new Vector<double>[m];
    Vector<double> *KV = new Vector<double>[m];
    Vector<double> *W = new Vector<double>[m];
    Vector<double> *KW = new Vector<double>[m];
    Vector<double> *MW = new Vector<double>[m];
    Vector<double> *Vnew = new Vector<double>[m];

    for (int i = 0; i < m; i++)
    {
        MV[i].SetSize(n);
        KV[i].SetSize(n);
        W[i].SetSize(n);
        KW[i].SetSize(n);
        MW[i].SetSize(n);
        Vnew[i].SetSize(n);
        Vnew[i].SetAll(0);
    }

    // MV_i = M*V_i & KV_i = K*V_i
    for (int i = 0; i < m; i++)
    {
        M.Apply(V[i], MV[i]);
        K.Apply(V[i], KV[i]);
    }

    double error_init = m;
    Vector<double> error(m);
    error.SetAll(1.);

    int step = 0;

    // set up eigenval problem, lapack uses col_major
    Vector<double> lambda_old(m);
    Matrix<double> Ms(2 * m, 2 * m, TMatrixStorage::COL_MAJOR), Ks(2 * m, 2 * m, TMatrixStorage::COL_MAJOR);

    while (step < 100 && error.Norm() > 1e-8 * error_init)
    {
        step++;

        // w_i (3.46)?
        // next update direction: w_i = v_i - lambda_i K^-1 M v_i
        for (int i = 0; i < m; i++)
        {
            solver.Apply(MV[i], W[i]);
            W[i] *= -lambda(i);
            W[i] += V[i];
        }

        // use residual w as error indicator
        for (int i = 0; i < m; i++)
        {
            error(i) = W[i].Norm();
        }

        if (step == 1)
            error_init = error.Norm();

        // calc KW and MW
        for (int i = 0; i < m; i++)
        {
            M.Apply(W[i], MW[i]);
            K.Apply(W[i], MW[i]);
        }

        for (int i = 0; i < m; i++)
        {
            for (int j = 0; j < m; j++)
            {
                // std::cout << V[i](0) << "\t" << KV[j](0) << "\t" << InnerProduct(V[i], KV[j]) << std::endl;
                Ks(i, j) = InnerProduct(V[i], KV[j]);
                Ks(i, m + j) = InnerProduct(W[j], KV[i]);
                Ks(m + j, i) = Ks(i, m + j);
                Ks(m + i, m + j) = InnerProduct(W[i], MW[j]);

                Ms(i, j) = InnerProduct(V[i], KV[j]);
                Ms(i, m + j) = InnerProduct(W[j], KV[i]);
                Ms(m + j, i) = Ms(i, m + j);
                Ms(m + i, m + j) = InnerProduct(W[i], MW[j]);
            }
        }

        // d) Kleines Standard-EWP lösen mit LAPACK: K_proj * Q = Q * D
        // ... (Code unverändert) ...
        int N = m;
        char jobz = 'V';
        char uplo = 'U';
        Vector<double> y_small(N);
        int lda = N;
        int info;
        double wkopt;
        int lwork = -1;
        dsyev_(&jobz, &uplo, &N, &Ks(0, 0), &lda, &y_small(0), &wkopt, &lwork, &info);
        lwork = static_cast<int>(wkopt);
        Vector<double> work(lwork);
        dsyev_(&jobz, &uplo, &N, &Ks(0, 0), &lda, &y_small(0), &work(0), &lwork, &info);
        if (info != 0)
        {
            std::cerr << "Fehler bei LAPACK dsyev: info = " << info << std::endl;
        }
        

        // for (int i = 0; i < 2 * m; i++)
        // {
        //     std::cout << "w(" << i + 1 << ") = " << w(i) << std::endl;
        // }
        // for (int i = 0; i < 2 * m; i++)
        // {
        //     std::cout << "z(" << i + 1 << ") = " << z(i) << std::endl;
        // }
        // std::cout << "------------------------------------" << std::endl;





        // compute Vnew based on y_small and V & W
        for (int s = 0; s < m; s++)
        {
            for (int z = 0; z < n; z++)
            {
                for (int i = 0; i < m; i++)
                {
                    Vnew[s](z) += V[i](z) * y_small(i, s) + W[i](z) * y_small(i + m, s);
                }
            }
        }
    }
}

int main()
{
    // Problem parameters
    const int m = 5;
    const float L = .5;   // length of flute
    const int c = 343;    // speed of sound
    const int n = 100;    // number of unknowns
    const int ne = n + 1; // number of elements

    const double tol = 1e-9;  // error tolerance
    const int maxIter = 1000; // max iterations

    const double h = L / ne; // mesh fineness

    TridiagSparseMatrix<double> K(n);
    FillTridiagSparseMatrix(K, 2.0 * pow(c, 2) / h, -pow(c, 2) / h);

    TridiagSparseMatrix<double> M(n);
    FillTridiagSparseMatrix(M, 4.0 * h / 6.0, h / 6.0);

    SC::Vector<double> V[m];
    SC::Vector<double> lambda(m);
    lambda.SetAll(1.);

    modifiedInvIterNRayleigh(K, M, lambda, V);
    return 0;
}