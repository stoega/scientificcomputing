
#define _USE_MATH_DEFINES
#include <math.h>
#include <complex>
#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include <vector>

#include <Eigen/Eigenvalues>
#include <Eigen/Sparse>
#include <Eigen/SparseLU>

using namespace std;

using Eigen::placeholders::all;

// use SPARSEMAT equivalently to Eigen::SparseMatrix<double, Eigen::RowMajor>
typedef Eigen::SparseMatrix<double, Eigen::RowMajor> SPARSEMAT;

template <typename VEC>
void PrintEigenfunctionToFile(ostream& out, int n_x, double lambda, VEC v)
{
    out << lambda << endl;
    for (int j=0; j<n_x+2; j++)
        out << "0 ";
    out << "\n";
    for (int i=0; i<n_x; i++)
    {
        out << "0 ";
        for (int j=0; j<n_x; j++)
            out << v[i*n_x+j] << " ";
        out << "0\n";
    }
    for (int j=0; j<n_x+2; j++)
        out << "0 ";
    out << endl;
}

void InverseIteration(SPARSEMAT& K, SPARSEMAT& M, double lambda_shift, double &lambda1, Eigen::VectorXd &v1)
{
    // pre-processing: LU decomposition for K - lambda_shift I = K - lambda_shift M
    int n = K.rows();
    for (int i=0; i<n; i++) K.coeffRef(i,i) -= lambda_shift;

    Eigen::SparseLU<SPARSEMAT> Kinv;
    Kinv.analyzePattern(K);
    Kinv.factorize(K);

    v1.setRandom(n);
    v1.normalize();
    lambda1 = 0;
    Eigen::VectorXd w(n), r(n), Mv1(n), Kv1(n);
    // initialize
    Mv1 = M*v1;
    Kv1 = K*v1;

    Eigen::VectorXd error;
    error = Kv1 - lambda1*Mv1;
    double error_init = error.norm();
    int step = 0;
    double Kvv, Mvv;
    double lambda_old;
    while (step<100 && error.norm() > 1e-10*error_init)
    {
        step++;
        // next iterate: w = K^-1 M v, v = w/normw
        w = Kinv.solve(Mv1);
        double normw = w.norm();

        error = Kv1 - lambda1*Mv1;
        // update in correct order:
        // Kv1 = K*v1 = K*w/normw = Mv1/normw;
        Kv1 = 1./normw*Mv1;
        v1 = w.normalized();
        Mv1 = M*v1;

        // compute eigenvalue through Rayleigh quotient
        Mvv = Mv1.dot(v1);
        Kvv = Kv1.dot(v1);
        lambda_old = lambda1;
        lambda1 = Kvv/Mvv;


        cout << "step " << step << " lam " << lambda1+lambda_shift << " error " << error.norm() << endl;
    }
    cout << "Inverse iteration: " << step <<  " steps, lambda = " << lambda1+lambda_shift << ", error " << error.norm() << endl;


    // shift back K + lambda_shift M
    for (int i=0; i<n; i++) K.coeffRef(i,i) += lambda_shift;
    lambda1 += lambda_shift;
    
}

// lambda1 must be initialized!
void InverseIterationRayleigh(SPARSEMAT& K, SPARSEMAT& M, double lambda_shift, double &lambda1, Eigen::VectorXd &v1)
{
    // pre-processing: LU decomposition for K - lambda_shift I
    int n = K.rows();
    for (int i=0; i<n; i++) K.coeffRef(i,i) -= lambda_shift;
    Eigen::SparseLU<SPARSEMAT> Kinv;
    Kinv.analyzePattern(K);
    Kinv.factorize(K);

    v1.setRandom(n);
    v1.normalize();
    Eigen::VectorXd w(n), r(n), Mv1(n), Kv1(n), Kw(n), Mw(n);
    // initialize
    Mv1 = M*v1;
    Kv1 = K*v1;

    double error_init = 1;
    double error = 1;
    int step = 0;
    double Kvv, Mvv;
    double lambda_old;

    Eigen::Matrix2d M2, K2;
    Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::Matrix2d> es;

    while (step<100 && error > 1e-10*error_init)
    {
        step++;
        // next update direction: w = v - lambda K^-1 M v
        w = Kinv.solve(Mv1);
        w *= -lambda1;
        w += v1;
        error = w.norm();
        if (step==1) error_init = error;

        // update Kw, Mw:
        Kw = K*w;
        Mw = M*w;

        M2(0,0) = v1.dot(Mv1);
        M2(0,1) = M2(1,0) = w.dot(Mv1);
        M2(1,1) = w.dot(Mw);
        
        K2(0,0) = v1.dot(Kv1);
        K2(0,1) = K2(1,0) = w.dot(Kv1);
        K2(1,1) = w.dot(Kw);
        
         
        es.compute(K2, M2);
        const Eigen::Vector2d& lam2 = es.eigenvalues();
        const Eigen::Matrix2d& y = es.eigenvectors();

        // cout << "eigenproblem " << K2 << endl << M2 << endl << lam2 << endl << y << endl;
        // eigenvalues are ordered by value (not by absolute value) -- choose the smaller one
        int i2 = abs(lam2[0])<abs(lam2[1]) ? 0 : 1;
        v1 *= y(0,i2);
        v1 += y(1,i2)*w;
        v1.normalize();

        // update Kv1, Mv1
        Mv1 = M*v1;
        Kv1 = K*v1;

        lambda_old = lambda1;
        lambda1 = lam2[i2];

        cout << "step " << step << " lam " << lambda1 << " error " << error << endl;
    }


    // shift back
    for (int i=0; i<n; i++) K.coeffRef(i,i) += lambda_shift;
    lambda1 += lambda_shift;

    cout << "Inverse iteration (Rayleigh): " << step <<  " steps, lambda = " << lambda1+lambda_shift << ", error " << error << endl;
    
}

void InverseIterationNRayleigh(SPARSEMAT& K, SPARSEMAT& M, Eigen::VectorXd &lambda1, Eigen::MatrixXd &V)
{
    // pre-processing: LU decomposition for K - lambda_shift I
    int n = K.rows();
    int m = lambda1.rows();
    cout << "m =  " << m << endl;

    Eigen::SparseLU<Eigen::SparseMatrix<double, Eigen::RowMajor>> Kinv;
    Kinv.analyzePattern(K);
    Kinv.factorize(K);

    V.setRandom(n,m);
    for (int i=0; i<m; i++)
        V(all,i).normalize();

    Eigen::MatrixXd W(n,m), R(n,n), MV(n,m), KV(n,m), KW(n,m), MW(n,m), Vnew(n,m);
    // initialize
    MV = M*V;
    KV = K*V;

    double error_init=m;
    Eigen::VectorXd error(m);
    error.setConstant(1.);

    int step = 0;
    Eigen::VectorXd lambda_old(m);

    Eigen::MatrixXd Msmall(2*m,2*m), Ksmall(2*m,2*m);


    Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXd> es;

    while (step<100 && error.norm() > 1e-8*error_init)
    {
        step++;
        // next update direction: w_i = v_i - lambda_i K^-1 M v_i

        W = Kinv.solve(MV);
        for (int i=0; i<m; i++)
            W(all, i) *= -lambda1(i);
        W += V;

        // use residual w as error indicator
        for (int i=0; i<m; i++)
            error(i) = W(all,i).norm();
        if (step==1) error_init = error.norm();

        KW = K*W;
        MW = M*W;

        for (int i=0; i<m; i++)
            for (int j=0; j<m; j++)
            {
                Ksmall(i,j) = V(all,i).dot(KV(all,j));
                Ksmall(i,m+j) = W(all,j).dot(KV(all,i));
                Ksmall(m+j,i) = Ksmall(i,m+j);
                Ksmall(m+i,m+j) = W(all,i).dot(KW(all,j));

                Msmall(i,j) = V(all,i).dot(MV(all,j));
                Msmall(i,m+j) = W(all,j).dot(MV(all,i));
                Msmall(m+j,i) = Msmall(i,m+j);
                Msmall(m+i,m+j) = W(all,i).dot(MW(all,j));

            }
        
        es.compute(Ksmall, Msmall);
        const auto& lam_small = es.eigenvalues();
        const auto& y_small = es.eigenvectors();

        // cout << "eigenproblem " << Ksmall << endl << Msmall << endl;
        // cout << "eigenvalues " << lam_small << endl << y_small << endl;

        Vnew = W*y_small(Eigen::seq(m,2*m-1),Eigen::seq(0,m-1));
        Vnew += V*y_small(Eigen::seq(0,m-1),Eigen::seq(0,m-1));
        for (int j=0; j<m; j++)
            V(all,j) = Vnew(all,j).normalized();

        KV = K*V;
        MV = M*V;


        lambda_old = lambda1;
        lambda1 = lam_small(Eigen::seq(0,m-1));

        cout << "step " << step << " error " << error.norm() << ":\nlam ";
        for (int i=0; i<m; i++) cout << lambda1(i) << " ";
        cout << endl;
    }
    cout << "Inverse iteration: " << step <<  " steps, error " << error.norm() << endl;
    cout << "Computed eigenvalues " << lambda1 << endl;


    
}




int main()
{
    int n_x = 100;
    int n = n_x*n_x;

    double meshsize = 1./(n_x+1);
    double c2_coeff = 1.;

    double meshsize2 = meshsize*meshsize;

    // stiffness matrix K
    Eigen::SparseMatrix<double, Eigen::RowMajor> K(n, n);
    Eigen::SparseMatrix<double, Eigen::RowMajor> M(n, n);
    vector<Eigen::Triplet<double> > tripletList, tripletListM;
    tripletList.reserve(5*n);
    tripletListM.reserve(n);
    double jj=4.;
    for(int i=0; i<n_x; i++)
    {
        for (int j=0; j<n_x; j++)
        {
            tripletList.push_back(Eigen::Triplet<double>(i+j*n_x, i+j*n_x, 4.*c2_coeff/meshsize2));
            if (j>0) tripletList.push_back(Eigen::Triplet<double>(i+(j-1)*n_x, i+j*n_x, -c2_coeff/meshsize2));
            if (j<n_x-1) tripletList.push_back(Eigen::Triplet<double>(i+(j+1)*n_x, i+j*n_x, -c2_coeff/meshsize2));
            if (i<n_x-1) tripletList.push_back(Eigen::Triplet<double>(i+1+(j)*n_x, i+j*n_x, -c2_coeff/meshsize2));
            if (i>0) tripletList.push_back(Eigen::Triplet<double>(i-1+(j)*n_x, i+j*n_x, -c2_coeff/meshsize2));
            tripletListM.push_back(Eigen::Triplet<double>(i+j*n_x, i+j*n_x, 1.));
            jj += 1.;
        }
    }

    K.setFromTriplets(tripletList.begin(), tripletList.end());
    M.setFromTriplets(tripletListM.begin(), tripletListM.end());

    Eigen::Matrix<double, Eigen::Dynamic, 1> b(n), u(n);

    // number of EV to be computed
    int m = 4;
    // shift parameter - compute eigenvector closest to lambda_shift
    // implemented only for single eigenvector computation
    double lambda_shift = 0.;


    if (m==1)
    {
        Eigen::VectorXd v1;
        double lambda1 = 1.;

        InverseIteration(K, M, lambda_shift, lambda1, v1);
        // InverseIterationRayleigh(K, M, lambda_shift, lambda1, v1);

        cout << "lambda = " << lambda1 << ", omega = " << sqrt(lambda1) << endl;
        stringstream filename;
        filename << "eigenproblem.txt";
        ofstream outfile(filename.str());
        PrintEigenfunctionToFile(outfile, n_x, lambda1, v1);
        outfile.close();
    }
    else
    {
        Eigen::MatrixXd V(n,m);
        Eigen::VectorXd lambda(m);
        lambda.setConstant(1.);
        InverseIterationNRayleigh(K, M, lambda, V);



        for (int i=0; i<m; i++)
        {
            cout << "lambda_" << i << " = " << lambda(i) << ", omega = " << sqrt(lambda(i)) << endl;
            stringstream filename;
            filename << "eigenproblem" << i << ".txt";
            ofstream outfile(filename.str());
            PrintEigenfunctionToFile(outfile, n_x, lambda(i), V(all,i));
            outfile.close();
        }

        // check -- are all eigenvectors orthogonal?
        Eigen::MatrixXd VTMV(m,m), &MV(V);
        // MV = M*V;
        VTMV = V.transpose()*MV;

        cout << "check orthogonality of eigenvectors:\n" << VTMV << endl;

    }
    

}