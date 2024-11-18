#define _USE_MATH_DEFINES
#include <math.h>

#include <Eigen/Sparse>
#include <SCvector.h>
#include <SCmatrix.h>
#include <SCSparseMatrix.h>
#include <iostream>

//using namespace Eigen;
using namespace std;
using namespace SC;

int main()
{
     // test matrix
     // [10    0    0    22]
     // [0     4    5    0]
     // [0     0    0    3]
     int m = 3;
     int n = 4;
     int nze =5;
     Vector<double> data_v(nze);
     Vector<int> colind_v(nze);
     Vector<int> rowptr_v(m+1);

     // set rowpointer vector [0, 2, 4, 5]
     rowptr_v(0) = 0;
     rowptr_v(1) = 2;
     rowptr_v(2) = 4;
     rowptr_v(3) = 5;
     
     // set column index vector and data vector
     // [0, 3, 1, 2, 3] and [10., 22., 4., 5., 3.]
     colind_v(0) = 0; data_v(0) = 10.;
     colind_v(1) = 3; data_v(1) = 22.;
     colind_v(2) = 1; data_v(2) = 4.;
     colind_v(3) = 2; data_v(3) = 5.;
     colind_v(4) = 3; data_v(4) = 3.;


     SparseMatrix<double> Amat(m, n, data_v, colind_v, rowptr_v);

     Vector<double> x(n); 
     for (int i=0; i<n; i++) x(i) = i;
     Vector<double> y;

     Amat.Apply(x, y);
     std::cout << "A x = " << y << endl;

     // compare to Eigen sparse matrix class
     Eigen::SparseMatrix<double, Eigen::RowMajor> eigenAmat(m,n);

     std::vector<Eigen::Triplet<double> > tripletList;
     tripletList.reserve(nze);
     tripletList.push_back(Eigen::Triplet<double>(0, 0, 10.));
     tripletList.push_back(Eigen::Triplet<double>(0, 3, 22.));
     tripletList.push_back(Eigen::Triplet<double>(1, 1, 4.));
     tripletList.push_back(Eigen::Triplet<double>(1, 2, 5.));
     tripletList.push_back(Eigen::Triplet<double>(2, 3, 3.));


     eigenAmat.setFromTriplets(tripletList.begin(), tripletList.end());

     Eigen::Matrix<double, Eigen::Dynamic, 1> eigenx(n), eigeny(m);
     for (int i=0; i<n; i++) eigenx[i] = i;

     eigeny = eigenAmat*eigenx;

     std::cout << "Eigen: A x = " << eigeny << endl;


}
