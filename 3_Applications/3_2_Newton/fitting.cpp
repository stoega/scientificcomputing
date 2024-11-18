#include <iostream>
#include <fstream>
#include <complex>

#include <SCmatrix.h>
#include <SCvector.h>
#include <SCQRSolver.h>

using namespace std;
using namespace SC;

class MeasurementData
{
    // measurement data
    int n_data;
    Vector<> omega;
    Vector<> E_storage;
    Vector<> E_loss;

public:
    MeasurementData() : n_data(0), E_storage(0), E_loss(0) { ; };
    MeasurementData(const MeasurementData&) = delete;
    MeasurementData& operator=(const MeasurementData&) = delete;
    ~MeasurementData() { ; }

    void ReadData(istream& infile)
    {
        infile >> n_data;

        omega.SetSize(n_data);
        E_storage.SetSize(n_data);
        E_loss.SetSize(n_data);

        for (int i=0; i<n_data; i++)
        {
            infile >> omega(i) >> E_storage(i) >> E_loss(i);
        }
    }

    int NData() const {return n_data; }
    Vector<>& E_l() {return E_loss; }
    const Vector<>& E_l() const {return E_loss; } 
    Vector<>& E_s() {return E_storage; }
    const Vector<>& E_s() const {return E_storage; }
    Vector<>& Omega() { return omega; }
    const Vector<>& Omega() const { return omega; }
    double& E_l(int i) {return E_loss(i); }
    const double& E_l(int i) const {return E_loss(i); }
    double& E_s(int i) {return E_storage(i); }
    const double& E_s(int i) const {return E_storage(i); }
    double& Omega(int i) { return omega(i); }
    const double& Omega(int i) const { return omega(i); }
};


double EvalEstorage(const Vector<double>& params, double omega)
{
    int n_param = (params.Size()-1)/2;
    double Estorage = fabs(params(0));
    double omega2 = omega*omega;
    for (int i=0; i<n_param; i++)
    {
        double E_i = fabs(params(2*i+1));
        double tau_i = fabs(params(2*i+2));
        Estorage += E_i * omega2*tau_i*tau_i/(1+omega2*tau_i*tau_i);
    }
    return Estorage;
}

double EvalEloss(const Vector<double>& params, double omega)
{
    int n_param = (params.Size()-1)/2;
    double Eloss = 0.;
    double omega2 = omega*omega;
    for (int i=0; i<n_param; i++)
    {
        double E_i = fabs(params(2*i+1));
        double tau_i = fabs(params(2*i+2));
        Eloss += E_i * omega*tau_i/(1+omega2*tau_i*tau_i);
    }
    return Eloss;
}


void EvalDEstorage(const Vector<double>& params, double omega, Vector<double>& DEs)
{
    int n_param = (params.Size()-1)/2;
    DEs.SetSize(params.Size());
    double Estorage = fabs(params(0));
    double DEstorage = params(0)/fabs(params(0));
    double omega2 = omega*omega;
    DEs(0) = DEstorage;
    for (int i=0; i<n_param; i++)
    {
        double E_i = fabs(params(2*i+1));
        double tau_i = fabs(params(2*i+2));
        double DE_i = params(2*i+1)/fabs(params(2*i+1));
        double Dtau_i = params(2*i+2)/fabs(params(2*i+2));
        double denominator = (1+omega2*tau_i*tau_i);
        DEs(1+2*i) = DE_i * omega2*tau_i*tau_i/denominator;
        DEs(2+2*i) = E_i * Dtau_i * (2*omega2*tau_i)/denominator/denominator;
    }
    return;
}

    
void EvalDEloss(const Vector<double>& params, double omega, Vector<double>& DEloss) 
{
    int n_param = (params.Size()-1)/2;
    DEloss.SetSize(params.Size());
    DEloss(0) = 0;
    double omega2 = omega*omega;
    for (int i=0; i<n_param; i++)
    {
        double E_i = fabs(params(2*i+1));
        double tau_i = fabs(params(2*i+2));
        double DE_i = params(2*i+1)/fabs(params(2*i+1));
        double Dtau_i = params(2*i+2)/fabs(params(2*i+2));
        double denominator = (1+omega2*tau_i*tau_i);
        DEloss(1+2*i) = DE_i * omega*tau_i/denominator;
        DEloss(2+2*i) = E_i *Dtau_i * omega*(1-omega2*tau_i*tau_i)/denominator/denominator;
    }
    return;
}


void ComputeF(const Vector<double>& params, const MeasurementData& data, Vector<double>& F)
{
    int n = data.NData();
    F.SetSize(2*n);
    for (int j=0; j<n; j++)
    {
        F(j) = EvalEstorage(params, data.Omega(j))-data.E_s(j);
        F(n+j) = EvalEloss(params, data.Omega(j))-data.E_l(j);
    }
}



void ComputeFprime_numdiff(const Vector<double>& params, const MeasurementData& data, Matrix<double>& Fprime)
{
    int n = data.NData();
    Fprime.SetSize(2*n,params.Size());
    Fprime.SetAll(0.);

    // numeric differentiation -> compare x, x+eps
    // dF/dxi = (F(x1, ... , xi+eps, ..) - F(x))/ eps
    Vector<double> F, Feps;
    double eps = 1e-8;
    Vector<double> params_eps(params);

    ComputeF(params, data, F);

    for (int i=0; i<params.Size(); i++)
    {
        params_eps(i) += eps;
        ComputeF(params_eps, data, Feps);
        for(int j=0; j<2*n; j++)
        {
            Fprime(j,i) = (Feps(j) - F(j))/eps;
        }
        params_eps(i) -= eps;
    }

}


void ComputeFprime(const Vector<double>& params, const MeasurementData& data, Matrix<double>& Fprime)
{
    int n = data.NData();
    Fprime.SetSize(2*n,params.Size());
    Fprime.SetAll(0.);
    Vector<double> DEstorage, DEloss;
    for (int j=0; j<n; j++)
    {
        EvalDEstorage(params, data.Omega(j), DEstorage);
        EvalDEloss(params, data.Omega(j), DEloss);
        for(int i=0; i<params.Size(); i++)
        {
            Fprime(j,i) = DEstorage(i);
            Fprime(n+j,i) = DEloss(i);
        }
    }

}


template<class T, class DATA>
void GaussNewtonSolver(Vector<T> &x, 
                       const DATA& data, 
                       void (*ComputeF)(const Vector<T>& x, const DATA& data, Vector<T>& fvec) , 
                       void (*ComputeFprime)(const Vector<T>& x, const DATA& data, Matrix<T>& fprime), 
                       int maxit = 50, 
                       double acc = 1e-8)
{
    Vector<T> F;
    Matrix<T> Fprime;
    Vector<T> deltax(x);

    Vector<T> FprimeTF(x.Size());

    double damping = 1.;

    // initialize
    ComputeF(x, data, F);
    double initial_dist = F.NormSqr();
    double dist;
    double newdist = initial_dist;
    cout << "initial distance " << initial_dist << endl;
    ComputeFprime(x, data, Fprime);

    // compute the defect F'^T F
    Fprime.ApplyH(F, FprimeTF);
    
    double initial_def = FprimeTF.Norm();
    cout << "initial defect ||F'^T F|| " << initial_def << endl;
    for (int step=0; step<maxit; step++)
    {
        ComputeFprime(x, data, Fprime);
        dist = newdist;

        QRSolver<T> invFprime(Fprime);
        invFprime.Apply(F, deltax);

        x -= damping*deltax;
        ComputeF(x, data, F);
        newdist = F.NormSqr();
        
        
        while (newdist > dist && damping > 1e-6)
        {
            damping *= 0.5;
            x += damping*deltax;
            ComputeF(x, data, F);
            newdist = F.NormSqr();
        }
        

        Fprime.ApplyH(F, FprimeTF);
        cout << "step " << step << ": defect ||F'^T F|| " << FprimeTF.Norm() <<
        ", damping factor " << damping << endl;

        damping = min(damping*2, 1.);

        if (FprimeTF.Norm() < acc*initial_def) break;
    }
}


int main()
{
    // path to data file:
    // on windows using MVSC, standard location ../../ (executable in build/Debug or build/Release)
    // on macos using makefiles, standard location ../ (executable in build)
    string path_to_data("../"); 
    // sting path_to_data("../../");
    string datafile("data_20C.txt");
    ifstream infile(path_to_data+datafile);
    MeasurementData data;
    data.ReadData(infile);
    infile.close();

    
    double scale = 1e-6;
    data.E_l() *= scale;
    data.E_s() *= scale;
    

    int n_params = 3;

    // nonlinear least-squares fit
    // fit E_infty, E_i, tau_i

    Vector<> params(1+2*n_params);
    double tau0 = 1.;
    params(0) = 2e6*scale;
    for (int i=0; i<n_params; i++)
    {
        params(2*i+1) = 1e5*scale;
        params(2*i+2) = tau0;
        tau0 *= 0.1;
    }
    
    GaussNewtonSolver(params, data, ComputeF, ComputeFprime);

    ofstream outfile(path_to_data+"computed.txt");
    int np = 100;
    double omega_min = 6., omega_max = 400.;
    for (int i=0; i<np; i++)
    {
        double om = omega_min + i*(omega_max-omega_min)/(np-1.);
        outfile << om << "\t"
        << EvalEstorage(params, om)/scale << "\t"
        << EvalEloss(params, om)/scale << endl;
    }
    outfile.close();
    

    Vector<> F;
    
    ComputeF(params, data, F);

    cout << "parameters nonlinear fit: " << endl;
    cout << " E_infty = " << fabs(params(0)) << endl;
    for (int i=0; i<n_params; i++)
        cout << "E/tau " << i << " : " << fabs(params(2*i+1)) << ", " << fabs(params(2*i+2)) << endl;
    cout << "distance: " << F.NormSqr() << endl;
    

}
