#include <iostream>

#include <SCvector.h>
#include <SCmatrix.h>

using namespace SC;
using namespace std;

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
            infile >> omega[i] >> E_storage[i] >> E_loss[i];
        }
    }

    int NData() {return n_data; }
    Vector<>& E_l() {return E_loss; }
    Vector<>& E_s() {return E_storage; }
    Vector<>& Omega() { return omega; }
    double& E_l(int i) {return E_loss(i); }
    double& E_s(int i) {return E_storage(i); }
    double& Omega(int i) { return omega(i); }
};

class ViscoelasticModel
{
protected:
    int n_param;

public:
    ViscoelasticModel() : n_param(0) { ; }
    ViscoelasticModel(int n_p) : n_param(n_p) { ; }
    ViscoelasticModel(const ViscoelasticModel&) = delete;
    ViscoelasticModel& operator=(const ViscoelasticModel&) = delete;
    ~ViscoelasticModel() { ; }


    virtual double GetEinfty(const Vector<double>& params) = 0;
    virtual double GetEi(int i, const Vector<double>& params) = 0;
    virtual double GetTaui(int i, const Vector<double>& params) = 0;
    
    virtual double evalEstorage(const Vector<double>& params, double omega)
    {
        double Estorage = GetEinfty(params);
        double omega2 = omega*omega;
        for (int i=0; i<n_param; i++)
        {
            double E_1 = GetEi(i, params);
            double tau_1 = GetTaui(i, params);
            Estorage += E_1 * omega2*tau_1*tau_1/(1+omega2*tau_1*tau_1);
        }
        return Estorage;
    }

    virtual double evalEloss(const Vector<double>& params, double omega)
    {
        double Eloss = 0.;
        double omega2 = omega*omega;
        for (int i=0; i<n_param; i++)
        {
            double E_1 = GetEi(i, params);
            double tau_1 = GetTaui(i, params);
            Eloss += E_1 * omega*tau_1/(1+omega2*tau_1*tau_1);
        }
        return Eloss;
    }

    virtual void evalDEstorage(const Vector<double>& params, double omega, Vector<double>& DEs)
    {
        throw "ViscoelasticModel::evalDEstorage not implemented in base class";
    }

        
    virtual void evalDEloss(const Vector<double>& params, double omega, Vector<double>& DEloss)
    {
        throw "ViscoelasticModel::evalDEloss not implemented in base class";
    }
};

class LinearViscoelasticModel : public ViscoelasticModel
{
    Vector<> taui;

public:
    using ViscoelasticModel::n_param;

    LinearViscoelasticModel() : ViscoelasticModel(), taui(0) { ; }
    LinearViscoelasticModel(int n_p, const Vector<> tau_i) : ViscoelasticModel(n_p), taui(tau_i) { ; }
    LinearViscoelasticModel(const LinearViscoelasticModel&) = delete;
    LinearViscoelasticModel& operator=(const LinearViscoelasticModel&) = delete;
    ~LinearViscoelasticModel() { ; }

    double GetEinfty(const Vector<double>& params) override { return (params(0)); }
    double GetEi(int i, const Vector<double>& params) override { return (params(i+1)); }
    double GetTaui(int i, const Vector<double>& params) override { return taui(i); }

    void evalDEstorage(const Vector<double>& params, double omega, Vector<double>& DEs) override 
    {
        DEs.SetSize(params.Size());
        double Estorage = GetEinfty(params);
        double omega2 = omega*omega;
        DEs(0) = 1.;
        for (int i=0; i<n_param; i++)
        {
            double E_1 = GetEi(i, params);
            double tau_1 = GetTaui(i, params);
            double denominator = (1+omega2*tau_1*tau_1);
            DEs(1+i) = omega2*tau_1*tau_1/denominator;
        }
        return;
    }

        
    void evalDEloss(const Vector<double>& params, double omega, Vector<double>& DEloss) override 
    {
        DEloss.SetSize(params.Size());
        DEloss(0) = 0;
        double omega2 = omega*omega;
        for (int i=0; i<n_param; i++)
        {
            double E_1 = GetEi(i, params);
            double tau_1 = GetTaui(i, params);
            double denominator = (1+omega2*tau_1*tau_1);
            DEloss(1+i) = omega*tau_1/denominator;
        }
        return;
    }
};

class NonlinearViscoelasticModel : public ViscoelasticModel
{

public:
    NonlinearViscoelasticModel() : ViscoelasticModel() { ; }
    NonlinearViscoelasticModel(int n_p) : ViscoelasticModel(n_p) { ; }
    NonlinearViscoelasticModel(const NonlinearViscoelasticModel&) = delete;
    NonlinearViscoelasticModel& operator=(const NonlinearViscoelasticModel&) = delete;
    ~NonlinearViscoelasticModel() { ; }

    // double GetEinfty(const Vector<double>& params) override { return params(0); }
    // double GetEi(int i, const Vector<double>& params) override { return params(2*i+1); }
    // double GetTaui(int i, const Vector<double>& params) override { return params(2*i+2); }
    
    // double GetDEinfty(const Vector<double>& params) { return 1.; }
    // double GetDEi(int i, const Vector<double>& params) { return 1.; }
    // double GetDTaui(int i, const Vector<double>& params) { return 1.; }
    
    double GetEinfty(const Vector<double>& params) override { return abs(params(0)); }
    double GetEi(int i, const Vector<double>& params) override { return abs(params(2*i+1)); }
    double GetTaui(int i, const Vector<double>& params) override { return abs(params(2*i+2)); }
    
    double GetDEinfty(const Vector<double>& params) { return params(0)/abs(params(0)); }
    double GetDEi(int i, const Vector<double>& params) { return params(2*i+1)/abs(params(2*i+1)); }
    double GetDTaui(int i, const Vector<double>& params) { return params(2*i+2)/abs(params(2*i+2));}

    void evalDEstorage(const Vector<double>& params, double omega, Vector<double>& DEs) override 
    {
        DEs.SetSize(params.Size());
        double Estorage = GetEinfty(params);
        double DEstorage = GetDEinfty(params);
        double omega2 = omega*omega;
        DEs(0) = DEstorage;
        for (int i=0; i<n_param; i++)
        {
            double E_1 = GetEi(i, params);
            double tau_1 = GetTaui(i, params);
            double DE_1 = GetDEi(i, params);
            double Dtau_1 = GetDTaui(i, params);
            double denominator = (1+omega2*tau_1*tau_1);
            DEs(1+2*i) = DE_1 * omega2*tau_1*tau_1/denominator;
            DEs(2+2*i) = E_1 * Dtau_1 * (2*omega2*tau_1)/denominator/denominator;
        }
        return;
    }

        
    void evalDEloss(const Vector<double>& params, double omega, Vector<double>& DEloss) override 
    {
        DEloss.SetSize(params.Size());
        DEloss(0) = 0;
        double omega2 = omega*omega;
        for (int i=0; i<n_param; i++)
        {
            double E_1 = GetEi(i, params);
            double tau_1 = GetTaui(i, params);
            double DE_1 = GetDEi(i, params);
            double Dtau_1 = GetDTaui(i, params);
            double denominator = (1+omega2*tau_1*tau_1);
            DEloss(1+2*i) = DE_1 * omega*tau_1/denominator;
            DEloss(2+2*i) = E_1 *Dtau_1 * omega*(1-omega2*tau_1*tau_1)/denominator/denominator;
        }
        return;
    }
};