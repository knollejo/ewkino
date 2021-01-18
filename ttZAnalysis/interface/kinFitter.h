#include <TGraphErrors.h>
#include <TLorentzVector.h>
#include <TMatrixD.h>
#include <TNamed.h>

class KinFitter {
public:

    struct Particle {
        TLorentzVector vec;
        TMatrixD cov;
        TNamed* ref;
        Particle() : cov(3,3), ref(nullptr) {}
    };

    virtual double JetEtUnc(double, double, double) = 0;
    virtual double JetEtaUnc(double, double, double) = 0;
    virtual double JetPhiUnc(double, double, double) = 0;
    Particle* Jet(TLorentzVector);

    virtual double BjetEtUnc(double et, double eta, double phi) { return JetEtUnc(et, eta, phi); };
    virtual double BjetEtaUnc(double et, double eta, double phi) { return JetEtaUnc(et, eta, phi); };
    virtual double BjetPhiUnc(double et, double eta, double phi) { return JetPhiUnc(et, eta, phi); };
    Particle* Bjet(TLorentzVector);

    virtual double ElectronEtUnc(double, double, double) = 0;
    virtual double ElectronEtaUnc(double, double, double) = 0;
    virtual double ElectronPhiUnc(double, double, double) = 0;
    Particle* Electron(TLorentzVector);

    virtual double MuonEtUnc(double, double, double) = 0;
    virtual double MuonEtaUnc(double, double, double) = 0;
    virtual double MuonPhiUnc(double, double, double) = 0;
    Particle* Muon(TLorentzVector);

    virtual double MissingEtUnc(double, double, double) = 0;
    virtual double MissingEtaUnc(double, double, double) = 0;
    virtual double MissingPhiUnc(double, double, double) = 0;
    Particle* Missing(TLorentzVector);

    struct MassConstraint {
        std::vector<Particle*> particles;
        double mass;
        TNamed* ref;
    };

    virtual double WMass() = 0;
    MassConstraint* WMassConstraint(Particle*, Particle*);

    virtual double TopMass() = 0;
    MassConstraint* TopMassConstraint(Particle*, Particle*, Particle*);

    struct Result {
        int status;
        double f, s;
        int ndf, iterations;
    };
    Result* fit(std::vector<Particle*>, std::vector<MassConstraint*>);

};


class ttZKinFitter : public KinFitter {
public:
    ttZKinFitter(bool, bool);
    virtual ~ttZKinFitter();

    double JetEtUnc(double, double, double);
    double JetEtaUnc(double, double, double);
    double JetPhiUnc(double, double, double);
    double BjetEtUnc(double, double, double);
    double BjetEtaUnc(double, double, double);
    double BjetPhiUnc(double, double, double);
    double ElectronEtUnc(double, double, double);
    double ElectronEtaUnc(double, double, double);
    double ElectronPhiUnc(double, double, double);
    double MuonEtUnc(double, double, double);
    double MuonEtaUnc(double, double, double);
    double MuonPhiUnc(double, double, double);
    double MissingEtUnc(double, double, double);
    double MissingEtaUnc(double, double, double);
    double MissingPhiUnc(double, double, double);

    double WMass();
    double TopMass();
    double BottomMass();

protected:
    TGraphErrors* resolElecEt=nullptr, * resolElecEta=nullptr, * resolElecPhi=nullptr,
                * resolMuonEt=nullptr, * resolMuonEta=nullptr, * resolMuonPhi=nullptr,
                * resolBjetEt=nullptr, * resolBjetEta=nullptr, * resolBjetPhi=nullptr,
                * resolLjetEt=nullptr, * resolLjetEta=nullptr, * resolLjetPhi=nullptr,
                * resolMetEt=nullptr,  * resolMetPhi=nullptr;

};
