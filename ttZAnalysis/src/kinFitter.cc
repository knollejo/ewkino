#include <iostream>
#include <string>

#include <TFile.h>
#include <TMath.h>

#include "../../kinFitter/interface/TFitConstraintM.h"
#include "../../kinFitter/interface/TFitParticleEtEtaPhi.h"
#include "../../kinFitter/interface/TKinFitter.h"

#include "../interface/kinFitter.h"

KinFitter::Particle* KinFitter::Jet(TLorentzVector lv) {
    Particle* p = new Particle;
    p->vec = lv;
    p->cov.Zero();
    p->cov(0,0) = JetEtUnc(lv.Et(), lv.Eta(), lv.Phi());
    p->cov(1,1) = JetEtaUnc(lv.Et(), lv.Eta(), lv.Phi());
    p->cov(2,2) = JetPhiUnc(lv.Et(), lv.Eta(), lv.Phi());
    return p;
}

KinFitter::Particle* KinFitter::Bjet(TLorentzVector lv) {
    Particle* p = new Particle;
    p->vec = lv;
    p->cov.Zero();
    p->cov(0,0) = BjetEtUnc(lv.Et(), lv.Eta(), lv.Phi());
    p->cov(1,1) = BjetEtaUnc(lv.Et(), lv.Eta(), lv.Phi());
    p->cov(2,2) = BjetPhiUnc(lv.Et(), lv.Eta(), lv.Phi());
    return p;
}

KinFitter::Particle* KinFitter::Electron(TLorentzVector lv) {
    Particle* p = new Particle;
    p->vec = lv;
    p->cov.Zero();
    p->cov(0,0) = ElectronEtUnc(lv.Et(), lv.Eta(), lv.Phi());
    p->cov(1,1) = ElectronEtaUnc(lv.Et(), lv.Eta(), lv.Phi());
    p->cov(2,2) = ElectronPhiUnc(lv.Et(), lv.Eta(), lv.Phi());
    return p;
}

KinFitter::Particle* KinFitter::Muon(TLorentzVector lv) {
    Particle* p = new Particle;
    p->vec = lv;
    p->cov.Zero();
    p->cov(0,0) = MuonEtUnc(lv.Et(), lv.Eta(), lv.Phi());
    p->cov(1,1) = MuonEtaUnc(lv.Et(), lv.Eta(), lv.Phi());
    p->cov(2,2) = MuonPhiUnc(lv.Et(), lv.Eta(), lv.Phi());
    return p;
}

KinFitter::Particle* KinFitter::Missing(TLorentzVector lv) {
    Particle* p = new Particle;
    p->vec = lv;
    p->cov.Zero();
    p->cov(0,0) = MissingEtUnc(lv.Et(), lv.Eta(), lv.Phi());
    p->cov(1,1) = MissingEtaUnc(lv.Et(), lv.Eta(), lv.Phi());
    p->cov(2,2) = MissingPhiUnc(lv.Et(), lv.Eta(), lv.Phi());
    return p;
}

KinFitter::MassConstraint* KinFitter::WMassConstraint(Particle* p1, Particle* p2) {
    MassConstraint* c = new MassConstraint;
    c->particles.push_back(p1);
    c->particles.push_back(p2);
    c->mass = WMass();
    return c;
}

KinFitter::MassConstraint* KinFitter::TopMassConstraint(Particle* p1, Particle* p2) {
    MassConstraint* c = new MassConstraint;
    c->particles.push_back(p1);
    c->particles.push_back(p2);
    c->mass = TopMass();
    return c;
}

KinFitter::MassConstraint* KinFitter::TopMassConstraint(Particle* p1, Particle* p2, Particle* p3) {
    MassConstraint* c = new MassConstraint;
    c->particles.push_back(p1);
    c->particles.push_back(p2);
    c->particles.push_back(p3);
    c->mass = TopMass();
    return c;
}

KinFitter::Result* KinFitter::fit(std::vector<Particle*> particles, std::vector<MassConstraint*> constraints) {
    // initialize TKinFitter
    TKinFitter fitter;
    fitter.setMaxNbIter(30);
    fitter.setMaxDeltaS(0.01);
    fitter.setMaxF(0.1);
    fitter.setVerbosity(verbosity);

    // create particles
    unsigned int nP = 0;
    for(std::vector<Particle*>::iterator particle=particles.begin(); particle!=particles.end(); ++particle) {
        (*particle)->ref = new TFitParticleEtEtaPhi(
            ("Particle"+std::to_string(++nP)).c_str(), "",
            &(*particle)->vec, &(*particle)->cov
        );
        fitter.addMeasParticle((TFitParticleEtEtaPhi*)(*particle)->ref);
    }

    // create constraints
    unsigned int nC = 0;
    for(std::vector<MassConstraint*>::iterator constraint=constraints.begin(); constraint!=constraints.end(); ++constraint) {
        TFitConstraintM* c = new TFitConstraintM();
        (*constraint)->ref = c;
        c->SetName(("MassConstraint"+std::to_string(++nC)).c_str());
        for(std::vector<Particle*>::iterator particle=(*constraint)->particles.begin(); particle!=(*constraint)->particles.end(); ++particle)
            c->addParticle1((TFitParticleEtEtaPhi*)(*particle)->ref);
        c->setMassConstraint((*constraint)->mass);
        fitter.addConstraint(c);
    }

    // perform the fit
    std::cerr.setstate(std::ios_base::failbit);
    fitter.fit();
    std::cerr.clear();

    // save the result
    Result* r = new Result;
    r->status = fitter.getStatus();
    r->f = fitter.getF();
    r->s = fitter.getS();
    r->ndf = fitter.getNDF();
    r->iterations = fitter.getNbIter();

    // store the fitted particles & free memory
    for(std::vector<Particle*>::iterator particle=particles.begin(); particle!=particles.end(); ++particle) {
        (*particle)->vec = *((TFitParticleEtEtaPhi*)(*particle)->ref)->getCurr4Vec();
        delete (*particle)->ref;
    }

    // free memory from the constraints
    for(std::vector<MassConstraint*>::iterator constraint=constraints.begin(); constraint!=constraints.end(); ++constraint) {
        delete (*constraint)->ref;
    }

    return r;

}

unsigned int KinFitter::findNeutrinoSolutions(TLorentzVector lep, TLorentzVector met, double& first, double& second) {
    const double a = TMath::Sq(lep.E())-TMath::Sq(lep.Z());
    const double ksq = 0.5*(TMath::Sq(WMass())-TMath::Sq(lep.M()))+lep.X()*met.X()+lep.Y()*met.Y();
    const double b = -2.0*ksq*lep.Z();
    const double c = TMath::Sq(lep.E()*met.X())+TMath::Sq(lep.E()*met.Y())-TMath::Sq(ksq);
    const double det = TMath::Sq(b)-4.0*a*c;
    first = 0.0, second = 0.0;
    if(a==0.0 || det<0.0) return 0;
    else if(det==0.0) {
        first = -0.5*b/a;
        return 1;
    } else {
        first = 0.5*(-1.0*b+TMath::Sqrt(det))/a;
        second = 0.5*(-1.0*b-TMath::Sqrt(det))/a;
        return 2;
    }

}

namespace ttZKinFitterHelpers {

    const int nBinsEt = 16;
    const float binsEt[nBinsEt+1] = { 0.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0, 90.0, 100.0, 120.0, 140.0, 160.0, 180.0, 200.0, 250.0 };
    unsigned int getEtBin(double et) {
        unsigned int iEt;
        for(iEt=0; iEt<nBinsEt; iEt++)
            if(et>=binsEt[iEt] && et<binsEt[iEt+1])
                break;
        return (iEt<nBinsEt-1 ? iEt : nBinsEt-1);
    }

    const int nBinsEtaForEt = 6;
    const float binsEtaForEt[nBinsEtaForEt+1] = { 0.0, 0.435, 0.783, 1.131, 1.479, 1.653, 2.5 };
    unsigned int getEtaBinForEt(double eta) {
        unsigned int iEta;
        eta = (eta<0 ? -1.0*eta : eta);
        for(iEta=0; iEta<nBinsEtaForEt; iEta++)
            if(eta>=binsEtaForEt[iEta] && eta<binsEtaForEt[iEta+1])
                break;
        return (iEta<nBinsEtaForEt-1 ? iEta : nBinsEtaForEt-1);
    }

    const int nBinsEtaForEtaPhi = 2;
    const float binsEtaForEtaPhi[nBinsEtaForEtaPhi+1] = { 0.0, 1.479, 2.5 };
    unsigned int getEtaBinForEtaPhi(double eta) {
        unsigned int iEta;
        eta = (eta<0 ? -1.0*eta : eta);
        for(iEta=0; iEta<nBinsEtaForEtaPhi; iEta++)
            if(eta>=binsEtaForEt[iEta] && eta<binsEtaForEt[iEta+1])
                break;
        return (iEta<nBinsEtaForEtaPhi-1 ? iEta : nBinsEtaForEtaPhi-1);
    }

}

ttZKinFitter::ttZKinFitter(bool is2016, bool is2017, int v) : KinFitter(v) {
    std::string year = is2016 ? "2016" : is2017 ? "2017" : "2018";
    TFile* tfile = TFile::Open(("resolutions/resolutions_"+year+"_V05_v3.root").c_str());
    resolElecEt = new TGraphErrors[ttZKinFitterHelpers::nBinsEtaForEt];
    resolElecEta = new TGraphErrors[ttZKinFitterHelpers::nBinsEtaForEtaPhi];
    resolElecPhi = new TGraphErrors[ttZKinFitterHelpers::nBinsEtaForEtaPhi];
    resolMuonEt = new TGraphErrors[ttZKinFitterHelpers::nBinsEtaForEt];
    resolMuonEta = new TGraphErrors[ttZKinFitterHelpers::nBinsEtaForEtaPhi];
    resolMuonPhi = new TGraphErrors[ttZKinFitterHelpers::nBinsEtaForEtaPhi];
    resolBjetEt = new TGraphErrors[ttZKinFitterHelpers::nBinsEtaForEt];
    resolBjetEta = new TGraphErrors[ttZKinFitterHelpers::nBinsEtaForEtaPhi];
    resolBjetPhi = new TGraphErrors[ttZKinFitterHelpers::nBinsEtaForEtaPhi];
    resolLjetEt = new TGraphErrors[ttZKinFitterHelpers::nBinsEtaForEt];
    resolLjetEta = new TGraphErrors[ttZKinFitterHelpers::nBinsEtaForEtaPhi];
    resolLjetPhi = new TGraphErrors[ttZKinFitterHelpers::nBinsEtaForEtaPhi];
    for(unsigned int iEta=0; iEta<ttZKinFitterHelpers::nBinsEtaForEt; iEta++) {
        resolElecEt[iEta] = *(TGraphErrors*)tfile->Get(("rmses_"+year+"_ttbar_eta"+std::to_string(iEta)+"_resolElecEt").c_str());
        resolMuonEt[iEta] = *(TGraphErrors*)tfile->Get(("rmses_"+year+"_ttbar_eta"+std::to_string(iEta)+"_resolMuonEt").c_str());
        resolBjetEt[iEta] = *(TGraphErrors*)tfile->Get(("rmses_"+year+"_ttbar_eta"+std::to_string(iEta)+"_resolBjetEt").c_str());
        resolLjetEt[iEta] = *(TGraphErrors*)tfile->Get(("rmses_"+year+"_ttbar_eta"+std::to_string(iEta)+"_resolLjetEt").c_str());
    }
    for(unsigned int iEta=0; iEta<ttZKinFitterHelpers::nBinsEtaForEtaPhi; iEta++) {
        resolElecEta[iEta] = *(TGraphErrors*)tfile->Get(("rmses_"+year+"_ttbar_eta"+std::to_string(iEta)+"_resolElecEta").c_str());
        resolElecPhi[iEta] = *(TGraphErrors*)tfile->Get(("rmses_"+year+"_ttbar_eta"+std::to_string(iEta)+"_resolElecPhi").c_str());
        resolMuonEta[iEta] = *(TGraphErrors*)tfile->Get(("rmses_"+year+"_ttbar_eta"+std::to_string(iEta)+"_resolMuonEta").c_str());
        resolMuonPhi[iEta] = *(TGraphErrors*)tfile->Get(("rmses_"+year+"_ttbar_eta"+std::to_string(iEta)+"_resolMuonPhi").c_str());
        resolBjetEta[iEta] = *(TGraphErrors*)tfile->Get(("rmses_"+year+"_ttbar_eta"+std::to_string(iEta)+"_resolBjetEta").c_str());
        resolBjetPhi[iEta] = *(TGraphErrors*)tfile->Get(("rmses_"+year+"_ttbar_eta"+std::to_string(iEta)+"_resolBjetPhi").c_str());
        resolLjetEta[iEta] = *(TGraphErrors*)tfile->Get(("rmses_"+year+"_ttbar_eta"+std::to_string(iEta)+"_resolLjetEta").c_str());
        resolLjetPhi[iEta] = *(TGraphErrors*)tfile->Get(("rmses_"+year+"_ttbar_eta"+std::to_string(iEta)+"_resolLjetPhi").c_str());
    }
    resolMetEt = (TGraphErrors*)tfile->Get(("rmses_"+year+"_ttbar_eta0_resolMetEt").c_str());
    resolMetPhi = (TGraphErrors*)tfile->Get(("rmses_"+year+"_ttbar_eta0_resolMetPhi").c_str());
    tfile->Close();
}

ttZKinFitter::~ttZKinFitter() {
    if(resolElecEt!=nullptr) delete [] resolElecEt;
    if(resolElecEta!=nullptr) delete [] resolElecEta;
    if(resolElecPhi!=nullptr) delete [] resolElecPhi;
    if(resolMuonEt!=nullptr) delete [] resolMuonEt;
    if(resolMuonEta!=nullptr) delete [] resolMuonEta;
    if(resolMuonPhi!=nullptr) delete [] resolMuonPhi;
    if(resolBjetEt!=nullptr) delete [] resolBjetEt;
    if(resolBjetEta!=nullptr) delete [] resolBjetEta;
    if(resolBjetPhi!=nullptr) delete [] resolBjetPhi;
    if(resolLjetEt!=nullptr) delete [] resolLjetEt;
    if(resolLjetEta!=nullptr) delete [] resolLjetEta;
    if(resolLjetPhi!=nullptr) delete [] resolLjetPhi;
    if(resolMetEt!=nullptr) delete resolMetEt;
    if(resolMetPhi!=nullptr) delete resolMetPhi;
}

double ttZKinFitter::JetEtUnc(double et, double eta, double) {
    const unsigned int iEt = ttZKinFitterHelpers::getEtBin(et);
    const unsigned int iEta = ttZKinFitterHelpers::getEtaBinForEt(eta);
    return et*resolLjetEt[iEta].GetY()[iEt];
}

double ttZKinFitter::JetEtaUnc(double et, double eta, double) {
    const unsigned int iEt = ttZKinFitterHelpers::getEtBin(et);
    const unsigned int iEta = ttZKinFitterHelpers::getEtaBinForEtaPhi(eta);
    return resolLjetEta[iEta].GetY()[iEt];
}

double ttZKinFitter::JetPhiUnc(double et, double eta, double) {
    const unsigned int iEt = ttZKinFitterHelpers::getEtBin(et);
    const unsigned int iEta = ttZKinFitterHelpers::getEtaBinForEtaPhi(eta);
    return resolLjetPhi[iEta].GetY()[iEt];
}

double ttZKinFitter::BjetEtUnc(double et, double eta, double) {
    const unsigned int iEt = ttZKinFitterHelpers::getEtBin(et);
    const unsigned int iEta = ttZKinFitterHelpers::getEtaBinForEt(eta);
    return et*resolBjetEt[iEta].GetY()[iEt];
}

double ttZKinFitter::BjetEtaUnc(double et, double eta, double) {
    const unsigned int iEt = ttZKinFitterHelpers::getEtBin(et);
    const unsigned int iEta = ttZKinFitterHelpers::getEtaBinForEtaPhi(eta);
    return resolBjetEta[iEta].GetY()[iEt];
}

double ttZKinFitter::BjetPhiUnc(double et, double eta, double) {
    const unsigned int iEt = ttZKinFitterHelpers::getEtBin(et);
    const unsigned int iEta = ttZKinFitterHelpers::getEtaBinForEtaPhi(eta);
    return resolBjetPhi[iEta].GetY()[iEt];
}

double ttZKinFitter::ElectronEtUnc(double et, double eta, double) {
    const unsigned int iEt = ttZKinFitterHelpers::getEtBin(et);
    const unsigned int iEta = ttZKinFitterHelpers::getEtaBinForEt(eta);
    return et*resolElecEt[iEta].GetY()[iEt];
}

double ttZKinFitter::ElectronEtaUnc(double et, double eta, double) {
    const unsigned int iEt = ttZKinFitterHelpers::getEtBin(et);
    const unsigned int iEta = ttZKinFitterHelpers::getEtaBinForEtaPhi(eta);
    return resolElecEta[iEta].GetY()[iEt];
}

double ttZKinFitter::ElectronPhiUnc(double et, double eta, double) {
    const unsigned int iEt = ttZKinFitterHelpers::getEtBin(et);
    const unsigned int iEta = ttZKinFitterHelpers::getEtaBinForEtaPhi(eta);
    return resolElecPhi[iEta].GetY()[iEt];
}

double ttZKinFitter::MuonEtUnc(double et, double eta, double) {
    const unsigned int iEt = ttZKinFitterHelpers::getEtBin(et);
    const unsigned int iEta = ttZKinFitterHelpers::getEtaBinForEt(eta);
    return et*resolMuonEt[iEta].GetY()[iEt];
}

double ttZKinFitter::MuonEtaUnc(double et, double eta, double) {
    const unsigned int iEt = ttZKinFitterHelpers::getEtBin(et);
    const unsigned int iEta = ttZKinFitterHelpers::getEtaBinForEtaPhi(eta);
    return resolMuonEta[iEta].GetY()[iEt];
}

double ttZKinFitter::MuonPhiUnc(double et, double eta, double) {
    const unsigned int iEt = ttZKinFitterHelpers::getEtBin(et);
    const unsigned int iEta = ttZKinFitterHelpers::getEtaBinForEtaPhi(eta);
    return resolMuonPhi[iEta].GetY()[iEt];
}

double ttZKinFitter::MissingEtUnc(double et, double, double) {
    const unsigned int iEt = ttZKinFitterHelpers::getEtBin(et);
    return resolMetEt->GetY()[iEt];
}

double ttZKinFitter::MissingEtaUnc(double, double, double) {
    return 9999.0;
}

double ttZKinFitter::MissingPhiUnc(double et, double, double) {
    const unsigned int iEt = ttZKinFitterHelpers::getEtBin(et);
    return resolMetPhi->GetY()[iEt];
}

double ttZKinFitter::WMass() {
    return 80.38;
}

double ttZKinFitter::TopMass() {
    return 172.5;
}

double ttZKinFitter::BottomMass() {
    return 4.92;
}
