#include "../interface/ttZHistograms.h"

//include ROOT classes
#include "TMath.h"
#include "TVector2.h"

//include other parts of framework
#include "../../constants/particleMasses.h"

// histogram definitions
const int NBINS_category3l = 4; const double MIN_category3l = -0.5; const double MAX_category3l = 3.5;
const int NBINS_category4l = 6; const double MIN_category4l = -0.5; const double MAX_category4l = 5.5;
const int NBINS_nJets = 6; const double MIN_nJets = 1.5; const double MAX_nJets = 7.5;
const int NBINS_nBjets = 3; const double MIN_nBjets = 0.5; const double MAX_nBjets = 3.5;
const int NBINS_dilepPt = 10; const double MIN_dilepPt = 0.0; const double MAX_dilepPt = 500.0;
const int NBINS_dilepEta = 12; const double MIN_dilepEta = -3.0; const double MAX_dilepEta = 3.0;
const int NBINS_dilepPhi = 10; const double MIN_dilepPhi = -TMath::Pi(); const double MAX_dilepPhi = TMath::Pi();
const int NBINS_dilepMass3l = 10; const double MIN_dilepMass3l = particle::mZ-10.0; const double MAX_dilepMass3l = particle::mZ+10.0;
const int NBINS_dilepMass4l = 10; const double MIN_dilepMass4l = particle::mZ-20.0; const double MAX_dilepMass4l = particle::mZ+20.0;
const int NBINS_missingEt = 11; const double MIN_missingEt = 0.0; const double MAX_missingEt = 220.0;
const int NBINS_missingPhi = 10; const double MIN_missingPhi = -TMath::Pi(); const double MAX_missingPhi = TMath::Pi();
const int NBINS_firstLepPt = 10; const double MIN_firstLepPt = 40.0; const double MAX_firstLepPt = 340.0;
const int NBINS_firstLepEta = 10; const double MIN_firstLepEta = -2.5; const double MAX_firstLepEta = 2.5;
const int NBINS_firstLepPhi = 10; const double MIN_firstLepPhi = -TMath::Pi(); const double MAX_firstLepPhi = -TMath::Pi();
const int NBINS_secondLepPt = 10; const double MIN_secondLepPt = 10.0; const double MAX_secondLepPt = 160.0;
const int NBINS_secondLepEta = 10; const double MIN_secondLepEta = -2.5; const double MAX_secondLepEta = 2.5;
const int NBINS_secondLepPhi = 10; const double MIN_secondLepPhi = -TMath::Pi(); const double MAX_secondLepPhi = -TMath::Pi();
const int NBINS_thirdLepPt = 10; const double MIN_thirdLepPt = 10.0; const double MAX_thirdLepPt = 110.0;
const int NBINS_thirdLepEta = 10; const double MIN_thirdLepEta = -2.5; const double MAX_thirdLepEta = 2.5;
const int NBINS_thirdLepPhi = 10; const double MIN_thirdLepPhi = -TMath::Pi(); const double MAX_thirdLepPhi = -TMath::Pi();
const int NBINS_fourthLepPt = 10; const double MIN_fourthLepPt = 10.0; const double MAX_fourthLepPt = 110.0;
const int NBINS_fourthLepEta = 10; const double MIN_fourthLepEta = -2.5; const double MAX_fourthLepEta = 2.5;
const int NBINS_fourthLepPhi = 10; const double MIN_fourthLepPhi = -TMath::Pi(); const double MAX_fourthLepPhi = -TMath::Pi();
const int NBINS_firstJetPt = 10; const double MIN_firstJetPt = 30.0; const double MAX_firstJetPt = 330.0;
const int NBINS_firstJetEta = 10; const double MIN_firstJetEta = -2.5; const double MAX_firstJetEta = 2.5;
const int NBINS_firstJetPhi = 10; const double MIN_firstJetPhi = -TMath::Pi(); const double MAX_firstJetPhi = -TMath::Pi();
const int NBINS_secondJetPt = 10; const double MIN_secondJetPt = 30.0; const double MAX_secondJetPt = 180.0;
const int NBINS_secondJetEta = 10; const double MIN_secondJetEta = -2.5; const double MAX_secondJetEta = 2.5;
const int NBINS_secondJetPhi = 10; const double MIN_secondJetPhi = -TMath::Pi(); const double MAX_secondJetPhi = -TMath::Pi();
const int NBINS_thirdJetPt = 10; const double MIN_thirdJetPt = 30.0; const double MAX_thirdJetPt = 130.0;
const int NBINS_thirdJetEta = 10; const double MIN_thirdJetEta = -2.5; const double MAX_thirdJetEta = 2.5;
const int NBINS_thirdJetPhi = 10; const double MIN_thirdJetPhi = -TMath::Pi(); const double MAX_thirdJetPhi = -TMath::Pi();
const int NBINS_fourthJetPt = 10; const double MIN_fourthJetPt = 30.0; const double MAX_fourthJetPt = 130.0;
const int NBINS_fourthJetEta = 10; const double MIN_fourthJetEta = -2.5; const double MAX_fourthJetEta = 2.5;
const int NBINS_fourthJetPhi = 10; const double MIN_fourthJetPhi = -TMath::Pi(); const double MAX_fourthJetPhi = -TMath::Pi();
const int NBINS_zbosonPt = 8; const double BINS_zbosonPt[NBINS_zbosonPt+1] = { 0.0, 37.5, 75.0, 112.5, 150.0, 200.0, 250.0, 375.0, 500.0 };
const int NBINS_ttzMass = 8; const double BINS_ttzMass[NBINS_ttzMass+1] = { 436.0, 580.0, 640.0, 720.0, 800.0, 890.0, 1000.0, 1150.0, 1500.0 };
const int NBINS_ttbarMass = 8; const double BINS_ttbarMass[NBINS_ttbarMass+1] = { 345.0, 400.0, 440.0, 485.0, 540.0, 610.0, 700.0, 850.0, 1200.0 };
const int NBINS_topPt = 8; const double BINS_topPt[NBINS_topPt+1] = { 0.0, 60.0, 90.0, 125.0, 160.0, 200.0, 250.0, 330.0, 500.0 };
const int NBINS_deltaPhiTtbar = 8; const double BINS_deltaPhiTtbar[NBINS_deltaPhiTtbar+1] = { 0.0, TMath::Pi()/4, TMath::Pi()/2, 5*TMath::Pi()/8, 3*TMath::Pi()/4, 13*TMath::Pi()/16, 7*TMath::Pi()/8, 15*TMath::Pi()/16, TMath::Pi() };
const int NBINS_deltaPhiTopZ = 8; const double BINS_deltaPhiTopZ[NBINS_deltaPhiTopZ+1] = { 0.0, TMath::Pi()/6, TMath::Pi()/3, TMath::Pi()/2, 2*TMath::Pi()/3, 3*TMath::Pi()/4, 5*TMath::Pi()/6, 11*TMath::Pi()/12, TMath::Pi() };
const int NBINS_deltaRapTtbar = 8; const double BINS_deltaRapTtbar[NBINS_deltaRapTtbar+1] = { 0.0, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 2.1, 3.0 };
const int NBINS_deltaRapTopZ = 8; const double BINS_deltaRapTopZ[NBINS_deltaRapTopZ+1] = { 0.0, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 2.0, 3.0 };
const int NBINS_lepTopMass = 10; const double MIN_lepTopMass = 170.0; const double MAX_lepTopMass = 175.0;
const int NBINS_hadTopMass = 10; const double MIN_hadTopMass = 170.0; const double MAX_hadTopMass = 175.0;
const int NBINS_topLeptonPt = 8; const double MIN_topLeptonPt = 0.0; const double MAX_topLeptonPt = 250.0;
const int NBINS_topLeptonsMass = 8; const double MIN_topLeptonsMass = 0.0; const double MAX_topLeptonsMass = 350.0;
const int NBINS_fourLeptonsMass = 8; const double MIN_fourLeptonsMass = 100.0; const double MAX_fourLeptonsMass = 600.0;
const int NBINS_deltaPhiTopLeptons = 8; const double MIN_deltaPhiTopLeptons = 0.0; const double MAX_deltaPhiTopLeptons = TMath::Pi();
const int NBINS_deltaPhiTopLeptonZ = 8; const double MIN_deltaPhiTopLeptonZ = 0.0; const double MAX_deltaPhiTopLeptonZ = TMath::Pi();
const int NBINS_deltaRapTopLeptons = 8; const double BINS_deltaRapTopLeptons[NBINS_deltaRapTopLeptons+1] = { 0.0, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 2.0, 3.0 };
const int NBINS_deltaRapTopLeptonZ = 8; const double BINS_deltaRapTopLeptonZ[NBINS_deltaRapTopLeptonZ+1] = { 0.0, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 2.0, 3.0 };

ttZHistograms::ttZHistograms(int nVar, int nSel, int nSys) :
    nVariants(nVar), nSelections(nSel), nSystematics(nSys),
    nHistograms(nVariants*nSelections*(1+nSystematics)),
    hists_category3l(new TH1F*[nHistograms]),
    hists_category4l(new TH1F*[nHistograms]),
    hists_nJets(new TH1F*[nHistograms]),
    hists_nBjets(new TH1F*[nHistograms]),
    hists_dilepPt(new TH1F*[nHistograms]),
    hists_dilepEta(new TH1F*[nHistograms]),
    hists_dilepPhi(new TH1F*[nHistograms]),
    hists_dilepMass3l(new TH1F*[nHistograms]),
    hists_dilepMass4l(new TH1F*[nHistograms]),
    hists_missingEt(new TH1F*[nHistograms]),
    hists_missingPhi(new TH1F*[nHistograms]),
    hists_firstLepPt(new TH1F*[nHistograms]),
    hists_firstLepEta(new TH1F*[nHistograms]),
    hists_firstLepPhi(new TH1F*[nHistograms]),
    hists_secondLepPt(new TH1F*[nHistograms]),
    hists_secondLepEta(new TH1F*[nHistograms]),
    hists_secondLepPhi(new TH1F*[nHistograms]),
    hists_thirdLepPt(new TH1F*[nHistograms]),
    hists_thirdLepEta(new TH1F*[nHistograms]),
    hists_thirdLepPhi(new TH1F*[nHistograms]),
    hists_fourthLepPt(new TH1F*[nHistograms]),
    hists_fourthLepEta(new TH1F*[nHistograms]),
    hists_fourthLepPhi(new TH1F*[nHistograms]),
    hists_firstJetPt(new TH1F*[nHistograms]),
    hists_firstJetEta(new TH1F*[nHistograms]),
    hists_firstJetPhi(new TH1F*[nHistograms]),
    hists_secondJetPt(new TH1F*[nHistograms]),
    hists_secondJetEta(new TH1F*[nHistograms]),
    hists_secondJetPhi(new TH1F*[nHistograms]),
    hists_thirdJetPt(new TH1F*[nHistograms]),
    hists_thirdJetEta(new TH1F*[nHistograms]),
    hists_thirdJetPhi(new TH1F*[nHistograms]),
    hists_fourthJetPt(new TH1F*[nHistograms]),
    hists_fourthJetEta(new TH1F*[nHistograms]),
    hists_fourthJetPhi(new TH1F*[nHistograms]),
    hists_zbosonPt(new TH1F*[nHistograms]),
    hists_ttzMass(new TH1F*[nHistograms]),
    hists_ttbarMass(new TH1F*[nHistograms]),
    hists_topPt(new TH1F*[nHistograms]),
    hists_deltaPhiTtbar(new TH1F*[nHistograms]),
    hists_deltaPhiTopZ(new TH1F*[nHistograms]),
    hists_deltaRapTtbar(new TH1F*[nHistograms]),
    hists_deltaRapTopZ(new TH1F*[nHistograms]),
    hists_lepTopMass(new TH1F*[nHistograms]),
    hists_hadTopMass(new TH1F*[nHistograms]),
    hists_topLeptonPt(new TH1F*[nHistograms]),
    hists_topLeptonsMass(new TH1F*[nHistograms]),
    hists_fourLeptonsMass(new TH1F*[nHistograms]),
    hists_deltaPhiTopLeptons(new TH1F*[nHistograms]),
    hists_deltaPhiTopLeptonZ(new TH1F*[nHistograms]),
    hists_deltaRapTopLeptons(new TH1F*[nHistograms]),
    hists_deltaRapTopLeptonZ(new TH1F*[nHistograms])
{
    for(int i=0; i<nHistograms; i++) {
        const std::string sI = std::to_string(i);
        hists_category3l[i] = new TH1F(("hist_category3l_var"+sI).c_str(), "", NBINS_category3l, MIN_category3l, MAX_category3l);
        hists_category4l[i] = new TH1F(("hist_category4l_var"+sI).c_str(), "", NBINS_category4l, MIN_category4l, MAX_category4l);
        hists_nJets[i] = new TH1F(("hist_nJets_var"+sI).c_str(), "", NBINS_nJets, MIN_nJets, MAX_nJets);
        hists_nBjets[i] = new TH1F(("hist_nBjets_var"+sI).c_str(), "", NBINS_nBjets, MIN_nBjets, MAX_nBjets);
        hists_dilepPt[i] = new TH1F(("hist_dilepPt_var"+sI).c_str(), "", NBINS_dilepPt, MIN_dilepPt, MAX_dilepPt);
        hists_dilepEta[i] = new TH1F(("hist_dilepEta_var"+sI).c_str(), "", NBINS_dilepEta, MIN_dilepEta, MAX_dilepEta);
        hists_dilepPhi[i] = new TH1F(("hist_dilepPhi_var"+sI).c_str(), "", NBINS_dilepPhi, MIN_dilepPhi, MAX_dilepPhi);
        hists_dilepMass3l[i] = new TH1F(("hist_dilepMass3l_var"+sI).c_str(), "", NBINS_dilepMass3l, MIN_dilepMass3l, MAX_dilepMass3l);
        hists_dilepMass4l[i] = new TH1F(("hist_dilepMass4l_var"+sI).c_str(), "", NBINS_dilepMass4l, MIN_dilepMass4l, MAX_dilepMass4l);
        hists_missingEt[i] = new TH1F(("hist_missingEt_var"+sI).c_str(), "", NBINS_missingEt, MIN_missingEt, MAX_missingEt);
        hists_missingPhi[i] = new TH1F(("hist_missingPhi_var"+sI).c_str(), "", NBINS_missingPhi, MIN_missingPhi, MAX_missingPhi);
        hists_firstLepPt[i] = new TH1F(("hist_firstLepPt_var"+sI).c_str(), "", NBINS_firstLepPt, MIN_firstLepPt, MAX_firstLepPt);
        hists_firstLepEta[i] = new TH1F(("hist_firstLepEta_var"+sI).c_str(), "", NBINS_firstLepEta, MIN_firstLepEta, MAX_firstLepEta);
        hists_firstLepPhi[i] = new TH1F(("hist_firstLepPhi_var"+sI).c_str(), "", NBINS_firstLepPhi, MIN_firstLepPhi, MAX_firstLepPhi);
        hists_secondLepPt[i] = new TH1F(("hist_secondLepPt_var"+sI).c_str(), "", NBINS_secondLepPt, MIN_secondLepPt, MAX_secondLepPt);
        hists_secondLepEta[i] = new TH1F(("hist_secondLepEta_var"+sI).c_str(), "", NBINS_secondLepEta, MIN_secondLepEta, MAX_secondLepEta);
        hists_secondLepPhi[i] = new TH1F(("hist_secondLepPhi_var"+sI).c_str(), "", NBINS_secondLepPhi, MIN_secondLepPhi, MAX_secondLepPhi);
        hists_thirdLepPt[i] = new TH1F(("hist_thirdLepPt_var"+sI).c_str(), "", NBINS_thirdLepPt, MIN_thirdLepPt, MAX_thirdLepPt);
        hists_thirdLepEta[i] = new TH1F(("hist_thirdLepEta_var"+sI).c_str(), "", NBINS_thirdLepEta, MIN_thirdLepEta, MAX_thirdLepEta);
        hists_thirdLepPhi[i] = new TH1F(("hist_thirdLepPhi_var"+sI).c_str(), "", NBINS_thirdLepPhi, MIN_thirdLepPhi, MAX_thirdLepPhi);
        hists_fourthLepPt[i] = new TH1F(("hist_fourthLepPt_var"+sI).c_str(), "", NBINS_fourthLepPt, MIN_fourthLepPt, MAX_fourthLepPt);
        hists_fourthLepEta[i] = new TH1F(("hist_fourthLepEta_var"+sI).c_str(), "", NBINS_fourthLepEta, MIN_fourthLepEta, MAX_fourthLepEta);
        hists_fourthLepPhi[i] = new TH1F(("hist_fourthLepPhi_var"+sI).c_str(), "", NBINS_fourthLepPhi, MIN_fourthLepPhi, MAX_fourthLepPhi);
        hists_firstJetPt[i] = new TH1F(("hist_firstJetPt_var"+sI).c_str(), "", NBINS_firstJetPt, MIN_firstJetPt, MAX_firstJetPt);
        hists_firstJetEta[i] = new TH1F(("hist_firstJetEta_var"+sI).c_str(), "", NBINS_firstJetEta, MIN_firstJetEta, MAX_firstJetEta);
        hists_firstJetPhi[i] = new TH1F(("hist_firstJetPhi_var"+sI).c_str(), "", NBINS_firstJetPhi, MIN_firstJetPhi, MAX_firstJetPhi);
        hists_secondJetPt[i] = new TH1F(("hist_secondJetPt_var"+sI).c_str(), "", NBINS_secondJetPt, MIN_secondJetPt, MAX_secondJetPt);
        hists_secondJetEta[i] = new TH1F(("hist_secondJetEta_var"+sI).c_str(), "", NBINS_secondJetEta, MIN_secondJetEta, MAX_secondJetEta);
        hists_secondJetPhi[i] = new TH1F(("hist_secondJetPhi_var"+sI).c_str(), "", NBINS_secondJetPhi, MIN_secondJetPhi, MAX_secondJetPhi);
        hists_thirdJetPt[i] = new TH1F(("hist_thirdJetPt_var"+sI).c_str(), "", NBINS_thirdJetPt, MIN_thirdJetPt, MAX_thirdJetPt);
        hists_thirdJetEta[i] = new TH1F(("hist_thirdJetEta_var"+sI).c_str(), "", NBINS_thirdJetEta, MIN_thirdJetEta, MAX_thirdJetEta);
        hists_thirdJetPhi[i] = new TH1F(("hist_thirdJetPhi_var"+sI).c_str(), "", NBINS_thirdJetPhi, MIN_thirdJetPhi, MAX_thirdJetPhi);
        hists_fourthJetPt[i] = new TH1F(("hist_fourthJetPt_var"+sI).c_str(), "", NBINS_fourthJetPt, MIN_fourthJetPt, MAX_fourthJetPt);
        hists_fourthJetEta[i] = new TH1F(("hist_fourthJetEta_var"+sI).c_str(), "", NBINS_fourthJetEta, MIN_fourthJetEta, MAX_fourthJetEta);
        hists_fourthJetPhi[i] = new TH1F(("hist_fourthJetPhi_var"+sI).c_str(), "", NBINS_fourthJetPhi, MIN_fourthJetPhi, MAX_fourthJetPhi);
        hists_zbosonPt[i] = new TH1F(("hist_zbosonPt_var"+sI).c_str(), "", NBINS_zbosonPt, BINS_zbosonPt);
        hists_ttzMass[i] = new TH1F(("hist_ttzMass_var"+sI).c_str(), "", NBINS_ttzMass, BINS_ttzMass);
        hists_ttbarMass[i] = new TH1F(("hist_ttbarMass_var"+sI).c_str(), "", NBINS_ttbarMass, BINS_ttbarMass);
        hists_topPt[i] = new TH1F(("hist_topPt_var"+sI).c_str(), "", NBINS_topPt, BINS_topPt);
        hists_deltaPhiTtbar[i] = new TH1F(("hist_deltaPhiTtbar_var"+sI).c_str(), "", NBINS_deltaPhiTtbar, BINS_deltaPhiTtbar);
        hists_deltaPhiTopZ[i] = new TH1F(("hist_deltaPhiTopZ_var"+sI).c_str(), "", NBINS_deltaPhiTopZ, BINS_deltaPhiTopZ);
        hists_deltaRapTtbar[i] = new TH1F(("hist_deltaRapTtbar_var"+sI).c_str(), "", NBINS_deltaRapTtbar, BINS_deltaRapTtbar);
        hists_deltaRapTopZ[i] = new TH1F(("hist_deltaRapTopZ_var"+sI).c_str(), "", NBINS_deltaRapTopZ, BINS_deltaRapTopZ);
        hists_lepTopMass[i] = new TH1F(("hist_lepTopMass_var"+sI).c_str(), "", NBINS_lepTopMass, MIN_lepTopMass, MAX_lepTopMass);
        hists_hadTopMass[i] = new TH1F(("hist_hadTopMass_var"+sI).c_str(), "", NBINS_hadTopMass, MIN_hadTopMass, MAX_hadTopMass);
        hists_topLeptonPt[i] = new TH1F(("hist_topLeptonPt_var"+sI).c_str(), "", NBINS_topLeptonPt, MIN_topLeptonPt, MAX_topLeptonPt);
        hists_topLeptonsMass[i] = new TH1F(("hist_topLeptonsMass_var"+sI).c_str(), "", NBINS_topLeptonsMass, MIN_topLeptonsMass, MAX_topLeptonsMass);
        hists_fourLeptonsMass[i] = new TH1F(("hist_fourLeptonsMass_var"+sI).c_str(), "", NBINS_fourLeptonsMass, MIN_fourLeptonsMass, MAX_fourLeptonsMass);
        hists_deltaPhiTopLeptons[i] = new TH1F(("hist_deltaPhiTopLeptons_var"+sI).c_str(), "", NBINS_deltaPhiTopLeptons, MIN_deltaPhiTopLeptons, MAX_deltaPhiTopLeptons);
        hists_deltaPhiTopLeptonZ[i] = new TH1F(("hist_deltaPhiTopLeptonZ_var"+sI).c_str(), "", NBINS_deltaPhiTopLeptonZ, MIN_deltaPhiTopLeptonZ, MAX_deltaPhiTopLeptonZ);
        hists_deltaRapTopLeptons[i] = new TH1F(("hist_deltaRapTopLeptons_var"+sI).c_str(), "", NBINS_deltaRapTopLeptons, BINS_deltaRapTopLeptons);
        hists_deltaRapTopLeptonZ[i] = new TH1F(("hist_deltaRapTopLeptonZ_var"+sI).c_str(), "", NBINS_deltaRapTopLeptonZ, BINS_deltaRapTopLeptonZ);
    }
}

ttZHistograms::~ttZHistograms() {
    for(int i=0; i<nHistograms; i++) {
        delete hists_category3l[i];
        delete hists_category4l[i];
        delete hists_nJets[i];
        delete hists_nBjets[i];
        delete hists_dilepPt[i];
        delete hists_dilepEta[i];
        delete hists_dilepPhi[i];
        delete hists_dilepMass3l[i];
        delete hists_dilepMass4l[i];
        delete hists_missingEt[i];
        delete hists_missingPhi[i];
        delete hists_firstLepPt[i];
        delete hists_firstLepEta[i];
        delete hists_firstLepPhi[i];
        delete hists_secondLepPt[i];
        delete hists_secondLepEta[i];
        delete hists_secondLepPhi[i];
        delete hists_thirdLepPt[i];
        delete hists_thirdLepEta[i];
        delete hists_thirdLepPhi[i];
        delete hists_fourthLepPt[i];
        delete hists_fourthLepEta[i];
        delete hists_fourthLepPhi[i];
        delete hists_firstJetPt[i];
        delete hists_firstJetEta[i];
        delete hists_firstJetPhi[i];
        delete hists_secondJetPt[i];
        delete hists_secondJetEta[i];
        delete hists_secondJetPhi[i];
        delete hists_thirdJetPt[i];
        delete hists_thirdJetEta[i];
        delete hists_thirdJetPhi[i];
        delete hists_fourthJetPt[i];
        delete hists_fourthJetEta[i];
        delete hists_fourthJetPhi[i];
        delete hists_zbosonPt[i];
        delete hists_ttzMass[i];
        delete hists_ttbarMass[i];
        delete hists_topPt[i];
        delete hists_deltaPhiTtbar[i];
        delete hists_deltaPhiTopZ[i];
        delete hists_deltaRapTtbar[i];
        delete hists_deltaRapTopZ[i];
        delete hists_lepTopMass[i];
        delete hists_hadTopMass[i];
        delete hists_topLeptonPt[i];
        delete hists_topLeptonsMass[i];
        delete hists_fourLeptonsMass[i];
        delete hists_deltaPhiTopLeptons[i];
        delete hists_deltaPhiTopLeptonZ[i];
        delete hists_deltaRapTopLeptons[i];
        delete hists_deltaRapTopLeptonZ[i];
    }
    delete [] hists_category3l;
    delete [] hists_category4l;
    delete [] hists_nJets;
    delete [] hists_nBjets;
    delete [] hists_dilepPt;
    delete [] hists_dilepEta;
    delete [] hists_dilepPhi;
    delete [] hists_dilepMass3l;
    delete [] hists_dilepMass4l;
    delete [] hists_missingEt;
    delete [] hists_missingPhi;
    delete [] hists_firstLepPt;
    delete [] hists_firstLepEta;
    delete [] hists_firstLepPhi;
    delete [] hists_secondLepPt;
    delete [] hists_secondLepEta;
    delete [] hists_secondLepPhi;
    delete [] hists_thirdLepPt;
    delete [] hists_thirdLepEta;
    delete [] hists_thirdLepPhi;
    delete [] hists_fourthLepPt;
    delete [] hists_fourthLepEta;
    delete [] hists_fourthLepPhi;
    delete [] hists_firstJetPt;
    delete [] hists_firstJetEta;
    delete [] hists_firstJetPhi;
    delete [] hists_secondJetPt;
    delete [] hists_secondJetEta;
    delete [] hists_secondJetPhi;
    delete [] hists_thirdJetPt;
    delete [] hists_thirdJetEta;
    delete [] hists_thirdJetPhi;
    delete [] hists_fourthJetPt;
    delete [] hists_fourthJetEta;
    delete [] hists_fourthJetPhi;
    delete [] hists_zbosonPt;
    delete [] hists_ttzMass;
    delete [] hists_ttbarMass;
    delete [] hists_topPt;
    delete [] hists_deltaPhiTtbar;
    delete [] hists_deltaPhiTopZ;
    delete [] hists_deltaRapTtbar;
    delete [] hists_deltaRapTopZ;
    delete [] hists_lepTopMass;
    delete [] hists_hadTopMass;
    delete [] hists_topLeptonPt;
    delete [] hists_topLeptonsMass;
    delete [] hists_fourLeptonsMass;
    delete [] hists_deltaPhiTopLeptons;
    delete [] hists_deltaPhiTopLeptonZ;
    delete [] hists_deltaRapTopLeptons;
    delete [] hists_deltaRapTopLeptonZ;
}

#define IN_RANGE(VAL, MIN, MAX) \
    (VAL<-900.0 ? -999.0 : TMath::Min(MAX-0.001*(MAX-MIN), TMath::Max(MIN+0.001*(MAX-MIN), VAL)))
#define PHI_IN_RANGE(VAL) \
    (VAL<-900.0 || isnan(VAL) ? -999.0 : TMath::Min(0.999*TMath::Pi(), TMath::Max(-0.999*TMath::Pi(), TVector2::Phi_mpi_pi(VAL))))

void ttZHistograms::SetValues(ttZ::leptonVariables leptonVars, ttZ::jetVariables jetVars, ttZ::reconstructedVariables reconstructedVars, ttZ::fourLeptonVariables fourLeptonVars) {
    value_category3l = leptonVars.isThreeLeptons ? leptonVars.category : -1;
    value_category4l = leptonVars.isThreeLeptons ? -1 : leptonVars.category;
    value_nJets = jetVars.nJets>MAX_nJets ? MAX_nJets-0.5 : jetVars.nJets;
    value_nBjets = jetVars.nBjets>MAX_nBjets ? MAX_nBjets-0.5 : jetVars.nBjets;
    value_dilepPt = IN_RANGE(leptonVars.dilepPt, MIN_dilepPt, MAX_dilepPt);
    value_dilepEta = IN_RANGE(leptonVars.dilepEta, MIN_dilepEta, MAX_dilepEta);
    value_dilepPhi = PHI_IN_RANGE(leptonVars.dilepPhi);
    value_dilepMass3l = IN_RANGE(leptonVars.dilepMass, MIN_dilepMass3l, MAX_dilepMass3l);
    value_dilepMass4l = IN_RANGE(leptonVars.dilepMass, MIN_dilepMass4l, MAX_dilepMass4l);
    value_missingEt = IN_RANGE(jetVars.missingEt, MIN_missingEt, MAX_missingEt);
    value_missingPhi = PHI_IN_RANGE(jetVars.missingPhi);
    value_firstLepPt = IN_RANGE(leptonVars.firstLepPt, MIN_firstLepPt, MAX_firstLepPt);
    value_firstLepEta = IN_RANGE(leptonVars.firstLepEta, MIN_firstLepEta, MAX_firstLepEta);
    value_firstLepPhi = PHI_IN_RANGE(leptonVars.firstLepPhi);
    value_secondLepPt = IN_RANGE(leptonVars.secondLepPt, MIN_secondLepPt, MAX_secondLepPt);
    value_secondLepEta = IN_RANGE(leptonVars.secondLepEta, MIN_secondLepEta, MAX_secondLepEta);
    value_secondLepPhi = PHI_IN_RANGE(leptonVars.secondLepPhi);
    value_thirdLepPt = IN_RANGE(leptonVars.thirdLepPt, MIN_thirdLepPt, MAX_thirdLepPt);
    value_thirdLepEta = IN_RANGE(leptonVars.thirdLepEta, MIN_thirdLepEta, MAX_thirdLepEta);
    value_thirdLepPhi = PHI_IN_RANGE(leptonVars.thirdLepPhi);
    value_fourthLepPt = IN_RANGE(leptonVars.fourthLepPt, MIN_fourthLepPt, MAX_fourthLepPt);
    value_fourthLepEta = IN_RANGE(leptonVars.fourthLepEta, MIN_fourthLepEta, MAX_fourthLepEta);
    value_fourthLepPhi = PHI_IN_RANGE(leptonVars.thirdLepPhi);
    value_firstJetPt = IN_RANGE(jetVars.firstJetPt, MIN_firstJetPt, MAX_firstJetPt);
    value_firstJetEta = IN_RANGE(jetVars.firstJetEta, MIN_firstJetEta, MAX_firstJetEta);
    value_firstJetPhi = PHI_IN_RANGE(jetVars.firstJetPhi);
    value_secondJetPt = IN_RANGE(jetVars.secondJetPt, MIN_secondJetPt, MAX_secondJetPt);
    value_secondJetEta = IN_RANGE(jetVars.secondJetEta, MIN_secondJetEta, MAX_secondJetEta);
    value_secondJetPhi = PHI_IN_RANGE(jetVars.secondJetPhi);
    value_thirdJetPt = IN_RANGE(jetVars.thirdJetPt, MIN_thirdJetPt, MAX_thirdJetPt);
    value_thirdJetEta = IN_RANGE(jetVars.thirdJetEta, MIN_thirdJetEta, MAX_thirdJetEta);
    value_thirdJetPhi = PHI_IN_RANGE(jetVars.thirdJetPhi);
    value_fourthJetPt = IN_RANGE(jetVars.fourthJetPt, MIN_fourthJetPt, MAX_fourthJetPt);
    value_fourthJetEta = IN_RANGE(jetVars.fourthJetEta, MIN_fourthJetEta, MAX_fourthJetEta);
    value_fourthJetPhi = PHI_IN_RANGE(jetVars.fourthJetPhi);
    value_zbosonPt = IN_RANGE(leptonVars.dilepPt, BINS_zbosonPt[0], BINS_zbosonPt[NBINS_zbosonPt]);
    value_ttzMass = IN_RANGE(reconstructedVars.ttzMass, BINS_ttzMass[0], BINS_ttzMass[NBINS_ttzMass]);
    value_ttbarMass = IN_RANGE(reconstructedVars.ttbarMass, BINS_ttbarMass[0], BINS_ttbarMass[NBINS_ttbarMass]);
    value_topPt = IN_RANGE(reconstructedVars.topPt, BINS_topPt[0], BINS_topPt[NBINS_topPt]);
    value_deltaPhiTtbar = IN_RANGE(reconstructedVars.deltaPhiTtbar, BINS_deltaPhiTtbar[0], BINS_deltaPhiTtbar[NBINS_deltaPhiTtbar]);
    value_deltaPhiTopZ = IN_RANGE(reconstructedVars.deltaPhiTopZ, BINS_deltaPhiTopZ[0], BINS_deltaPhiTopZ[NBINS_deltaPhiTopZ]);
    value_deltaRapTtbar = IN_RANGE(reconstructedVars.deltaRapTtbar, BINS_deltaRapTtbar[0], BINS_deltaRapTtbar[NBINS_deltaRapTtbar]);
    value_deltaRapTopZ = IN_RANGE(reconstructedVars.deltaRapTopZ, BINS_deltaRapTopZ[0], BINS_deltaRapTopZ[NBINS_deltaRapTopZ]);
    value_lepTopMass = IN_RANGE(reconstructedVars.lepTopMass, MIN_lepTopMass, MAX_lepTopMass);
    value_hadTopMass = IN_RANGE(reconstructedVars.hadTopMass, MIN_hadTopMass, MAX_hadTopMass);
    value_topLeptonPt = IN_RANGE(fourLeptonVars.topLeptonPt, MIN_topLeptonPt, MAX_topLeptonPt);
    value_topLeptonsMass = IN_RANGE(fourLeptonVars.topLeptonsMass, MIN_topLeptonsMass, MAX_topLeptonsMass);
    value_fourLeptonsMass = IN_RANGE(fourLeptonVars.fourLeptonsMass, MIN_fourLeptonsMass, MAX_fourLeptonsMass);
    value_deltaPhiTopLeptons = IN_RANGE(fourLeptonVars.deltaPhiTopLeptons, MIN_deltaPhiTopLeptons, MAX_deltaPhiTopLeptons);
    value_deltaPhiTopLeptonZ = IN_RANGE(fourLeptonVars.deltaPhiTopLeptonZ, MIN_deltaPhiTopLeptonZ, MAX_deltaPhiTopLeptonZ);
    value_deltaRapTopLeptons = IN_RANGE(fourLeptonVars.deltaRapTopLeptons, BINS_deltaRapTopLeptons[0], BINS_deltaRapTopLeptons[NBINS_deltaRapTopLeptons]);
    value_deltaRapTopLeptonZ = IN_RANGE(fourLeptonVars.deltaRapTopLeptonZ, BINS_deltaRapTopLeptonZ[0], BINS_deltaRapTopLeptonZ[NBINS_deltaRapTopLeptonZ]);
}

void ttZHistograms::Fill(double weight, int variant, int selection, int systematic) {
    Fill(weight, variant+nVariants*selection+nVariants*nSelections*(1+systematic));
    if(systematic==0) Fill(weight>0.0 ? 1.0 : -1.0, variant+nVariants*selection);
}

void ttZHistograms::Fill(double weight, int histogram) {
    hists_category3l[histogram]->Fill(value_category3l, weight);
    hists_category4l[histogram]->Fill(value_category4l, weight);
    hists_nJets[histogram]->Fill(value_nJets, weight);
    hists_nBjets[histogram]->Fill(value_nBjets, weight);
    hists_dilepPt[histogram]->Fill(value_dilepPt, weight);
    hists_dilepEta[histogram]->Fill(value_dilepEta, weight);
    hists_dilepPhi[histogram]->Fill(value_dilepPhi, weight);
    hists_dilepMass3l[histogram]->Fill(value_dilepMass3l, weight);
    hists_dilepMass4l[histogram]->Fill(value_dilepMass4l, weight);
    hists_missingEt[histogram]->Fill(value_missingEt, weight);
    hists_missingPhi[histogram]->Fill(value_missingPhi, weight);
    hists_firstLepPt[histogram]->Fill(value_firstLepPt, weight);
    hists_firstLepEta[histogram]->Fill(value_firstLepEta, weight);
    hists_firstLepPhi[histogram]->Fill(value_firstLepPhi, weight);
    hists_secondLepPt[histogram]->Fill(value_secondLepPt, weight);
    hists_secondLepEta[histogram]->Fill(value_secondLepEta, weight);
    hists_secondLepPhi[histogram]->Fill(value_secondLepPhi, weight);
    hists_thirdLepPt[histogram]->Fill(value_thirdLepPt, weight);
    hists_thirdLepEta[histogram]->Fill(value_thirdLepEta, weight);
    hists_thirdLepPhi[histogram]->Fill(value_thirdLepPhi, weight);
    hists_fourthLepPt[histogram]->Fill(value_fourthLepPt, weight);
    hists_fourthLepEta[histogram]->Fill(value_fourthLepEta, weight);
    hists_fourthLepPhi[histogram]->Fill(value_fourthLepPhi, weight);
    hists_firstJetPt[histogram]->Fill(value_firstJetPt, weight);
    hists_firstJetEta[histogram]->Fill(value_firstJetEta, weight);
    hists_firstJetPhi[histogram]->Fill(value_firstJetPhi, weight);
    hists_secondJetPt[histogram]->Fill(value_secondJetPt, weight);
    hists_secondJetEta[histogram]->Fill(value_secondJetEta, weight);
    hists_secondJetPhi[histogram]->Fill(value_secondJetPhi, weight);
    hists_thirdJetPt[histogram]->Fill(value_thirdJetPt, weight);
    hists_thirdJetEta[histogram]->Fill(value_thirdJetEta, weight);
    hists_thirdJetPhi[histogram]->Fill(value_thirdJetPhi, weight);
    hists_fourthJetPt[histogram]->Fill(value_fourthJetPt, weight);
    hists_fourthJetEta[histogram]->Fill(value_fourthJetEta, weight);
    hists_fourthJetPhi[histogram]->Fill(value_fourthJetPhi, weight);
    hists_zbosonPt[histogram]->Fill(value_zbosonPt, weight);
    hists_ttzMass[histogram]->Fill(value_ttzMass, weight);
    hists_ttbarMass[histogram]->Fill(value_ttbarMass, weight);
    hists_topPt[histogram]->Fill(value_topPt, weight);
    hists_deltaPhiTtbar[histogram]->Fill(value_deltaPhiTtbar, weight);
    hists_deltaPhiTopZ[histogram]->Fill(value_deltaPhiTopZ, weight);
    hists_deltaRapTtbar[histogram]->Fill(value_deltaRapTtbar, weight);
    hists_deltaRapTopZ[histogram]->Fill(value_deltaRapTopZ, weight);
    hists_lepTopMass[histogram]->Fill(value_lepTopMass, weight);
    hists_hadTopMass[histogram]->Fill(value_hadTopMass, weight);
    hists_topLeptonPt[histogram]->Fill(value_topLeptonPt, weight);
    hists_topLeptonsMass[histogram]->Fill(value_topLeptonsMass, weight);
    hists_fourLeptonsMass[histogram]->Fill(value_fourLeptonsMass, weight);
    hists_deltaPhiTopLeptons[histogram]->Fill(value_deltaPhiTopLeptons, weight);
    hists_deltaPhiTopLeptonZ[histogram]->Fill(value_deltaPhiTopLeptonZ, weight);
    hists_deltaRapTopLeptons[histogram]->Fill(value_deltaRapTopLeptons, weight);
    hists_deltaRapTopLeptonZ[histogram]->Fill(value_deltaRapTopLeptonZ, weight);
}

void ttZHistograms::Write(TDirectory* directory) {
    Write(directory, "control_category3l", hists_category3l);
    Write(directory, "control_category4l", hists_category4l);
    Write(directory, "control_nJets", hists_nJets);
    Write(directory, "control_nBjets", hists_nBjets);
    Write(directory, "control_dilepPt", hists_dilepPt);
    Write(directory, "control_dilepEta", hists_dilepEta);
    Write(directory, "control_dilepPhi", hists_dilepPhi);
    Write(directory, "control_dilepMass3l", hists_dilepMass3l);
    Write(directory, "control_dilepMass4l", hists_dilepMass4l);
    Write(directory, "control_missingEt", hists_missingEt);
    Write(directory, "control_missingPhi", hists_missingPhi);
    Write(directory, "control_firstLepPt", hists_firstLepPt);
    Write(directory, "control_firstLepEta", hists_firstLepEta);
    Write(directory, "control_firstLepPhi", hists_firstLepPhi);
    Write(directory, "control_secondLepPt", hists_secondLepPt);
    Write(directory, "control_secondLepEta", hists_secondLepEta);
    Write(directory, "control_secondLepPhi", hists_secondLepPhi);
    Write(directory, "control_thirdLepPt", hists_thirdLepPt);
    Write(directory, "control_thirdLepEta", hists_thirdLepEta);
    Write(directory, "control_thirdLepPhi", hists_thirdLepPhi);
    Write(directory, "control_fourthLepPt", hists_fourthLepPt);
    Write(directory, "control_fourthLepEta", hists_fourthLepEta);
    Write(directory, "control_fourthLepPhi", hists_fourthLepPhi);
    Write(directory, "control_firstJetPt", hists_firstJetPt);
    Write(directory, "control_firstJetEta", hists_firstJetEta);
    Write(directory, "control_firstJetPhi", hists_firstJetPhi);
    Write(directory, "control_secondJetPt", hists_secondJetPt);
    Write(directory, "control_secondJetEta", hists_secondJetEta);
    Write(directory, "control_secondJetPhi", hists_secondJetPhi);
    Write(directory, "control_thirdJetPt", hists_thirdJetPt);
    Write(directory, "control_thirdJetEta", hists_thirdJetEta);
    Write(directory, "control_thirdJetPhi", hists_thirdJetPhi);
    Write(directory, "control_fourthJetPt", hists_fourthJetPt);
    Write(directory, "control_fourthJetEta", hists_fourthJetEta);
    Write(directory, "control_fourthJetPhi", hists_fourthJetPhi);
    Write(directory, "reco_zbosonPt", hists_zbosonPt);
    Write(directory, "reco3l_ttzMass", hists_ttzMass);
    Write(directory, "reco3l_ttbarMass", hists_ttbarMass);
    Write(directory, "reco3l_topPt", hists_topPt);
    Write(directory, "reco3l_deltaPhiTtbar", hists_deltaPhiTtbar);
    Write(directory, "reco3l_deltaPhiTopZ", hists_deltaPhiTopZ);
    Write(directory, "reco3l_deltaRapTtbar", hists_deltaRapTtbar);
    Write(directory, "reco3l_deltaRapTopZ", hists_deltaRapTopZ);
    Write(directory, "reco3l_lepTopMass", hists_lepTopMass);
    Write(directory, "reco3l_hadTopMass", hists_hadTopMass);
    Write(directory, "reco4l_topLeptonPt", hists_topLeptonPt);
    Write(directory, "reco4l_topLeptonsMass", hists_topLeptonsMass);
    Write(directory, "reco4l_fourLeptonsMass", hists_fourLeptonsMass);
    Write(directory, "reco4l_deltaPhiTopLeptons", hists_deltaPhiTopLeptons);
    Write(directory, "reco4l_deltaPhiTopLeptonZ", hists_deltaPhiTopLeptonZ);
    Write(directory, "reco4l_deltaRapTopLeptons", hists_deltaRapTopLeptons);
    Write(directory, "reco4l_deltaRapTopLeptonZ", hists_deltaRapTopLeptonZ);
}

void ttZHistograms::Write(TDirectory* directory, std::string name, TH1F** hists) {
    directory->mkdir(name.c_str())->cd();
    for(int iVar=0; iVar<nVariants; iVar++) {
        const std::string sVar = "_var"+std::to_string(iVar);
        for(int iSel=0; iSel<nSelections; iSel++) {
            const std::string sSel = "_sel"+std::to_string(iSel);
            hists[iVar+nVariants*iSel]->Write((name+sVar+sSel+"_stat").c_str());
            for(int iSys=0; iSys<nSystematics; iSys++) {
                const std::string sSys = "_sys"+std::to_string(iSys);
                hists[iVar+nVariants*iSel+nVariants*nSelections*(1+iSys)]->Write((name+sVar+sSel+sSys).c_str());
            }
        }
    }
}
