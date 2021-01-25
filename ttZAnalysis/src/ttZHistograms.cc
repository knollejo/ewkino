#include "../interface/ttZHistograms.h"

//include ROOT classes
#include "TMath.h"
#include "TVector2.h"

//include other parts of framework

// histogram definitions
const int NBINS_category = 4; const double MIN_category = -0.5; const double MAX_category = 3.5;
const int NBINS_nJets = 5; const double MIN_nJets = 2.5; const double MAX_nJets = 7.5;
const int NBINS_nBjets = 3; const double MIN_nBjets = 0.5; const double MAX_nBjets = 7.5;
const int NBINS_dilepPt = 10; const double MIN_dilepPt = 0.0; const double MAX_dilepPt = 500.0;
const int NBINS_dilepEta = 12; const double MIN_dilepEta = -3.0; const double MAX_dilepEta = 3.0;
const int NBINS_dilepPhi = 10; const double MIN_dilepPhi = -TMath::Pi(); const double MAX_dilepPhi = TMath::Pi();
const int NBINS_dilepMass = 10; const double MIN_dilepMass = 81.0; const double MAX_dilepMass = 101.0;
const int NBINS_missingEt = 11; const double MIN_missingEt = 0.0; const double MAX_missingEt = 220.0;
const int NBINS_missingPhi = 10; const double MIN_missingPhi = -TMath::Pi(); const double MAX_missingPhi = TMath::Pi();
const int NBINS_leadLepPt = 10; const double MIN_leadLepPt = 40.0; const double MAX_leadLepPt = 330.0;
const int NBINS_leadLepEta = 10; const double MIN_leadLepEta = -2.5; const double MAX_leadLepEta = 2.5;
const int NBINS_leadLepPhi = 10; const double MIN_leadLepPhi = -TMath::Pi(); const double MAX_leadLepPhi = -TMath::Pi();
const int NBINS_sublLepPt = 10; const double MIN_sublLepPt = 20.0; const double MAX_sublLepPt = 165.0;
const int NBINS_sublLepEta = 10; const double MIN_sublLepEta = -2.5; const double MAX_sublLepEta = 2.5;
const int NBINS_sublLepPhi = 10; const double MIN_sublLepPhi = -TMath::Pi(); const double MAX_sublLepPhi = -TMath::Pi();
const int NBINS_trailLepPt = 10; const double MIN_trailLepPt = 10.0; const double MAX_trailLepPt = 110.0;
const int NBINS_trailLepEta = 10; const double MIN_trailLepEta = -2.5; const double MAX_trailLepEta = 2.5;
const int NBINS_trailLepPhi = 10; const double MIN_trailLepPhi = -TMath::Pi(); const double MAX_trailLepPhi = -TMath::Pi();
const int NBINS_firstJetPt = 10; const double MIN_firstJetPt = 30.0; const double MAX_firstJetPt = 330.0;
const int NBINS_firstJetEta = 10; const double MIN_firstJetEta = -2.5; const double MAX_firstJetEta = 2.5;
const int NBINS_firstJetPhi = 10; const double MIN_firstJetPhi = -TMath::Pi(); const double MAX_firstJetPhi = -TMath::Pi();
const int NBINS_secondJetPt = 10; const double MIN_secondJetPt = 30.0; const double MAX_secondJetPt = 330.0;
const int NBINS_secondJetEta = 10; const double MIN_secondJetEta = -2.5; const double MAX_secondJetEta = 2.5;
const int NBINS_secondJetPhi = 10; const double MIN_secondJetPhi = -TMath::Pi(); const double MAX_secondJetPhi = -TMath::Pi();
const int NBINS_thirdJetPt = 10; const double MIN_thirdJetPt = 30.0; const double MAX_thirdJetPt = 330.0;
const int NBINS_thirdJetEta = 10; const double MIN_thirdJetEta = -2.5; const double MAX_thirdJetEta = 2.5;
const int NBINS_thirdJetPhi = 10; const double MIN_thirdJetPhi = -TMath::Pi(); const double MAX_thirdJetPhi = -TMath::Pi();
const int NBINS_fourthJetPt = 10; const double MIN_fourthJetPt = 30.0; const double MAX_fourthJetPt = 330.0;
const int NBINS_fourthJetEta = 10; const double MIN_fourthJetEta = -2.5; const double MAX_fourthJetEta = 2.5;
const int NBINS_fourthJetPhi = 10; const double MIN_fourthJetPhi = -TMath::Pi(); const double MAX_fourthJetPhi = -TMath::Pi();
const int NBINS_ttzMass = 8; const double BINS_ttzMass[NBINS_ttzMass+1] = { 436.0, 580.0, 640.0, 720.0, 800.0, 890.0, 1000.0, 1150.0, 1500.0 };
const int NBINS_ttbarMass = 8; const double BINS_ttbarMass[NBINS_ttbarMass+1] = { 345.0, 400.0, 440.0, 485.0, 540.0, 610.0, 700.0, 850.0, 1200.0 };
const int NBINS_topPt = 8; const double BINS_topPt[NBINS_topPt+1] = { 0.0, 60.0, 90.0, 125.0, 160.0, 200.0, 250.0, 330.0, 500.0 };
const int NBINS_deltaPhiTtbar = 8; const double BINS_deltaPhiTtbar[NBINS_deltaPhiTtbar+1] = { 0.0, TMath::Pi()/4, TMath::Pi()/2, 5*TMath::Pi()/8, 3*TMath::Pi()/4, 13*TMath::Pi()/16, 7*TMath::Pi()/8, 15*TMath::Pi()/16, TMath::Pi() };
const int NBINS_deltaPhiTopZ = 8; const double BINS_deltaPhiTopZ[NBINS_deltaPhiTopZ+1] = { 0.0, TMath::Pi()/6, TMath::Pi()/3, TMath::Pi()/2, 2*TMath::Pi()/3, 3*TMath::Pi()/4, 5*TMath::Pi()/6, 11*TMath::Pi()/12, TMath::Pi() };
const int NBINS_deltaRapTtbar = 8; const double BINS_deltaRapTtbar[NBINS_deltaRapTtbar+1] = { 0.0, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 2.1, 3.0 };
const int NBINS_deltaRapTopZ = 8; const double BINS_deltaRapTopZ[NBINS_deltaRapTopZ+1] = { 0.0, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 2.0, 3.0 };

ttZHistograms::ttZHistograms(int nVar, int nSel, int nSys) :
    nVariants(nVar), nSelections(nSel), nSystematics(nSys),
    nHistograms(nVariants*nSelections*(1+nSystematics)),
    hists_category(new TH1F*[nHistograms]),
    hists_nJets(new TH1F*[nHistograms]),
    hists_nBjets(new TH1F*[nHistograms]),
    hists_dilepPt(new TH1F*[nHistograms]),
    hists_dilepEta(new TH1F*[nHistograms]),
    hists_dilepPhi(new TH1F*[nHistograms]),
    hists_dilepMass(new TH1F*[nHistograms]),
    hists_missingEt(new TH1F*[nHistograms]),
    hists_missingPhi(new TH1F*[nHistograms]),
    hists_leadLepPt(new TH1F*[nHistograms]),
    hists_leadLepEta(new TH1F*[nHistograms]),
    hists_leadLepPhi(new TH1F*[nHistograms]),
    hists_sublLepPt(new TH1F*[nHistograms]),
    hists_sublLepEta(new TH1F*[nHistograms]),
    hists_sublLepPhi(new TH1F*[nHistograms]),
    hists_trailLepPt(new TH1F*[nHistograms]),
    hists_trailLepEta(new TH1F*[nHistograms]),
    hists_trailLepPhi(new TH1F*[nHistograms]),
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
    hists_ttzMass(new TH1F*[nHistograms]),
    hists_ttbarMass(new TH1F*[nHistograms]),
    hists_topPt(new TH1F*[nHistograms]),
    hists_deltaPhiTtbar(new TH1F*[nHistograms]),
    hists_deltaPhiTopZ(new TH1F*[nHistograms]),
    hists_deltaRapTtbar(new TH1F*[nHistograms]),
    hists_deltaRapTopZ(new TH1F*[nHistograms])
{
    for(int i=0; i<nHistograms; i++) {
        const std::string sI = std::to_string(i);
        hists_category[i] = new TH1F(("hist_category_var"+sI).c_str(), "", NBINS_category, MIN_category, MAX_category);
        hists_nJets[i] = new TH1F(("hist_nJets_var"+sI).c_str(), "", NBINS_nJets, MIN_nJets, MAX_nJets);
        hists_nBjets[i] = new TH1F(("hist_nBjets_var"+sI).c_str(), "", NBINS_nBjets, MIN_nBjets, MAX_nBjets);
        hists_dilepPt[i] = new TH1F(("hist_dilepPt_var"+sI).c_str(), "", NBINS_dilepPt, MIN_dilepPt, MAX_dilepPt);
        hists_dilepEta[i] = new TH1F(("hist_dilepEta_var"+sI).c_str(), "", NBINS_dilepEta, MIN_dilepEta, MAX_dilepEta);
        hists_dilepPhi[i] = new TH1F(("hist_dilepPhi_var"+sI).c_str(), "", NBINS_dilepPhi, MIN_dilepPhi, MAX_dilepPhi);
        hists_dilepMass[i] = new TH1F(("hist_dilepMass_var"+sI).c_str(), "", NBINS_dilepMass, MIN_dilepMass, MAX_dilepMass);
        hists_missingEt[i] = new TH1F(("hist_missingEt_var"+sI).c_str(), "", NBINS_missingEt, MIN_missingEt, MAX_missingEt);
        hists_missingPhi[i] = new TH1F(("hist_missingPhi_var"+sI).c_str(), "", NBINS_missingPhi, MIN_missingPhi, MAX_missingPhi);
        hists_leadLepPt[i] = new TH1F(("hist_leadLepPt_var"+sI).c_str(), "", NBINS_leadLepPt, MIN_leadLepPt, MAX_leadLepPt);
        hists_leadLepEta[i] = new TH1F(("hist_leadLepEta_var"+sI).c_str(), "", NBINS_leadLepEta, MIN_leadLepEta, MAX_leadLepEta);
        hists_leadLepPhi[i] = new TH1F(("hist_leadLepPhi_var"+sI).c_str(), "", NBINS_leadLepPhi, MIN_leadLepPhi, MAX_leadLepPhi);
        hists_sublLepPt[i] = new TH1F(("hist_sublLepPt_var"+sI).c_str(), "", NBINS_sublLepPt, MIN_sublLepPt, MAX_sublLepPt);
        hists_sublLepEta[i] = new TH1F(("hist_sublLepEta_var"+sI).c_str(), "", NBINS_sublLepEta, MIN_sublLepEta, MAX_sublLepEta);
        hists_sublLepPhi[i] = new TH1F(("hist_sublLepPhi_var"+sI).c_str(), "", NBINS_sublLepPhi, MIN_sublLepPhi, MAX_sublLepPhi);
        hists_trailLepPt[i] = new TH1F(("hist_trailLepPt_var"+sI).c_str(), "", NBINS_trailLepPt, MIN_trailLepPt, MAX_trailLepPt);
        hists_trailLepEta[i] = new TH1F(("hist_trailLepEta_var"+sI).c_str(), "", NBINS_trailLepEta, MIN_trailLepEta, MAX_trailLepEta);
        hists_trailLepPhi[i] = new TH1F(("hist_trailLepPhi_var"+sI).c_str(), "", NBINS_trailLepPhi, MIN_trailLepPhi, MAX_trailLepPhi);
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
        hists_ttzMass[i] = new TH1F(("hist_ttzMass_var"+sI).c_str(), "", NBINS_ttzMass, BINS_ttzMass);
        hists_ttbarMass[i] = new TH1F(("hist_ttbarMass_var"+sI).c_str(), "", NBINS_ttbarMass, BINS_ttbarMass);
        hists_topPt[i] = new TH1F(("hist_topPt_var"+sI).c_str(), "", NBINS_topPt, BINS_topPt);
        hists_deltaPhiTtbar[i] = new TH1F(("hist_deltaPhiTtbar_var"+sI).c_str(), "", NBINS_deltaPhiTtbar, BINS_deltaPhiTtbar);
        hists_deltaPhiTopZ[i] = new TH1F(("hist_deltaPhiTopZ_var"+sI).c_str(), "", NBINS_deltaPhiTopZ, BINS_deltaPhiTopZ);
        hists_deltaRapTtbar[i] = new TH1F(("hist_deltaRapTtbar_var"+sI).c_str(), "", NBINS_deltaRapTtbar, BINS_deltaRapTtbar);
        hists_deltaRapTopZ[i] = new TH1F(("hist_deltaRapTopZ_var"+sI).c_str(), "", NBINS_deltaRapTopZ, BINS_deltaRapTopZ);
    }
}

ttZHistograms::~ttZHistograms() {
    for(int i=0; i<nHistograms; i++) {
        delete hists_category[i];
        delete hists_nJets[i];
        delete hists_nBjets[i];
        delete hists_dilepPt[i];
        delete hists_dilepEta[i];
        delete hists_dilepPhi[i];
        delete hists_dilepMass[i];
        delete hists_missingEt[i];
        delete hists_missingPhi[i];
        delete hists_leadLepPt[i];
        delete hists_leadLepEta[i];
        delete hists_leadLepPhi[i];
        delete hists_sublLepPt[i];
        delete hists_sublLepEta[i];
        delete hists_sublLepPhi[i];
        delete hists_trailLepPt[i];
        delete hists_trailLepEta[i];
        delete hists_trailLepPhi[i];
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
        delete hists_ttzMass[i];
        delete hists_ttbarMass[i];
        delete hists_topPt[i];
        delete hists_deltaPhiTtbar[i];
        delete hists_deltaPhiTopZ[i];
        delete hists_deltaRapTtbar[i];
        delete hists_deltaRapTopZ[i];
    }
    delete [] hists_category;
    delete [] hists_nJets;
    delete [] hists_nBjets;
    delete [] hists_dilepPt;
    delete [] hists_dilepEta;
    delete [] hists_dilepPhi;
    delete [] hists_dilepMass;
    delete [] hists_missingEt;
    delete [] hists_missingPhi;
    delete [] hists_leadLepPt;
    delete [] hists_leadLepEta;
    delete [] hists_leadLepPhi;
    delete [] hists_sublLepPt;
    delete [] hists_sublLepEta;
    delete [] hists_sublLepPhi;
    delete [] hists_trailLepPt;
    delete [] hists_trailLepEta;
    delete [] hists_trailLepPhi;
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
    delete [] hists_ttzMass;
    delete [] hists_ttbarMass;
    delete [] hists_topPt;
    delete [] hists_deltaPhiTtbar;
    delete [] hists_deltaPhiTopZ;
    delete [] hists_deltaRapTtbar;
    delete [] hists_deltaRapTopZ;
}

#define IN_RANGE(VAL, MIN, MAX) \
    (VAL<-900.0 ? -999.0 : TMath::Min(MAX-0.001*(MAX-MIN), TMath::Max(MIN+0.001*(MAX-MIN), VAL)))
#define PHI_IN_RANGE(VAL) \
    (VAL<-900.0 || isnan(VAL) ? -999.0 : TMath::Min(0.999*TMath::Pi(), TMath::Max(-0.999*TMath::Pi(), TVector2::Phi_mpi_pi(VAL))))

void ttZHistograms::SetValues(ttZ::leptonVariables leptonVars, ttZ::jetVariables jetVars, ttZ::reconstructedVariables reconstructedVars) {
    value_category = leptonVars.category;
    value_nJets = jetVars.nJets>MAX_nJets ? MAX_nJets-0.5 : jetVars.nJets;
    value_nBjets = jetVars.nBjets>MAX_nBjets ? MAX_nBjets-0.5 : jetVars.nBjets;
    value_dilepPt = IN_RANGE(leptonVars.dilepPt, MIN_dilepPt, MAX_dilepPt);
    value_dilepEta = IN_RANGE(leptonVars.dilepEta, MIN_dilepEta, MAX_dilepEta);
    value_dilepPhi = PHI_IN_RANGE(leptonVars.dilepPhi);
    value_dilepMass = IN_RANGE(leptonVars.dilepMass, MIN_dilepMass, MAX_dilepMass);
    value_missingEt = IN_RANGE(jetVars.missingEt, MIN_missingEt, MAX_missingEt);
    value_missingPhi = PHI_IN_RANGE(jetVars.missingPhi);
    value_leadLepPt = IN_RANGE(leptonVars.leadLepPt, MIN_leadLepPt, MAX_leadLepPt);
    value_leadLepEta = IN_RANGE(leptonVars.leadLepEta, MIN_leadLepEta, MAX_leadLepEta);
    value_leadLepPhi = PHI_IN_RANGE(leptonVars.leadLepPhi);
    value_sublLepPt = IN_RANGE(leptonVars.sublLepPt, MIN_sublLepPt, MAX_sublLepPt);
    value_sublLepEta = IN_RANGE(leptonVars.sublLepEta, MIN_sublLepEta, MAX_sublLepEta);
    value_sublLepPhi = PHI_IN_RANGE(leptonVars.sublLepPhi);
    value_trailLepPt = IN_RANGE(leptonVars.trailLepPt, MIN_trailLepPt, MAX_trailLepPt);
    value_trailLepEta = IN_RANGE(leptonVars.trailLepEta, MIN_trailLepEta, MAX_trailLepEta);
    value_trailLepPhi = PHI_IN_RANGE(leptonVars.trailLepPhi);
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
    value_ttzMass = IN_RANGE(reconstructedVars.ttzMass, BINS_ttzMass[0], BINS_ttzMass[NBINS_ttzMass]);
    value_ttbarMass = IN_RANGE(reconstructedVars.ttbarMass, BINS_ttbarMass[0], BINS_ttbarMass[NBINS_ttbarMass]);
    value_topPt = IN_RANGE(reconstructedVars.topPt, BINS_topPt[0], BINS_topPt[NBINS_topPt]);
    value_deltaPhiTtbar = IN_RANGE(reconstructedVars.deltaPhiTtbar, BINS_deltaPhiTtbar[0], BINS_deltaPhiTtbar[NBINS_deltaPhiTtbar]);
    value_deltaPhiTopZ = IN_RANGE(reconstructedVars.deltaPhiTopZ, BINS_deltaPhiTopZ[0], BINS_deltaPhiTopZ[NBINS_deltaPhiTopZ]);
    value_deltaRapTtbar = IN_RANGE(reconstructedVars.deltaRapTtbar, BINS_deltaRapTtbar[0], BINS_deltaRapTtbar[NBINS_deltaRapTtbar]);
    value_deltaRapTopZ = IN_RANGE(reconstructedVars.deltaRapTopZ, BINS_deltaRapTopZ[0], BINS_deltaRapTopZ[NBINS_deltaRapTopZ]);
}

void ttZHistograms::Fill(double weight, int variant, int selection, int systematic) {
    Fill(weight, variant+nVariants*selection+nVariants*nSelections*(1+systematic));
    if(systematic==0) Fill(weight>0.0 ? 1.0 : -1.0, variant+nVariants*selection);
}

void ttZHistograms::Fill(double weight, int histogram) {
    hists_category[histogram]->Fill(value_category, weight);
    hists_nJets[histogram]->Fill(value_nJets, weight);
    hists_nBjets[histogram]->Fill(value_nBjets, weight);
    hists_dilepPt[histogram]->Fill(value_dilepPt, weight);
    hists_dilepEta[histogram]->Fill(value_dilepEta, weight);
    hists_dilepPhi[histogram]->Fill(value_dilepPhi, weight);
    hists_dilepMass[histogram]->Fill(value_dilepMass, weight);
    hists_missingEt[histogram]->Fill(value_missingEt, weight);
    hists_missingPhi[histogram]->Fill(value_missingPhi, weight);
    hists_leadLepPt[histogram]->Fill(value_leadLepPt, weight);
    hists_leadLepEta[histogram]->Fill(value_leadLepEta, weight);
    hists_leadLepPhi[histogram]->Fill(value_leadLepPhi, weight);
    hists_sublLepPt[histogram]->Fill(value_sublLepPt, weight);
    hists_sublLepEta[histogram]->Fill(value_sublLepEta, weight);
    hists_sublLepPhi[histogram]->Fill(value_sublLepPhi, weight);
    hists_trailLepPt[histogram]->Fill(value_trailLepPt, weight);
    hists_trailLepEta[histogram]->Fill(value_trailLepEta, weight);
    hists_trailLepPhi[histogram]->Fill(value_trailLepPhi, weight);
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
    hists_ttzMass[histogram]->Fill(value_ttzMass, weight);
    hists_ttbarMass[histogram]->Fill(value_ttbarMass, weight);
    hists_topPt[histogram]->Fill(value_topPt, weight);
    hists_deltaPhiTtbar[histogram]->Fill(value_deltaPhiTtbar, weight);
    hists_deltaPhiTopZ[histogram]->Fill(value_deltaPhiTopZ, weight);
    hists_deltaRapTtbar[histogram]->Fill(value_deltaRapTtbar, weight);
    hists_deltaRapTopZ[histogram]->Fill(value_deltaRapTopZ, weight);
}

void ttZHistograms::Write(TDirectory* directory) {
    Write(directory, "category", hists_category);
    Write(directory, "nJets", hists_nJets);
    Write(directory, "nBjets", hists_nBjets);
    Write(directory, "dilepPt", hists_dilepPt);
    Write(directory, "dilepEta", hists_dilepEta);
    Write(directory, "dilepPhi", hists_dilepPhi);
    Write(directory, "dilepMass", hists_dilepMass);
    Write(directory, "missingEt", hists_missingEt);
    Write(directory, "missingPhi", hists_missingPhi);
    Write(directory, "leadLepPt", hists_leadLepPt);
    Write(directory, "leadLepEta", hists_leadLepEta);
    Write(directory, "leadLepPhi", hists_leadLepPhi);
    Write(directory, "sublLepPt", hists_sublLepPt);
    Write(directory, "sublLepEta", hists_sublLepEta);
    Write(directory, "sublLepPhi", hists_sublLepPhi);
    Write(directory, "trailLepPt", hists_trailLepPt);
    Write(directory, "trailLepEta", hists_trailLepEta);
    Write(directory, "trailLepPhi", hists_trailLepPhi);
    Write(directory, "firstJetPt", hists_firstJetPt);
    Write(directory, "firstJetEta", hists_firstJetEta);
    Write(directory, "firstJetPhi", hists_firstJetPhi);
    Write(directory, "secondJetPt", hists_secondJetPt);
    Write(directory, "secondJetEta", hists_secondJetEta);
    Write(directory, "secondJetPhi", hists_secondJetPhi);
    Write(directory, "thirdJetPt", hists_thirdJetPt);
    Write(directory, "thirdJetEta", hists_thirdJetEta);
    Write(directory, "thirdJetPhi", hists_thirdJetPhi);
    Write(directory, "fourthJetPt", hists_fourthJetPt);
    Write(directory, "fourthJetEta", hists_fourthJetEta);
    Write(directory, "fourthJetPhi", hists_fourthJetPhi);
    Write(directory, "ttzMass", hists_ttzMass);
    Write(directory, "ttbarMass", hists_ttbarMass);
    Write(directory, "topPt", hists_topPt);
    Write(directory, "deltaPhiTtbar", hists_deltaPhiTtbar);
    Write(directory, "deltaPhiTopZ", hists_deltaPhiTopZ);
    Write(directory, "deltaRapTtbar", hists_deltaRapTtbar);
    Write(directory, "deltaRapTopZ", hists_deltaRapTopZ);
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
