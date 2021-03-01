#ifndef ttZHistograms_H
#define ttZHistograms_H

//include c++ library classes
#include <string>

//include ROOT classes
#include "TDirectory.h"
#include "TMath.h"
#include "TH1.h"

//include ttZ specific code
#include "ttZVariables.h"

class ttZHistograms {
public:
    ttZHistograms(int, int, int);
    virtual ~ttZHistograms();

    void SetValues(ttZ::leptonVariables, ttZ::jetVariables, ttZ::reconstructedVariables, ttZ::fourLeptonVariables);
    void Fill(double, int, int, int);

    void Write(TDirectory*);

protected:
    void Fill(double, int);
    void Write(TDirectory*, std::string, TH1F**);

    const int nVariants, nSelections, nSystematics;
    const int nHistograms;
    TH1F** hists_category3l=nullptr, ** hists_category4l=nullptr,
        ** hists_nJets=nullptr, ** hists_nBjets=nullptr,
        ** hists_dilepPt=nullptr, ** hists_dilepEta=nullptr, ** hists_dilepPhi=nullptr, ** hists_dilepMass3l=nullptr, ** hists_dilepMass4l=nullptr,
        ** hists_missingEt=nullptr, ** hists_missingPhi=nullptr,
        ** hists_firstLepPt=nullptr, ** hists_firstLepEta=nullptr, ** hists_firstLepPhi=nullptr,
        ** hists_secondLepPt=nullptr, ** hists_secondLepEta=nullptr, ** hists_secondLepPhi=nullptr,
        ** hists_thirdLepPt=nullptr, ** hists_thirdLepEta=nullptr, ** hists_thirdLepPhi=nullptr,
        ** hists_fourthLepPt=nullptr, ** hists_fourthLepEta=nullptr, ** hists_fourthLepPhi=nullptr,
        ** hists_firstJetPt=nullptr, ** hists_firstJetEta=nullptr, ** hists_firstJetPhi=nullptr,
        ** hists_secondJetPt=nullptr, ** hists_secondJetEta=nullptr, ** hists_secondJetPhi=nullptr,
        ** hists_thirdJetPt=nullptr, ** hists_thirdJetEta=nullptr, ** hists_thirdJetPhi=nullptr,
        ** hists_fourthJetPt=nullptr, ** hists_fourthJetEta=nullptr, ** hists_fourthJetPhi=nullptr,
        ** hists_ttzMass=nullptr, ** hists_ttbarMass=nullptr, ** hists_topPt=nullptr,
        ** hists_deltaPhiTtbar=nullptr, ** hists_deltaPhiTopZ=nullptr, ** hists_deltaRapTtbar=nullptr, ** hists_deltaRapTopZ=nullptr,
        ** hists_lepTopMass=nullptr, ** hists_hadTopMass=nullptr,
        ** hists_topLeptonPt=nullptr, ** hists_topLeptonsMass=nullptr, ** hists_fourLeptonsMass=nullptr,
        ** hists_deltaPhiTopLeptons=nullptr, ** hists_deltaPhiTopLeptonZ=nullptr, ** hists_deltaRapTopLeptons=nullptr, ** hists_deltaRapTopLeptonZ=nullptr;
    int value_category3l, value_category4l,
        value_nJets, value_nBjets;
    double value_dilepPt, value_dilepEta, value_dilepPhi, value_dilepMass3l, value_dilepMass4l,
           value_missingEt, value_missingPhi,
           value_firstLepPt, value_firstLepEta, value_firstLepPhi,
           value_secondLepPt, value_secondLepEta, value_secondLepPhi,
           value_thirdLepPt, value_thirdLepEta, value_thirdLepPhi,
           value_fourthLepPt, value_fourthLepEta, value_fourthLepPhi,
           value_firstJetPt, value_firstJetEta, value_firstJetPhi,
           value_secondJetPt, value_secondJetEta, value_secondJetPhi,
           value_thirdJetPt, value_thirdJetEta, value_thirdJetPhi,
           value_fourthJetPt, value_fourthJetEta, value_fourthJetPhi,
           value_ttzMass, value_ttbarMass, value_topPt,
           value_deltaPhiTtbar, value_deltaPhiTopZ, value_deltaRapTtbar, value_deltaRapTopZ,
           value_lepTopMass, value_hadTopMass,
           value_topLeptonPt, value_topLeptonsMass, value_fourLeptonsMass,
           value_deltaPhiTopLeptons, value_deltaPhiTopLeptonZ, value_deltaRapTopLeptons, value_deltaRapTopLeptonZ;
};

#endif
