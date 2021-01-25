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

    void SetValues(ttZ::leptonVariables, ttZ::jetVariables, ttZ::reconstructedVariables);
    void Fill(double, int, int, int);

    void Write(TDirectory*);

protected:
    void Fill(double, int);
    void Write(TDirectory*, std::string, TH1F**);

    const int nVariants, nSelections, nSystematics;
    const int nHistograms;
    TH1F** hists_category=nullptr,
        ** hists_nJets=nullptr, ** hists_nBjets=nullptr,
        ** hists_dilepPt=nullptr, ** hists_dilepEta=nullptr, ** hists_dilepPhi=nullptr, ** hists_dilepMass=nullptr,
        ** hists_missingEt=nullptr, ** hists_missingPhi=nullptr,
        ** hists_leadLepPt=nullptr, ** hists_leadLepEta=nullptr, ** hists_leadLepPhi=nullptr,
        ** hists_sublLepPt=nullptr, ** hists_sublLepEta=nullptr, ** hists_sublLepPhi=nullptr,
        ** hists_trailLepPt=nullptr, ** hists_trailLepEta=nullptr, ** hists_trailLepPhi=nullptr,
        ** hists_firstJetPt=nullptr, ** hists_firstJetEta=nullptr, ** hists_firstJetPhi=nullptr,
        ** hists_secondJetPt=nullptr, ** hists_secondJetEta=nullptr, ** hists_secondJetPhi=nullptr,
        ** hists_thirdJetPt=nullptr, ** hists_thirdJetEta=nullptr, ** hists_thirdJetPhi=nullptr,
        ** hists_fourthJetPt=nullptr, ** hists_fourthJetEta=nullptr, ** hists_fourthJetPhi=nullptr,
        ** hists_ttzMass=nullptr, ** hists_ttbarMass=nullptr, ** hists_topPt=nullptr,
        ** hists_deltaPhiTtbar=nullptr, ** hists_deltaPhiTopZ=nullptr, ** hists_deltaRapTtbar=nullptr, ** hists_deltaRapTopZ=nullptr;
    int value_category,
        value_nJets, value_nBjets;
    double value_dilepPt, value_dilepEta, value_dilepPhi, value_dilepMass,
           value_missingEt, value_missingPhi,
           value_leadLepPt, value_leadLepEta, value_leadLepPhi,
           value_sublLepPt, value_sublLepEta, value_sublLepPhi,
           value_trailLepPt, value_trailLepEta, value_trailLepPhi,
           value_firstJetPt, value_firstJetEta, value_firstJetPhi,
           value_secondJetPt, value_secondJetEta, value_secondJetPhi,
           value_thirdJetPt, value_thirdJetEta, value_thirdJetPhi,
           value_fourthJetPt, value_fourthJetEta, value_fourthJetPhi,
           value_ttzMass, value_ttbarMass, value_topPt,
           value_deltaPhiTtbar, value_deltaPhiTopZ, value_deltaRapTtbar, value_deltaRapTopZ;
};

#endif
