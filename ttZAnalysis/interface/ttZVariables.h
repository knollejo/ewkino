#ifndef ttZVariables_H
#define ttZVariables_H

//include other parts of framework
#include "../../Event/interface/Event.h"
#include "../interface/kinFitter.h"

namespace ttZ {
    struct leptonVariables {
        int category;
        double leadLepPt, leadLepEta, leadLepPhi,
               sublLepPt, sublLepEta, sublLepPhi,
               trailLepPt, trailLepEta, trailLepPhi,
               dilepPt, dilepEta, dilepPhi, dilepMass;
        leptonVariables(int _category,  double _leadLepPt, double _leadLepEta, double _leadLepPhi, double _sublLepPt, double _sublLepEta, double _sublLepPhi, double _trailLepPt, double _trailLepEta, double _trailLepPhi, double _dilepPt, double _dilepEta, double _dilepPhi, double _dilepMass) :
            category(_category),  leadLepPt(_leadLepPt), leadLepEta(_leadLepEta), leadLepPhi(_leadLepPhi), sublLepPt(_sublLepPt), sublLepEta(_sublLepEta), sublLepPhi(_sublLepPhi), trailLepPt(_trailLepPt), trailLepEta(_trailLepEta), trailLepPhi(_trailLepPhi), dilepPt(_dilepPt), dilepEta(_dilepEta), dilepPhi(_dilepPhi), dilepMass(_dilepMass)
            {}
    };
    leptonVariables computeLeptonVariables(Event& event);

    struct jetVariables {
        int nJets, nBjets;
        double missingEt, missingPhi,
               firstJetPt, firstJetEta, firstJetPhi,
               secondJetPt, secondJetEta, secondJetPhi,
               thirdJetPt, thirdJetEta, thirdJetPhi,
               fourthJetPt, fourthJetEta, fourthJetPhi;
        jetVariables(int _nJets, int _nBjets, double _missingEt, double _missingPhi, double _firstJetPt, double _firstJetEta, double _firstJetPhi, double _secondJetPt, double _secondJetEta, double _secondJetPhi, double _thirdJetPt, double _thirdJetEta, double _thirdJetPhi, double _fourthJetPt, double _fourthJetEta, double _fourthJetPhi) :
            nJets(_nJets), nBjets(_nBjets), missingEt(_missingEt), missingPhi(_missingPhi), firstJetPt(_firstJetPt), firstJetEta(_firstJetEta), firstJetPhi(_firstJetPhi), secondJetPt(_secondJetPt), secondJetEta(_secondJetEta), secondJetPhi(_secondJetPhi), thirdJetPt(_thirdJetPt), thirdJetEta(_thirdJetEta), thirdJetPhi(_thirdJetPhi), fourthJetPt(_fourthJetPt), fourthJetEta(_fourthJetEta), fourthJetPhi(_fourthJetPhi)
            {}
    };
    jetVariables computeJetVariables(Event& event, const std::string& unc);

    struct reconstructedVariables {
        double ttzMass, ttbarMass, topPt,
               deltaPhiTtbar, deltaPhiTopZ, deltaRapTtbar, deltaRapTopZ;
        reconstructedVariables(double _ttzMass, double _ttbarMass, double _topPt, double _deltaPhiTtbar, double _deltaPhiTopZ, double _deltaRapTtbar, double _deltaRapTopZ) :
            ttzMass(_ttzMass), ttbarMass(_ttbarMass), topPt(_topPt), deltaPhiTtbar(_deltaPhiTtbar), deltaPhiTopZ(_deltaPhiTopZ), deltaRapTtbar(_deltaRapTtbar), deltaRapTopZ(_deltaRapTopZ)
            {}
    };
    reconstructedVariables performKinematicReconstruction(Event& event, const std::string& unc, KinFitter* fitter);
}

#endif
