#ifndef ttZVariables_H
#define ttZVariables_H

//include other parts of framework
#include "../../Event/interface/Event.h"
#include "../interface/kinFitter.h"

namespace ttZ {
    struct leptonVariables {
        bool isThreeLeptons;
        int category;
        double firstLepPt, firstLepEta, firstLepPhi,
               secondLepPt, secondLepEta, secondLepPhi,
               thirdLepPt, thirdLepEta, thirdLepPhi,
               fourthLepPt, fourthLepEta, fourthLepPhi,
               dilepPt, dilepEta, dilepPhi, dilepMass;
        leptonVariables(int _category,  double _firstLepPt, double _firstLepEta, double _firstLepPhi, double _secondLepPt, double _secondLepEta, double _secondLepPhi, double _thirdLepPt, double _thirdLepEta, double _thirdLepPhi, double _dilepPt, double _dilepEta, double _dilepPhi, double _dilepMass) :
            isThreeLeptons(true), category(_category),  firstLepPt(_firstLepPt), firstLepEta(_firstLepEta), firstLepPhi(_firstLepPhi), secondLepPt(_secondLepPt), secondLepEta(_secondLepEta), secondLepPhi(_secondLepPhi), thirdLepPt(_thirdLepPt), thirdLepEta(_thirdLepEta), thirdLepPhi(_thirdLepPhi), fourthLepPt(-999.0), fourthLepEta(-999.0), fourthLepPhi(-999.0), dilepPt(_dilepPt), dilepEta(_dilepEta), dilepPhi(_dilepPhi), dilepMass(_dilepMass)
            {}
        leptonVariables(int _category,  double _firstLepPt, double _firstLepEta, double _firstLepPhi, double _secondLepPt, double _secondLepEta, double _secondLepPhi, double _thirdLepPt, double _thirdLepEta, double _thirdLepPhi, double _fourthLepPt, double _fourthLepEta, double _fourthLepPhi, double _dilepPt, double _dilepEta, double _dilepPhi, double _dilepMass) :
            isThreeLeptons(false), category(_category),  firstLepPt(_firstLepPt), firstLepEta(_firstLepEta), firstLepPhi(_firstLepPhi), secondLepPt(_secondLepPt), secondLepEta(_secondLepEta), secondLepPhi(_secondLepPhi), thirdLepPt(_thirdLepPt), thirdLepEta(_thirdLepEta), thirdLepPhi(_thirdLepPhi), fourthLepPt(_fourthLepPt), fourthLepEta(_fourthLepEta), fourthLepPhi(_fourthLepPhi), dilepPt(_dilepPt), dilepEta(_dilepEta), dilepPhi(_dilepPhi), dilepMass(_dilepMass)
            {}
        leptonVariables() :
            isThreeLeptons(false), category(-999),  firstLepPt(-999.0), firstLepEta(-999.0), firstLepPhi(-999.0), secondLepPt(-999.0), secondLepEta(-999.0), secondLepPhi(-999.0), thirdLepPt(-999.0), thirdLepEta(-999.0), thirdLepPhi(-999.0), fourthLepPt(-999.0), fourthLepEta(-999.0), fourthLepPhi(-999.0), dilepPt(-999.0), dilepEta(-999.0), dilepPhi(-999.0), dilepMass(-999.0)
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
        jetVariables(int _nJets, int _nBjets, double _missingEt, double _missingPhi, double _firstJetPt, double _firstJetEta, double _firstJetPhi, double _secondJetPt, double _secondJetEta, double _secondJetPhi) :
            nJets(_nJets), nBjets(_nBjets), missingEt(_missingEt), missingPhi(_missingPhi), firstJetPt(_firstJetPt), firstJetEta(_firstJetEta), firstJetPhi(_firstJetPhi), secondJetPt(_secondJetPt), secondJetEta(_secondJetEta), secondJetPhi(_secondJetPhi), thirdJetPt(-999.0), thirdJetEta(-999.0), thirdJetPhi(-999.0), fourthJetPt(-999.0), fourthJetEta(-999.0), fourthJetPhi(-999.0)
            {}
        jetVariables(int _nJets, int _nBjets, double _missingEt, double _missingPhi, double _firstJetPt, double _firstJetEta, double _firstJetPhi, double _secondJetPt, double _secondJetEta, double _secondJetPhi, double _thirdJetPt, double _thirdJetEta, double _thirdJetPhi) :
            nJets(_nJets), nBjets(_nBjets), missingEt(_missingEt), missingPhi(_missingPhi), firstJetPt(_firstJetPt), firstJetEta(_firstJetEta), firstJetPhi(_firstJetPhi), secondJetPt(_secondJetPt), secondJetEta(_secondJetEta), secondJetPhi(_secondJetPhi), thirdJetPt(_thirdJetPt), thirdJetEta(_thirdJetEta), thirdJetPhi(_thirdJetPhi), fourthJetPt(-999.0), fourthJetEta(-999.0), fourthJetPhi(-999.0)
            {}
        jetVariables(int _nJets, int _nBjets, double _missingEt, double _missingPhi, double _firstJetPt, double _firstJetEta, double _firstJetPhi, double _secondJetPt, double _secondJetEta, double _secondJetPhi, double _thirdJetPt, double _thirdJetEta, double _thirdJetPhi, double _fourthJetPt, double _fourthJetEta, double _fourthJetPhi) :
            nJets(_nJets), nBjets(_nBjets), missingEt(_missingEt), missingPhi(_missingPhi), firstJetPt(_firstJetPt), firstJetEta(_firstJetEta), firstJetPhi(_firstJetPhi), secondJetPt(_secondJetPt), secondJetEta(_secondJetEta), secondJetPhi(_secondJetPhi), thirdJetPt(_thirdJetPt), thirdJetEta(_thirdJetEta), thirdJetPhi(_thirdJetPhi), fourthJetPt(_fourthJetPt), fourthJetEta(_fourthJetEta), fourthJetPhi(_fourthJetPhi)
            {}
        jetVariables() :
            nJets(-999), nBjets(-999), missingEt(-999.0), missingPhi(-999.0), firstJetPt(-999.0), firstJetEta(-999.0), firstJetPhi(-999.0), secondJetPt(-999.0), secondJetEta(-999.0), secondJetPhi(-999.0), thirdJetPt(-999.0), thirdJetEta(-999.0), thirdJetPhi(-999.0), fourthJetPt(-999.0), fourthJetEta(-999.0), fourthJetPhi(-999.0)
            {}
    };
    jetVariables computeJetVariables(Event& event, const std::string& unc);

    struct reconstructedVariables {
        bool hasReco;
        double ttzMass, ttbarMass, topPt,
               deltaPhiTtbar, deltaPhiTopZ, deltaRapTtbar, deltaRapTopZ,
               lepTopMass, hadTopMass;
        reconstructedVariables(double _ttzMass, double _ttbarMass, double _topPt, double _deltaPhiTtbar, double _deltaPhiTopZ, double _deltaRapTtbar, double _deltaRapTopZ, double _lepTopMass, double _hadTopMass) :
            hasReco(true), ttzMass(_ttzMass), ttbarMass(_ttbarMass), topPt(_topPt), deltaPhiTtbar(_deltaPhiTtbar), deltaPhiTopZ(_deltaPhiTopZ), deltaRapTtbar(_deltaRapTtbar), deltaRapTopZ(_deltaRapTopZ), lepTopMass(_lepTopMass), hadTopMass(_hadTopMass)
            {}
        reconstructedVariables() :
            hasReco(false), ttzMass(-999.0), ttbarMass(-999.0), topPt(-999.0), deltaPhiTtbar(-999.0), deltaPhiTopZ(-999.0), deltaRapTtbar(-999.0), deltaRapTopZ(-999.0), lepTopMass(-999.0), hadTopMass(-999.0)
            {}
    };
    reconstructedVariables performKinematicReconstruction(Event& event, const std::string& unc, KinFitter* fitter);

    struct fourLeptonVariables {
        bool hasFourLep;
        double topLeptonPt, topLeptonsMass, fourLeptonsMass,
               deltaPhiTopLeptons, deltaPhiTopLeptonZ, deltaRapTopLeptons, deltaRapTopLeptonZ;
        fourLeptonVariables(double _topLeptonPt, double _topLeptonsMass, double _fourLeptonsMass, double _deltaPhiTopLeptons, double _deltaPhiTopLeptonZ, double _deltaRapTopLeptons, double _deltaRapTopLeptonZ) :
            hasFourLep(true), topLeptonPt(_topLeptonPt), topLeptonsMass(_topLeptonsMass), fourLeptonsMass(_fourLeptonsMass), deltaPhiTopLeptons(_deltaPhiTopLeptons), deltaPhiTopLeptonZ(_deltaPhiTopLeptonZ), deltaRapTopLeptons(_deltaRapTopLeptons), deltaRapTopLeptonZ(_deltaRapTopLeptonZ)
        {}
        fourLeptonVariables() :
            hasFourLep(false), topLeptonPt(-999.0), topLeptonsMass(-999.0), fourLeptonsMass(-999.0), deltaPhiTopLeptons(-999.0), deltaPhiTopLeptonZ(-999.0), deltaRapTopLeptons(-999.0), deltaRapTopLeptonZ(-999.0)
        {}
    };
    fourLeptonVariables performFourLeptonsComputation(Event& event);
}

#endif
