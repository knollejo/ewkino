#include "../interface/ttZVariables.h"

//include c++ library classes
#include <cmath>

//include ROOT classes
#include "TMath.h"
#include "TVector2.h"

//include other parts of framework
#include "../interface/ttZSelection.h"


#define IN_RANGE(VAL, MIN, MAX) \
    TMath::Max(MAX-0.001*(MAX-MIN), TMath::Min(MIN+0.001*(MAX-MIN), VAL))
#define PHI_IN_RANGE(VAL) \
    (isnan(VAL) ? -999.0 : IN_RANGE(TVector2::Phi_mpi_pi(VAL), -TMath::Pi(), TMath::Pi()))

std::map<std::string, double> ttZ::computeLeptonVariables(Event& event) {
    const double category = 1.0*event.numberOfElectrons();

    Lepton& lep1 = event.leptonCollection()[event.bestZBosonCandidateIndices().first];
    Lepton& lep2 = event.leptonCollection()[event.bestZBosonCandidateIndices().second];
    auto dilep = lep1+lep2;
    const double dilepPt = dilep.pt();
    const double dilepEta = dilep.eta();
    const double dilepPhi = dilep.phi();
    const double dilepMass = dilep.mass();

    std::map<std::string, double> ret = {
        {"category", category},
        {"dilepPt", IN_RANGE(dilepPt, 0.0, 500.0)},
        {"dilepEta", IN_RANGE(dilepEta, -3.0, 3.0)},
        {"dilepPhi", PHI_IN_RANGE(dilepPhi)},
        {"dilepMass", IN_RANGE(dilepMass, 81.0, 101.0)},
    };
    return ret;
}

std::map<std::string, double> ttZ::computeJetVariables(Event& event, const std::string& unc) {
    JetCollection variedJetCollection = ttZ::variedJetCollection(event, unc);
    Met variedMet = ttZ::variedMet(event, unc);

    const double nJets = 1.0*variedJetCollection.size();
    const double nBjets = 1.0*variedJetCollection.numberOfMediumBTaggedJets();

    const double missingEt = variedMet.pt();
    const double missingPhi = variedMet.phi();

    const bool has_first_jet = variedJetCollection.size()>=1;
    const double firstJetPt = has_first_jet ? variedJetCollection[0].pt() : -999.0;
    const double firstJetEta = has_first_jet ? variedJetCollection[0].eta() : -999.0;
    const double firstJetPhi = has_first_jet ? variedJetCollection[0].phi() : -999.0;

    const bool has_second_jet = variedJetCollection.size()>=2;
    const double secondJetPt = has_second_jet ? variedJetCollection[1].pt() : -999.0;
    const double secondJetEta = has_second_jet ? variedJetCollection[1].eta() : -999.0;
    const double secondJetPhi = has_second_jet ? variedJetCollection[1].phi() : -999.0;

    const bool has_third_jet = variedJetCollection.size()>=3;
    const double thirdJetPt = has_third_jet ? variedJetCollection[2].pt() : -999.0;
    const double thirdJetEta = has_third_jet ? variedJetCollection[2].eta() : -999.0;
    const double thirdJetPhi = has_third_jet ? variedJetCollection[2].phi() : -999.0;

    const bool has_fourth_jet = variedJetCollection.size()>=4;
    const double fourthJetPt = has_fourth_jet ? variedJetCollection[3].pt() : -999.0;
    const double fourthJetEta = has_fourth_jet ? variedJetCollection[3].eta() : -999.0;
    const double fourthJetPhi = has_fourth_jet ? variedJetCollection[3].phi() : -999.0;

    std::map<std::string, double> ret = {
        {"nJets", IN_RANGE(nJets, -0.5, 7.5)},
        {"nBjets", IN_RANGE(nBjets, -0.5, 3.5)},
        {"missingEt", IN_RANGE(missingEt, 0.0, 220.0)},
        {"missingPhi", PHI_IN_RANGE(missingPhi)},
        {"firstJetPt", IN_RANGE(firstJetPt, 30.0, 330.0)},
        {"firstJetEta", IN_RANGE(firstJetEta, -2.5, 2.5)},
        {"firstJetPhi", PHI_IN_RANGE(firstJetPhi)},
        {"secondJetPt", IN_RANGE(secondJetPt, 30.0, 330.0)},
        {"secondJetEta", IN_RANGE(secondJetEta, -2.5, 2.5)},
        {"secondJetPhi", PHI_IN_RANGE(secondJetPhi)},
        {"thirdJetPt", IN_RANGE(thirdJetPt, 30.0, 330.0)},
        {"thirdJetEta", IN_RANGE(thirdJetEta, -2.5, 2.5)},
        {"thirdJetPhi", PHI_IN_RANGE(thirdJetPhi)},
        {"fourthJetPt", IN_RANGE(fourthJetPt, 30.0, 330.0)},
        {"fourthJetEta", IN_RANGE(fourthJetEta, -2.5, 2.5)},
        {"fourthJetPhi", PHI_IN_RANGE(fourthJetPhi)},
    };
    return ret;
}

// std::map<std::string, double> ttZ::performKinematicReconstruction(Event& event, const std::string& unc) {
//     JetCollection variedJetCollection = ttZ::variedJetCollection(event, unc);
//     Met variedMet = ttZ::variedMet(event, unc);
//     Lepton& thirdLep = event.WLepton();
//
//
// }
