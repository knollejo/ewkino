#ifndef PlotAxes_H
#define PlotAxes_H

#include <map>
#include <string>

const std::string LABEL_GEV = "[GeV]";
const std::string LABEL_PT = "p_{T} "+LABEL_GEV;
const std::string LABEL_ETA = "η";
const std::string LABEL_PHI = "φ";
const std::string LABEL_MASS = "mass "+LABEL_GEV;

const std::map<std::string, std::string> XAXIS_LABELS({
    {"category3l", "flavor category"},
    {"category4l", "flavor category"},
    {"nJets", "jet multiplicity"},
    {"nBjets", "b-jet multiplicity"},
    {"dilepPt", "dilepton "+LABEL_PT},
    {"dilepEta", "dilepton "+LABEL_ETA},
    {"dilepPhi", "dilepton "+LABEL_PHI},
    {"dilepMass3l", "dilepton mass "+LABEL_MASS},
    {"dilepMass4l", "dilepton mass "+LABEL_MASS},
    {"missingEt", "missing momentum "+LABEL_GEV},
    {"missingPhi", "missing momentum "+LABEL_PHI},
    {"firstLepPt", "1st lepton "+LABEL_PT},
    {"firstLepEta", "1st lepton "+LABEL_ETA},
    {"firstLepPhi", "1st lepton"+LABEL_PHI},
    {"secondLepPt", "2nd lepton "+LABEL_PT},
    {"secondLepEta", "2nd lepton "+LABEL_ETA},
    {"secondLepPhi", "2nd lepton"+LABEL_PHI},
    {"thirdLepPt", "3rd lepton "+LABEL_PT},
    {"thirdLepEta", "3rd lepton "+LABEL_ETA},
    {"thirdLepPhi", "3rd lepton"+LABEL_PHI},
    {"fourthLepPt", "4th lepton "+LABEL_PT},
    {"fourthLepEta", "4th lepton "+LABEL_ETA},
    {"fourthLepPhi", "4th lepton"+LABEL_PHI},
    {"firstJetPt", "1st jet "+LABEL_PT},
    {"firstJetEta", "1st jet "+LABEL_ETA},
    {"firstJetPhi", "1st jet"+LABEL_PHI},
    {"secondJetPt", "2nd jet "+LABEL_PT},
    {"secondJetEta", "2nd jet "+LABEL_ETA},
    {"secondJetPhi", "2nd jet"+LABEL_PHI},
    {"thirdJetPt", "3rd jet "+LABEL_PT},
    {"thirdJetEta", "3rd jet "+LABEL_ETA},
    {"thirdJetPhi", "3rd jet"+LABEL_PHI},
    {"fourthJetPt", "4th jet "+LABEL_PT},
    {"fourthJetEta", "4th jet "+LABEL_ETA},
    {"fourthJetPhi", "4th jet"+LABEL_PHI},
    {"zbosonPt", "p_{T}(Z) "+LABEL_GEV},
    {"ttzMass", "m(ttZ) "+LABEL_GEV},
    {"ttbarMass", "m(tt) "+LABEL_GEV},
    {"topPt", "p_{T}(t) "+LABEL_GEV},
    {"deltaPhiTtbar", "|Δφ(t,t)|"},
    {"deltaPhiTopZ", "|Δφ(t,Z)|"},
    {"deltaRapTtbar", "|Δy(t,t)|"},
    {"deltaRapTopZ", "|Δy(t,Z)|"},
    {"lepTopMass", "m(leptonic t) "+LABEL_GEV},
    {"hadTopMass", "m(hadronic t) "+LABEL_GEV},
    {"topLeptonPt", "p_{T}(ℓ_{t}) "+LABEL_GEV},
    {"topLeptonsMass", "m(ℓ_{t}ℓ_{t}) "+LABEL_GEV},
    {"fourLeptonsMass", "m(4ℓ) "+LABEL_GEV},
    {"deltaPhiTopLeptons", "|Δφ(ℓ_{t},ℓ_{t})|"},
    {"deltaPhiTopLeptonZ", "|Δφ(ℓ_{t},Z)|"},
    {"deltaRapTopLeptons", "|Δy(ℓ_{t},ℓ_{t})|"},
    {"deltaRapTopLeptonZ", "|Δy(ℓ_{t},Z)|"},
});
std::string get_xaxis_label(std::string obs) {
    auto found = XAXIS_LABELS.find(obs);
    return (found==XAXIS_LABELS.end()) ? obs : found->second;
}

const std::string XAXIS_ADD_FIRST_PHI = "set xtics 0.5*pi;\\\n";
const std::map<std::string, std::string> XAXIS_ADD_FIRST({
    {"category3l", (
        "set xtics -0.5,1;\\\n"
        "set mxtics 1;\\\n"
    )},
    {"category3l", (
        "set xtics -0.5,1;\\\n"
        "set mxtics 1;\\\n"
    )},
    {"nJets", (
        "set xtics 2.5,1;\\\n"
        "set mxtics 1;\\\n"
    )},
    {"nBjets", (
        "set xtics 0.5,1;\\\n"
        "set mxtics 1;\\\n"
    )},
    {"dilepPhi", XAXIS_ADD_FIRST_PHI},
    {"missingPhi", XAXIS_ADD_FIRST_PHI},
    {"firstLepPhi", XAXIS_ADD_FIRST_PHI},
    {"secondLepPhi", XAXIS_ADD_FIRST_PHI},
    {"thirdLepPhi", XAXIS_ADD_FIRST_PHI},
    {"fourthLepPhi", XAXIS_ADD_FIRST_PHI},
    {"firstJetPhi", XAXIS_ADD_FIRST_PHI},
    {"secondJetPhi", XAXIS_ADD_FIRST_PHI},
    {"thirdJetPhi", XAXIS_ADD_FIRST_PHI},
    {"fourthJetPhi", XAXIS_ADD_FIRST_PHI},
    {"deltaPhiTopLeptons": "set xtics 0.25*pi;\\\n"},
    {"deltaPhiTopLeptonZ": "set xtics 0.25*pi;\\\n"},
});
std::string get_xaxis_add_first(std::string obs) {
    auto found = XAXIS_ADD_FIRST.find(obs);
    return (found==XAXIS_ADD_FIRST.end()) ? "" : found->second;
}

const std::string XAXIS_ADD_SECOND_PHI = "set xtics add ('–π' -pi, '–π/2' -0.5*pi, 'π/2' 0.5*pi, 'π' pi);\\\n";
const std::map<std::string, std::string> XAXIS_ADD_SECOND({
    {"category3l", (
        "set xtics format '';\\\n"
        "set label 'µµµ' center at 0,screen 0.09 font ',20pt';\\\n"
        "set label 'µµe' center at 1,screen 0.09 font ',20pt';\\\n"
        "set label 'µee' center at 2,screen 0.09 font ',20pt';\\\n"
        "set label 'eee' center at 3,screen 0.09 font ',20pt';\\\n"
    )},
    {"category4l", (
        "set xtics format '';\\\n"
        "set label '(µµ)µµ' center at 0,screen 0.1 font ',20pt';\\\n"
        "set label '(µµ)µe' center at 1,screen 0.1 font ',20pt';\\\n"
        "set label '(µµ)ee' center at 2,screen 0.1 font ',20pt';\\\n"
        "set label '(ee)µµ' center at 3,screen 0.1 font ',20pt';\\\n"
        "set label '(ee)µe' center at 4,screen 0.1 font ',20pt';\\\n"
        "set label '(ee)ee' center at 5,screen 0.1 font ',20pt';\\\n"
    )},
    {"nJets", (
        "set xtics format '';\\\n"
        "set label '2' center at 2,screen 0.1 font ',20pt';\\\n"
        "set label '3' center at 3,screen 0.1 font ',20pt';\\\n"
        "set label '4' center at 4,screen 0.1 font ',20pt';\\\n"
        "set label '5' center at 5,screen 0.1 font ',20pt';\\\n"
        "set label '6' center at 6,screen 0.1 font ',20pt';\\\n"
        "set label '≥7' center at 7,screen 0.1 font ',20pt';\\\n"
    )},
    {"nJets", (
        "set xtics format '';\\\n"
        "set label '1' center at 1,screen 0.1 font ',20pt';\\\n"
        "set label '2' center at 2,screen 0.1 font ',20pt';\\\n"
        "set label '≥3' center at 3,screen 0.1 font ',20pt';\\\n"
    )},
    {"dilepPhi", XAXIS_ADD_SECOND_PHI},
    {"missingPhi", XAXIS_ADD_SECOND_PHI},
    {"firstLepPhi", XAXIS_ADD_SECOND_PHI},
    {"secondLepPhi", XAXIS_ADD_SECOND_PHI},
    {"thirdLepPhi", XAXIS_ADD_SECOND_PHI},
    {"fourthLepPhi", XAXIS_ADD_SECOND_PHI},
    {"firstJetPhi", XAXIS_ADD_SECOND_PHI},
    {"secondJetPhi", XAXIS_ADD_SECOND_PHI},
    {"thirdJetPhi", XAXIS_ADD_SECOND_PHI},
    {"fourthJetPhi", XAXIS_ADD_SECOND_PHI},
    {"deltaPhiTopLeptons", "set xtics add ('π/4' 0.25*pi, 'π/2' 0.5*pi, '3π/4' 0.75*pi, 'π' pi);\\\n"},
    {"deltaPhiTopLeptonZ", "set xtics add ('π/4' 0.25*pi, 'π/2' 0.5*pi, '3π/4' 0.75*pi, 'π' pi);\\\n"},
});
std::string get_xaxis_add_second(std::string obs) {
    auto found = XAXIS_ADD_SECOND.find(obs);
    return (found==XAXIS_ADD_SECOND.end()) ? "" : found->second;
}

#endif
