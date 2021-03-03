#ifndef ttZTruth_H
#define ttZTruth_H

#include "TLorentzVector.h"

//include other parts of framework
#include "../../TreeReader/interface/TreeReader.h"

namespace ttZ {

    struct ttzTruth {
        bool isOnShell, hasLeptonicZ, hasTaus;
        int nLeptons;
        TLorentzVector lvZ, lvTop, lvAntitop;
        ttzTruth(bool _isOnShell, bool _hasLeptonicZ, bool _hasTaus, int _nLeptons, TLorentzVector _lvZ, TLorentzVector _lvTop, TLorentzVector _lvAntitop) :
            isOnShell(_isOnShell), hasLeptonicZ(_hasLeptonicZ), hasTaus(_hasTaus), nLeptons(_nLeptons), lvZ(_lvZ), lvTop(_lvTop), lvAntitop(_lvAntitop)
            {}
        ttzTruth() : isOnShell(false), hasLeptonicZ(false), hasTaus(false), nLeptons(false) {}
    };
    ttzTruth evaluateTruthStatus(TreeReader& reader);

}

#endif
