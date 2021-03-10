#ifndef ttZTruth_H
#define ttZTruth_H

//include other parts of framework
#include "../../TreeReader/interface/TreeReader.h"

namespace ttZ {

    struct ttzTruth {
        bool isOnShell, hasLeptonicZ, hasTausInZ, hasTausInTops;
        int nLeptons;
        double zbosonPt,
               ttzMass, ttbarMass, topPt,
               deltaPhiTtbar, deltaPhiTopZ, deltaRapTtbar, deltaRapTopZ,
               topLeptonPt, topLeptonsMass, fourLeptonsMass,
               deltaPhiTopLeptons, deltaPhiTopLeptonZ, deltaRapTopLeptons, deltaRapTopLeptonZ;
        ttzTruth(bool _isOnShell, bool _hasLeptonicZ, bool _hasTausInZ, bool _hasTausInTops, int _nLeptons, double _zbosonPt, double _ttzMass, double _ttbarMass, double _topPt, double _deltaPhiTtbar, double _deltaPhiTopZ, double _deltaRapTtbar, double _deltaRapTopZ, double _topLeptonPt, double _topLeptonsMass, double _fourLeptonsMass, double _deltaPhiTopLeptons, double _deltaPhiTopLeptonZ, double _deltaRapTopLeptons, double _deltaRapTopLeptonZ) :
            isOnShell(_isOnShell), hasLeptonicZ(_hasLeptonicZ), hasTausInZ(_hasTausInZ), hasTausInTops(_hasTausInTops), nLeptons(_nLeptons), zbosonPt(_zbosonPt), ttzMass(_ttzMass), ttbarMass(_ttbarMass), topPt(_topPt), deltaPhiTtbar(_deltaPhiTtbar), deltaPhiTopZ(_deltaPhiTopZ), deltaRapTtbar(_deltaRapTtbar), deltaRapTopZ(_deltaRapTopZ), topLeptonPt(_topLeptonPt), topLeptonsMass(_topLeptonsMass), fourLeptonsMass(_fourLeptonsMass), deltaPhiTopLeptons(_deltaPhiTopLeptons), deltaPhiTopLeptonZ(_deltaPhiTopLeptonZ), deltaRapTopLeptons(_deltaRapTopLeptons), deltaRapTopLeptonZ(_deltaRapTopLeptonZ)
            {}
        ttzTruth() : isOnShell(false), hasLeptonicZ(false), hasTausInZ(false), hasTausInTops(false), nLeptons(0), zbosonPt(-999.0), ttzMass(-999.0), ttbarMass(-999.0), topPt(-999.0), deltaPhiTtbar(-999.0), deltaPhiTopZ(-999.0), deltaRapTtbar(-999.0), deltaRapTopZ(-999.0), topLeptonPt(-999.0), topLeptonsMass(-999.0), fourLeptonsMass(-999.0), deltaPhiTopLeptons(-999.0), deltaPhiTopLeptonZ(-999.0), deltaRapTopLeptons(-999.0), deltaRapTopLeptonZ(-999.0) {}
    };
    ttzTruth evaluateTruthStatus(TreeReader& reader);

}

#endif
