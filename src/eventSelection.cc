#include "../interface/treeReader.h"


void treeReader::orderByPt(std::vector<unsigned>& ind, const double* pt, const unsigned count){
    std::vector<std::pair<double, unsigned>> ptMap;
    for(unsigned p = 0; p < count; ++p){
        ptMap.push_back({pt[ind[p]], ind[p]});
    }
    std::sort(ptMap.begin(), ptMap.end(), [](std::pair<double, unsigned>& p1, std::pair<double, unsigned>& p2){return p1.first > p2.first;} );
    for(unsigned p = 0; p < count; ++p){
        ind[p] = ptMap[p].second;
    }
}

unsigned treeReader::dilFlavorComb(const std::vector<unsigned>& ind){
    unsigned flavCount[3] = {0,0,0};
    for(unsigned l = 0; l < 2; ++l) ++flavCount[_lFlavor[ind[l]]];
    if(flavCount[2] == 0){
        if(flavCount[1] == 0) return 0; //ee
        else if(flavCount[1] == 1) return 1; //em
        else return 2; //mm
    } else if(flavCount[2] == 1){
        if(flavCount[1] == 0) return 3; //et
        else return 4; //mt
    }
    return 5; //tt
}

double treeReader::coneCorr(const unsigned ind){
    return 1. + std::max(_relIso[ind] - 0.1, 0.);
}

void treeReader::setConePt(){
    for(unsigned l = 0; l < _nLight; ++l){
        double coneC = coneCorr(l);
        _lPt[l] *= coneC;
        _lE[l] *= coneC;
    }
}

bool treeReader::lepIsGood(const unsigned l){
    //temporary selection for dileptonCR, update later
    if(!_lEwkLoose[l]) return false;
    if(_lFlavor[l] == 0 && !_lElectronPassEmu[l]) return false;
    if(_lFlavor[l] == 1 && !_lPOGMedium[l]) return false;
    return true;
}

bool treeReader::lepIsTight(const unsigned l){
    return _lEwkTight[l];
}

unsigned treeReader::selectLep(std::vector<unsigned>& ind){
    //setConePt(); REMOVE CONE CORRECTION UNTIL MOVING TO FR
    ind.clear();
    unsigned lCount = 0;
    for(unsigned l = 0; l < _nLight; ++l){
        if(lepIsGood(l)){
            ++lCount;
            ind.push_back(l);
        }
    }
    if(lCount < 2) return 0;
    orderByPt(ind, _lPt, lCount);
    return lCount;	
}

unsigned treeReader::tightLepCount(const std::vector<unsigned>& ind, const unsigned lCount){
    unsigned tightC = 0; 
    for(unsigned l = 0; l < lCount; ++l){
        if(lepIsTight(l)) ++tightC;
        else return tightC;
    }
    return tightC;
}

bool treeReader::passPtCuts(const std::vector<unsigned>& ind){
    if(_lPt[ind[0]] < 25) return false;
    if(_lPt[ind[1]] < 15) return false;
    return true;
}

bool treeReader::jetIsClean(const unsigned ind){
    TLorentzVector jet;	
    jet.SetPtEtaPhiE(_jetPt[ind], _jetEta[ind], _jetPhi[ind], _jetE[ind]);
    for(unsigned l = 0; l < _nLight; ++l){
        if(lepIsGood(l)){
            TLorentzVector lep;
            lep.SetPtEtaPhiE(_lPt[l], _lEta[l], _lPhi[l], _lE[l]);
            if(lep.DeltaR(jet) < 0.4) return false;
        }
    }
    return true;
}

bool treeReader::jetIsGood(const unsigned ind, const unsigned ptCut, const unsigned unc, const bool clean){
    if(fabs(_jetEta[ind]) > 2.4) return false;
    switch(unc){
        case 0: if(_jetPt[ind] < ptCut) return false;
        case 1: if(_jetPt_JECDown[ind] < ptCut) return false;
        case 2: if(_jetPt_JECUp[ind] < ptCut) return false;
        case 3: if(_jetPt_JERDown[ind] < ptCut) return false;
        case 4: if(_jetPt_JERUp[ind] < ptCut) return true;
        default: ;
    }
    return !clean || jetIsClean(ind);
}

unsigned treeReader::nJets(const unsigned unc, const bool clean){
    unsigned nJets = 0;
    for(unsigned j = 0; j < _nJets; ++j){
        if(jetIsGood(j, 30, unc, clean)) ++nJets;
    }
    return nJets;
}

unsigned treeReader::nJets(std::vector<unsigned>& jetInd, const unsigned unc, const bool clean){
    unsigned nJets = 0;
    jetInd.clear();
    for(unsigned j = 0; j < _nJets; ++j){
        if(jetIsGood(j, 30, unc, clean)){
            ++nJets;
            jetInd.push_back(nJets);
        }
    }
    orderByPt(jetInd, _jetPt, nJets);
    return nJets;
}

bool treeReader::bTaggedDeepCSV(const unsigned ind, const unsigned wp){
    static const double bTagWP[3] = {0.2219, 0.6324,  0.8958};
    return (_jetDeepCsv_b[ind] + _jetDeepCsv_bb[ind]) > bTagWP[wp];
}

bool treeReader::bTaggedCSVv2(const unsigned ind, const unsigned wp){
    static const double bTagWP[3] = {0.5426, 0.8484, 0.9535};
    return _jetCsvV2[ind] > bTagWP[wp];
}

unsigned treeReader::nBJets(const unsigned unc, const bool deepCSV, const bool clean, const unsigned wp){
    unsigned nbJets = 0;
    for(unsigned j = 0; j < _nJets; ++j){
        if(jetIsGood(j, 25, unc, clean)){
            if(deepCSV && bTaggedDeepCSV(j, wp)) ++nbJets;
            else if(bTaggedCSVv2(j, wp)) ++nbJets;
        }
    }
    return nbJets;
}
