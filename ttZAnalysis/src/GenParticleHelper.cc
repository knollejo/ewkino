#include "../interface/GenParticleHelper.h"

#include <unordered_set>

int GenParticleHelper::find_first_copy(int i) {
    while(pdgId(i)==pdgId(motherIndex(i))) {
        i = motherIndex(i);
    }
    return i;
}

int GenParticleHelper::find_last_copy(int i) {
    bool found;
    do {
        found = false;
        for(int j=0; j<daughter_n(i); j++) {
            const int daughter = daughterIndex(i, j);
            if(pdgId(i)==pdgId(daughter)) {
                i = daughter;
                found = true;
                break;
            }
        }
    } while(found);
    return i;
}

bool GenParticleHelper::is_first_copy(int i) {
    return pdgId(i)!=pdgId(motherIndex(i));
}

bool GenParticleHelper::is_last_copy(int i) {
    for(int j=0; j<daughter_n(i); j++) {
        if(pdgId(i)==pdgId(daughterIndex(i, j))) return false;
    }
    return true;
}

std::vector<int> GenParticleHelper::find_particles(GenParticleHelper::Criterion passes_criterion1, GenParticleHelper::Criterion passes_criterion2) {
    std::vector<int> found;
    for(unsigned int i=0; i<n(); i++) {
        if((this->*passes_criterion1)(i) && (this->*passes_criterion2)(i)) found.push_back(i);
    }
    return found;
}

std::vector<int> GenParticleHelper::find_me_dilepton_pair() {
    std::vector<int> candidates = find_particles(&GenParticleHelper::is_lepton_or_neutrino, &GenParticleHelper::is_stable_prompt_or_tau);
    std::unordered_set<int> found;
    for(auto const& candidate: candidates) {
        int first = find_first_copy(candidate);
        if(isDirectPromptTauDecayProductFinalState(first)) first = find_first_copy(motherIndex(first));
        if(!is_wboson(motherIndex(first))) found.insert(first);
    }
    return std::vector<int>(found.begin(), found.end());
}

GenParticleHelper::TopDecay GenParticleHelper::find_top_decay(int i) {
    const int topquark = find_last_copy(i);
    int wboson = -1, bquark = -1;
    for(int j=0; j<daughter_n(topquark); j++) {
        const int daughter = daughterIndex(topquark, j);
        if(is_wboson(daughter)) wboson = daughter;
        else if(is_bottom(daughter)) bquark = daughter;
    }
    if(wboson<0) return TopDecay(topquark, bquark, -1, -1, -1);
    const int wid = pdgId(wboson);
    int wfordaughters = find_last_copy(wboson);
    int wdaughter1 = -1, wdaughter2 = -1;
    for(int j=0; j<daughter_n(wfordaughters); j++) {
        const int daughter = daughterIndex(wfordaughters, j);
        if(!is_lepton_or_neutrino(daughter) && !is_quark(daughter)) continue;
        if(pdgId(daughter)*wid>0) wdaughter1 = daughter;
        else wdaughter2 = daughter;
    }
    return TopDecay(topquark, bquark, wboson, wdaughter1, wdaughter2);
}
