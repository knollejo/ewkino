#ifndef GenParticleHelper_H
#define GenParticleHelper_H

#include <vector>

#include "TLorentzVector.h"

#include "../../TreeReader/interface/TreeReader.h"

class GenParticleHelper {
protected:
    unsigned int* gen_n;
    double* gen_pt, * gen_eta, * gen_phi, * gen_E;
    int* gen_pdgId, * gen_charge, * gen_status;
    bool* gen_isPromptFinalState, * gen_isDirectPromptTauDecayProductFinalState, * gen_isLastCopy;
    int* gen_index, * gen_motherIndex, * gen_daughter_n;
    int(* gen_daughterIndex)[10];

public:
    // constructor
    GenParticleHelper(TreeReader& reader) :
        gen_n(&reader._gen_n),
        gen_pt(reader._gen_pt),
        gen_eta(reader._gen_eta),
        gen_phi(reader._gen_phi),
        gen_E(reader._gen_E),
        gen_pdgId(reader._gen_pdgId),
        gen_charge(reader._gen_charge),
        gen_status(reader._gen_status),
        gen_isPromptFinalState(reader._gen_isPromptFinalState),
        gen_isDirectPromptTauDecayProductFinalState(reader._gen_isDirectPromptTauDecayProductFinalState),
        gen_isLastCopy(reader._gen_isLastCopy),
        gen_index(reader._gen_index),
        gen_motherIndex(reader._gen_motherIndex),
        gen_daughter_n(reader._gen_daughter_n),
        gen_daughterIndex(reader._gen_daughterIndex)
    {}

    // access data fields
    unsigned int n() { return *gen_n; }
    double pt(int i) { return gen_pt[i]; }
    double eta(int i) { return gen_eta[i]; }
    double phi(int i) { return gen_phi[i]; }
    double E(int i) { return gen_E[i]; }
    TLorentzVector lv(int i) { TLorentzVector lv; lv.SetPtEtaPhiE(pt(i), eta(i), phi(i), E(i)); return lv; }
    int pdgId(int i) { return gen_pdgId[i]; }
    int charge(int i) { return gen_charge[i]; }
    int status(int i) { return gen_status[i]; }
    bool isPromptFinalState(int i) { return gen_isPromptFinalState[i]; }
    bool isDirectPromptTauDecayProductFinalState(int i) { return gen_isDirectPromptTauDecayProductFinalState[i]; }
    bool isLastCopy(int i) { return gen_isLastCopy[i]; }
    int index(int i) { return gen_index[i]; }
    int motherIndex(int i) { return gen_motherIndex[i]; }
    int daughter_n(int i) { return gen_daughter_n[i]; }
    int daughterIndex(int i, int j) { return gen_daughterIndex[i][j]; }

    // helpers: evaluation of particle properties
    typedef bool(GenParticleHelper::*Criterion)(int);
    bool is_any(int) { return true; }
    bool is_top(int i) { return pdgId(i)==6 || pdgId(i)==-6; }
    bool is_bottom(int i) { return pdgId(i)==5 || pdgId(i)==-5; }
    bool is_lightquark(int i) { return pdgId(i)==1 || pdgId(i)==-1 || pdgId(i)==2 || pdgId(i)==-2 || pdgId(i)==3 || pdgId(i)==-3 || pdgId(i)==4 || pdgId(i)==-4; }
    bool is_bottomorlight(int i) { return is_bottom(i) || is_lightquark(i); }
    bool is_quark(int i) { return is_top(i) || is_bottomorlight(i); }
    bool is_electron(int i) { return pdgId(i)==11 || pdgId(i)==-11; }
    bool is_muon(int i) { return pdgId(i)==13 || pdgId(i)==-13; }
    bool is_tau(int i) { return pdgId(i)==15 || pdgId(i)==-15; }
    bool is_lepton(int i) { return is_electron(i) || is_muon(i) || is_tau(i); }
    bool is_neutrino(int i) { return pdgId(i)==12 || pdgId(i)==-12 || pdgId(i)==14 || pdgId(i)==-14 || pdgId(i)==16 || pdgId(i)==-16; }
    bool is_lepton_or_neutrino(int i) { return is_lepton(i) || is_neutrino(i); }
    bool is_zboson(int i) { return pdgId(i)==23; }
    bool is_wboson(int i) { return pdgId(i)==24 || pdgId(i)==-24; }

    // helpers: routines
    int find_first_copy(int i);
    int find_last_copy(int i);
    bool is_first_copy(int i);
    bool is_last_copy(int i);

    // helpers: find sets of particles
    std::vector<int> find_particles(Criterion passes_criterion1, Criterion passes_criterion2=&GenParticleHelper::is_any);
    bool is_primal_top(int i) { return !is_top(motherIndex(i)); }
    std::vector<int> find_primal_topquarks() { return find_particles(&GenParticleHelper::is_top, &GenParticleHelper::is_primal_top); }
    std::vector<int> find_decaying_topquarks() { return find_particles(&GenParticleHelper::is_top, &GenParticleHelper::is_last_copy); }
    bool is_status62(int i) { return status(i)==62; }
    std::vector<int> find_status62_zbosons() { return find_particles(&GenParticleHelper::is_zboson, &GenParticleHelper::is_status62); }
    bool is_stable_prompt_or_tau(int i) { return status(i)==1 && (isPromptFinalState(i) || isDirectPromptTauDecayProductFinalState(i)); }
    std::vector<int> find_me_dilepton_pair();
    struct TopDecay {
        int t, b, w, w1, w2;
        TopDecay() : t(-1), b(-1), w(-1), w1(-1), w2(-1) {}
        TopDecay(int _t, int _b, int _w, int _w1, int _w2) : t(_t), b(_b), w(_w), w1(_w1), w2(_w2) {}
    };
    TopDecay find_top_decay(int i);

};

#endif
