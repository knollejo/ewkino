#include "../interface/ttZTruth.h"

#include <cmath>
#include <vector>

#include "TLorentzVector.h"

#include "../../constants/particleMasses.h"
#include "../interface/GenParticleHelper.h"
#include "../interface/ttZObservables.h"

ttZ::ttzTruth ttZ::evaluateTruthStatus(TreeReader& reader) {
    // initialize GenParticleHelper
    GenParticleHelper geninfo(reader);

    // find gen-level Z boson
    std::vector<int> zleptons = geninfo.find_me_dilepton_pair();
    const bool has_zboson = (zleptons.size()==2);
    const bool zboson_decays_leptonic = has_zboson ? geninfo.is_lepton(zleptons.at(0)) : false;
    const bool zboson_has_taus = zboson_decays_leptonic ? geninfo.is_tau(zleptons.at(0)) || geninfo.is_tau(zleptons.at(1)) : false;
    TLorentzVector lvZlep1 = has_zboson ? geninfo.lv(zleptons.at(geninfo.pdgId(zleptons.at(0))>0 ? 0 : 1)) : TLorentzVector();
    TLorentzVector lvZlep2 = has_zboson ? geninfo.lv(zleptons.at(geninfo.pdgId(zleptons.at(0))>0 ? 1 : 0)) : TLorentzVector();
    TLorentzVector lvZ = lvZlep1+lvZlep2;
    const bool is_onshell_zboson = has_zboson ? std::abs(lvZ.M()-particle::mZ)<=15.0 : false;

    // find gen-level top quarks
    std::vector<int> topquarks = geninfo.find_decaying_topquarks();
    const bool has_topquarks = (topquarks.size()==2);
    GenParticleHelper::TopDecay topquark = has_topquarks ? geninfo.find_top_decay(topquarks.at(geninfo.pdgId(topquarks.at(0))>0 ? 0 : 1)) : GenParticleHelper::TopDecay();
    GenParticleHelper::TopDecay topantiquark = has_topquarks ? geninfo.find_top_decay(topquarks.at(geninfo.pdgId(topquarks.at(0))>0 ? 1 : 0)) : GenParticleHelper::TopDecay();
    const bool topquarks_have_taus = has_topquarks ? geninfo.is_tau(topquark.w2) || geninfo.is_tau(topantiquark.w2) : false;
    const bool topquarks_are_dileptonic = has_topquarks ? geninfo.is_lepton(topquark.w2) && geninfo.is_lepton(topantiquark.w2) : false;
    const bool topquarks_are_semileptonic = has_topquarks ? geninfo.is_lepton(topquark.w2) != geninfo.is_lepton(topantiquark.w2) : false;
    TLorentzVector lvTop = has_topquarks ? geninfo.lv(topquark.t) : TLorentzVector();
    TLorentzVector lvAntitop = has_topquarks ? geninfo.lv(topantiquark.t) : TLorentzVector();
    TLorentzVector lvTopLepton = topquarks_are_dileptonic ? geninfo.lv(topquark.w2) : TLorentzVector();
    TLorentzVector lvAntitopLepton = topquarks_are_dileptonic ? geninfo.lv(topantiquark.w2) : TLorentzVector();

    // compute observables
    const double zbosonPt = lvZ.Pt();
    const double ttzMass = (lvZ+lvTop+lvAntitop).M();
    const double ttbarMass = (lvTop+lvAntitop).M();
    const double topPt = lvTop.Pt();
    const double deltaPhiTtbar = deltaPhi(lvTop.Phi(), lvAntitop.Phi());
    const double deltaPhiTopZ = deltaPhi(lvTop.Phi(), lvZ.Phi());
    const double deltaRapTtbar = deltaRap(lvTop.Rapidity(), lvAntitop.Rapidity());
    const double deltaRapTopZ = deltaRap(lvTop.Rapidity(), lvZ.Rapidity());
    const double topLeptonPt = topquarks_are_dileptonic ? lvTopLepton.Pt() : -999.0;
    const double topLeptonsMass = topquarks_are_dileptonic ? (lvTopLepton+lvAntitopLepton).M() : -999.0;
    const double fourLeptonsMass = topquarks_are_dileptonic ? (lvTopLepton+lvAntitopLepton+lvZ).Pt() : -999.0;
    const double deltaPhiTopLeptons = topquarks_are_dileptonic ? deltaPhi(lvTopLepton.Phi(), lvAntitopLepton.Phi()) : -999.0;
    const double deltaPhiTopLeptonZ = topquarks_are_dileptonic ? deltaPhi(lvTopLepton.Phi(), lvZ.Phi()) : -999.0;
    const double deltaRapTopLeptons = topquarks_are_dileptonic ? deltaRap(lvTopLepton.Rapidity(), lvAntitopLepton.Rapidity()) : -999.0;
    const double deltaRapTopLeptonZ = topquarks_are_dileptonic ? deltaRap(lvTopLepton.Rapidity(), lvZ.Rapidity()) : -999.0;

    // prepare return value
    return ttzTruth(
        is_onshell_zboson,
        zboson_decays_leptonic,
        zboson_has_taus,
        topquarks_have_taus,
        (zboson_decays_leptonic ? 2 : 0)+(topquarks_are_dileptonic ? 2 : topquarks_are_semileptonic ? 1 : 0),
        zbosonPt,
        ttzMass, ttbarMass, topPt,
        deltaPhiTtbar, deltaPhiTopZ, deltaRapTtbar, deltaRapTopZ,
        topLeptonPt, topLeptonsMass, fourLeptonsMass,
        deltaPhiTopLeptons, deltaPhiTopLeptonZ, deltaRapTopLeptons, deltaRapTopLeptonZ
    );
}
