#include "../interface/ttZVariables.h"

//include c++ library classes
#include <algorithm>
#include <cmath>

//include ROOT classes
#include "TLorentzVector.h"
#include "TMath.h"
#include "TVector2.h"

//include other parts of framework
#include "../interface/ttZSelection.h"

ttZ::leptonVariables ttZ::computeLeptonVariables(Event& event) {
    const int category = event.numberOfElectrons();

    const double leadLepPt = event.leptonCollection()[0].pt();
    const double leadLepEta = event.leptonCollection()[0].eta();
    const double leadLepPhi = event.leptonCollection()[0].phi();
    const double sublLepPt = event.leptonCollection()[1].pt();
    const double sublLepEta = event.leptonCollection()[1].eta();
    const double sublLepPhi = event.leptonCollection()[1].phi();
    const double trailLepPt = event.leptonCollection()[2].pt();
    const double trailLepEta = event.leptonCollection()[2].eta();
    const double trailLepPhi = event.leptonCollection()[2].phi();

    Lepton& lep1 = event.leptonCollection()[event.bestZBosonCandidateIndices().first];
    Lepton& lep2 = event.leptonCollection()[event.bestZBosonCandidateIndices().second];
    auto dilep = lep1+lep2;
    const double dilepPt = dilep.pt();
    const double dilepEta = dilep.eta();
    const double dilepPhi = dilep.phi();
    const double dilepMass = dilep.mass();

    return leptonVariables(
        category,
        leadLepPt, leadLepEta, leadLepPhi,
        sublLepPt, sublLepEta, sublLepPhi,
        trailLepPt, trailLepEta, trailLepPhi,
        dilepPt, dilepEta, dilepPhi, dilepMass
    );
}

ttZ::jetVariables ttZ::computeJetVariables(Event& event, const std::string& unc) {
    JetCollection variedJetCollection = ttZ::variedJetCollection(event, unc);
    Met variedMet = ttZ::variedMet(event, unc);

    const int nJets = variedJetCollection.size();
    const int nBjets = variedJetCollection.numberOfMediumBTaggedJets();

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

    return jetVariables(
        nJets, nBjets,
        missingEt, missingPhi,
        firstJetPt, firstJetEta, firstJetPhi,
        secondJetPt, secondJetEta, secondJetPhi,
        thirdJetPt, thirdJetEta, thirdJetPhi,
        fourthJetPt, fourthJetEta, fourthJetPhi
    );
}

struct FitResult {
    KinFitter::Result result;
    unsigned int iJets[4];
    TLorentzVector lvLepton, lvNeutrino, lvJets[4];
};

double deltaPhi(double phi1, double phi2) {
    const double dphi = phi1-phi2;
    if(isnan(dphi)) return -999.0;
    else return TMath::Abs(TVector2::Phi_mpi_pi(dphi));
}

double deltaRap(double rap1, double rap2) {
    return TMath::Abs(rap1-rap2);
}

ttZ::reconstructedVariables ttZ::performKinematicReconstruction(Event& event, const std::string& unc, KinFitter* fitter) {
    // prepare third lepton
    Lepton& thirdLep = event.WLepton();
    TLorentzVector lvLepton;
    lvLepton.SetPtEtaPhiM(thirdLep.pt(), thirdLep.eta(), thirdLep.phi(), 0.0);
    const bool leptonIsMuon = thirdLep.isMuon();

    // prepare neutrino candidates
    Met variedMet = ttZ::variedMet(event, unc);
    TLorentzVector lvMet;
    lvMet.SetXYZM(variedMet.px(), variedMet.py(), 0.0, 0.0);
    double neutrinoSolutions[2];
    const unsigned int nNeutrinos = fitter->findNeutrinoSolutions(lvLepton, lvMet, neutrinoSolutions[0], neutrinoSolutions[1]);
    std::vector<TLorentzVector> neutrinos;
    if(nNeutrinos==0) {
        neutrinos.push_back(lvMet);
    } else for(unsigned int iNeutrino=0; iNeutrino<nNeutrinos; iNeutrino++) {
        TLorentzVector neutrino;
        neutrino.SetXYZM(lvMet.X(), lvMet.Y(), neutrinoSolutions[iNeutrino], 0.0);
        neutrinos.push_back(neutrino);
    }

    // prepare jets
    JetCollection variedJetCollection = ttZ::variedJetCollection(event, unc);
    variedJetCollection.sortByPt();
    const unsigned int nJets = variedJetCollection.size();
    std::vector<TLorentzVector> jets, bjets;
    std::vector<unsigned int> iBjets;
    for(unsigned int iJet=0; iJet<nJets; iJet++) {
        Jet& jet = variedJetCollection[iJet];
        TLorentzVector lvJet, lvBjet;
        lvJet.SetPtEtaPhiM(jet.pt(), jet.eta(), jet.phi(), 0.0);
        lvBjet.SetPtEtaPhiM(jet.pt(), jet.eta(), jet.phi(), fitter->BottomMass());
        jets.push_back(lvJet);
        bjets.push_back(lvBjet);
        if(jet.isBTaggedMedium()) iBjets.push_back(iJet);
    }
    const unsigned int nBjets = iBjets.size();

    // loop over all possible combinations
    std::vector<FitResult> results;
    // loop over neutrino solutions
    for(auto const& lvNeutrino: neutrinos) {
        // loop over first b quark candidates
        for(auto const& iBjet1: iBjets) {
            TLorentzVector lvBjet1 = bjets.at(iBjet1);
            std::vector<unsigned int> iBjet2Candidates;
            if(nBjets>=2) for(auto const& iBjet2: iBjets) { if(iBjet2!=iBjet1) iBjet2Candidates.push_back(iBjet2); }
            else for(unsigned int iBjet2=0; iBjet2<nJets; iBjet2++) { if(iBjet2!=iBjet1) iBjet2Candidates.push_back(iBjet2); }
            // loop over second b quark candidates
            for(auto const& iBjet2: iBjet2Candidates) {
                TLorentzVector lvBjet2 = bjets.at(iBjet2);
                // loop over first light quark candidates
                for(unsigned int iLight1=0; iLight1<nJets; iLight1++) {
                    if(iLight1==iBjet1 || iLight1==iBjet2) continue;
                    TLorentzVector lvLight1 = jets.at(iLight1);
                    // reconstruction with three jets
                    if(nJets==3) {
                        // initialize particles
                        std::vector<KinFitter::Particle*> particles;
                        particles.push_back(leptonIsMuon ? fitter->Muon(lvLepton) : fitter->Electron(lvLepton));
                        particles.push_back(fitter->Missing(lvNeutrino));
                        particles.push_back(fitter->Bjet(lvBjet1));
                        particles.push_back(fitter->Bjet(lvBjet2));
                        particles.push_back(fitter->Jet(lvLight1));

                        // initialize constraints
                        std::vector<KinFitter::MassConstraint*> constraints;
                        constraints.push_back(fitter->WMassConstraint(particles.at(0), particles.at(1)));
                        constraints.push_back(fitter->TopMassConstraint(particles.at(0), particles.at(1), particles.at(2)));
                        constraints.push_back(fitter->TopMassConstraint(particles.at(3), particles.at(4)));

                        // perform the fit
                        KinFitter::Result result = *fitter->fit(particles, constraints);

                        // store the fit results
                        FitResult res;
                        res.result = result;
                        res.iJets[0] = iBjet1;
                        res.iJets[1] = iBjet2;
                        res.iJets[2] = iLight1;
                        res.lvLepton = particles.at(0)->vec;
                        res.lvNeutrino = particles.at(1)->vec;
                        for(unsigned int i=0; i<3; i++) res.lvJets[i] = particles.at(2+i)->vec;
                        results.push_back(res);
                    }
                    // reconstruction with four jets
                    // loop over second light quark candidates
                    else for(unsigned int iLight2=iLight1+1; iLight2<nJets; iLight2++) {
                        TLorentzVector lvLight2 = jets.at(iLight2);

                        // initialize particles
                        std::vector<KinFitter::Particle*> particles;
                        particles.push_back(leptonIsMuon ? fitter->Muon(lvLepton) : fitter->Electron(lvLepton));
                        particles.push_back(fitter->Missing(lvNeutrino));
                        particles.push_back(fitter->Bjet(lvBjet1));
                        particles.push_back(fitter->Bjet(lvBjet2));
                        particles.push_back(fitter->Jet(lvLight1));
                        particles.push_back(fitter->Jet(lvLight2));

                        // initialize constraints
                        std::vector<KinFitter::MassConstraint*> constraints;
                        constraints.push_back(fitter->WMassConstraint(particles.at(0), particles.at(1)));
                        constraints.push_back(fitter->WMassConstraint(particles.at(4), particles.at(5)));
                        constraints.push_back(fitter->TopMassConstraint(particles.at(0), particles.at(1), particles.at(2)));
                        constraints.push_back(fitter->TopMassConstraint(particles.at(3), particles.at(4), particles.at(5)));

                        // perform the fit
                        KinFitter::Result result = *fitter->fit(particles, constraints);

                        // store the fit results
                        FitResult res;
                        res.result = result;
                        res.iJets[0] = iBjet1;
                        res.iJets[1] = iBjet2;
                        res.iJets[2] = iLight1;
                        res.iJets[3] = iLight2;
                        res.lvLepton = particles.at(0)->vec;
                        res.lvNeutrino = particles.at(1)->vec;
                        for(unsigned int i=0; i<4; i++) res.lvJets[i] = particles.at(2+i)->vec;
                        results.push_back(res);
                    }
                }
            }
        }
    }

    // select best fit result
    std::stable_sort(results.begin(), results.end(),
        [](const FitResult& a, const FitResult& b) {
            if(a.result.status==0 && b.result.status!=0) return true;
            else if(a.result.status!=0 && b.result.status==0) return false;
            else return a.result.s/a.result.ndf<b.result.s/b.result.ndf;
        }
    );
    FitResult bestfit = results.at(0);

    // reconstruct Z boson
    Lepton& lep1 = event.leptonCollection()[event.bestZBosonCandidateIndices().first];
    Lepton& lep2 = event.leptonCollection()[event.bestZBosonCandidateIndices().second];
    auto dilep = lep1+lep2;
    TLorentzVector zboson;
    zboson.SetPtEtaPhiM(dilep.pt(), dilep.eta(), dilep.phi(), dilep.mass());

    // reconstruct tops
    TLorentzVector leptop = bestfit.lvLepton+bestfit.lvNeutrino+bestfit.lvJets[0];
    TLorentzVector hadtop = bestfit.lvJets[1]+bestfit.lvJets[2]+bestfit.lvJets[3];
    TLorentzVector ttbar = leptop+hadtop;
    TLorentzVector ttz = ttbar+zboson;
    TLorentzVector top = thirdLep.charge()>0 ? leptop : hadtop;
    TLorentzVector antitop = thirdLep.charge()>0 ? hadtop : leptop;

    // calculate observables
    const double ttzMass = ttz.M();
    const double ttbarMass = ttbar.M();
    const double topPt = top.Pt();
    const double deltaPhiTtbar = deltaPhi(top.Phi(), antitop.Phi());
    const double deltaPhiTopZ = deltaPhi(top.Phi(), zboson.Phi());
    const double deltaRapTtbar = deltaRap(top.Rapidity(), antitop.Rapidity());
    const double deltaRapTopZ = deltaRap(top.Rapidity(), zboson.Rapidity());

    return reconstructedVariables(
        ttzMass, ttbarMass, topPt,
        deltaPhiTtbar, deltaPhiTopZ, deltaRapTtbar, deltaRapTopZ
    );
}
