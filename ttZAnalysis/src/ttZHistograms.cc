#include "../interface/ttZHistograms.h"

//include ROOT classes
#include "TMath.h"
#include "TVector2.h"

//include other parts of framework
#include "../../constants/particleMasses.h"
#include "../interface/ttZObservables.h"

// histogram definitions
const int NBINS_category3l = 4; const double MIN_category3l = -0.5; const double MAX_category3l = 3.5;
const int NBINS_category4l = 6; const double MIN_category4l = -0.5; const double MAX_category4l = 5.5;
const int NBINS_nJets = 6; const double MIN_nJets = 1.5; const double MAX_nJets = 7.5;
const int NBINS_nBjets = 3; const double MIN_nBjets = 0.5; const double MAX_nBjets = 3.5;
const int NBINS_dilepPt = 10; const double MIN_dilepPt = 0.0; const double MAX_dilepPt = 500.0;
const int NBINS_dilepEta = 12; const double MIN_dilepEta = -3.0; const double MAX_dilepEta = 3.0;
const int NBINS_dilepPhi = 10; const double MIN_dilepPhi = -TMath::Pi(); const double MAX_dilepPhi = TMath::Pi();
const int NBINS_dilepMass3l = 10; const double MIN_dilepMass3l = particle::mZ-10.0; const double MAX_dilepMass3l = particle::mZ+10.0;
const int NBINS_dilepMass4l = 10; const double MIN_dilepMass4l = particle::mZ-20.0; const double MAX_dilepMass4l = particle::mZ+20.0;
const int NBINS_missingEt = 11; const double MIN_missingEt = 0.0; const double MAX_missingEt = 220.0;
const int NBINS_missingPhi = 10; const double MIN_missingPhi = -TMath::Pi(); const double MAX_missingPhi = TMath::Pi();
const int NBINS_firstLepPt = 10; const double MIN_firstLepPt = 40.0; const double MAX_firstLepPt = 340.0;
const int NBINS_firstLepEta = 10; const double MIN_firstLepEta = -2.5; const double MAX_firstLepEta = 2.5;
const int NBINS_firstLepPhi = 10; const double MIN_firstLepPhi = -TMath::Pi(); const double MAX_firstLepPhi = -TMath::Pi();
const int NBINS_secondLepPt = 10; const double MIN_secondLepPt = 10.0; const double MAX_secondLepPt = 160.0;
const int NBINS_secondLepEta = 10; const double MIN_secondLepEta = -2.5; const double MAX_secondLepEta = 2.5;
const int NBINS_secondLepPhi = 10; const double MIN_secondLepPhi = -TMath::Pi(); const double MAX_secondLepPhi = -TMath::Pi();
const int NBINS_thirdLepPt = 10; const double MIN_thirdLepPt = 10.0; const double MAX_thirdLepPt = 110.0;
const int NBINS_thirdLepEta = 10; const double MIN_thirdLepEta = -2.5; const double MAX_thirdLepEta = 2.5;
const int NBINS_thirdLepPhi = 10; const double MIN_thirdLepPhi = -TMath::Pi(); const double MAX_thirdLepPhi = -TMath::Pi();
const int NBINS_fourthLepPt = 10; const double MIN_fourthLepPt = 10.0; const double MAX_fourthLepPt = 110.0;
const int NBINS_fourthLepEta = 10; const double MIN_fourthLepEta = -2.5; const double MAX_fourthLepEta = 2.5;
const int NBINS_fourthLepPhi = 10; const double MIN_fourthLepPhi = -TMath::Pi(); const double MAX_fourthLepPhi = -TMath::Pi();
const int NBINS_firstJetPt = 10; const double MIN_firstJetPt = 30.0; const double MAX_firstJetPt = 330.0;
const int NBINS_firstJetEta = 10; const double MIN_firstJetEta = -2.5; const double MAX_firstJetEta = 2.5;
const int NBINS_firstJetPhi = 10; const double MIN_firstJetPhi = -TMath::Pi(); const double MAX_firstJetPhi = -TMath::Pi();
const int NBINS_secondJetPt = 10; const double MIN_secondJetPt = 30.0; const double MAX_secondJetPt = 180.0;
const int NBINS_secondJetEta = 10; const double MIN_secondJetEta = -2.5; const double MAX_secondJetEta = 2.5;
const int NBINS_secondJetPhi = 10; const double MIN_secondJetPhi = -TMath::Pi(); const double MAX_secondJetPhi = -TMath::Pi();
const int NBINS_thirdJetPt = 10; const double MIN_thirdJetPt = 30.0; const double MAX_thirdJetPt = 130.0;
const int NBINS_thirdJetEta = 10; const double MIN_thirdJetEta = -2.5; const double MAX_thirdJetEta = 2.5;
const int NBINS_thirdJetPhi = 10; const double MIN_thirdJetPhi = -TMath::Pi(); const double MAX_thirdJetPhi = -TMath::Pi();
const int NBINS_fourthJetPt = 10; const double MIN_fourthJetPt = 30.0; const double MAX_fourthJetPt = 130.0;
const int NBINS_fourthJetEta = 10; const double MIN_fourthJetEta = -2.5; const double MAX_fourthJetEta = 2.5;
const int NBINS_fourthJetPhi = 10; const double MIN_fourthJetPhi = -TMath::Pi(); const double MAX_fourthJetPhi = -TMath::Pi();
const int NBINS_lepTopMass = 10; const double MIN_lepTopMass = 170.0; const double MAX_lepTopMass = 175.0;
const int NBINS_hadTopMass = 10; const double MIN_hadTopMass = 170.0; const double MAX_hadTopMass = 175.0;
const int NBINS_zbosonPt = 8; const double BINS_zbosonPt[NBINS_zbosonPt+1] = { 0.0, 37.5, 75.0, 112.5, 150.0, 200.0, 250.0, 375.0, 500.0 };
const int NBINS_ttzMass = 8; const double BINS_ttzMass[NBINS_ttzMass+1] = { 436.0, 580.0, 640.0, 720.0, 800.0, 890.0, 1000.0, 1150.0, 1500.0 };
const int NBINS_ttbarMass = 8; const double BINS_ttbarMass[NBINS_ttbarMass+1] = { 345.0, 400.0, 440.0, 485.0, 540.0, 610.0, 700.0, 850.0, 1200.0 };
const int NBINS_topPt = 8; const double BINS_topPt[NBINS_topPt+1] = { 0.0, 60.0, 90.0, 125.0, 160.0, 200.0, 250.0, 330.0, 500.0 };
const int NBINS_deltaPhiTtbar = 8; const double BINS_deltaPhiTtbar[NBINS_deltaPhiTtbar+1] = { 0.0, TMath::Pi()/4, TMath::Pi()/2, 5*TMath::Pi()/8, 3*TMath::Pi()/4, 13*TMath::Pi()/16, 7*TMath::Pi()/8, 15*TMath::Pi()/16, TMath::Pi() };
const int NBINS_deltaPhiTopZ = 8; const double BINS_deltaPhiTopZ[NBINS_deltaPhiTopZ+1] = { 0.0, TMath::Pi()/6, TMath::Pi()/3, TMath::Pi()/2, 2*TMath::Pi()/3, 3*TMath::Pi()/4, 5*TMath::Pi()/6, 11*TMath::Pi()/12, TMath::Pi() };
const int NBINS_deltaRapTtbar = 8; const double BINS_deltaRapTtbar[NBINS_deltaRapTtbar+1] = { 0.0, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 2.1, 3.0 };
const int NBINS_deltaRapTopZ = 8; const double BINS_deltaRapTopZ[NBINS_deltaRapTopZ+1] = { 0.0, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 2.0, 3.0 };
const int NBINS_topLeptonPt = 8; const double MIN_topLeptonPt = 0.0; const double MAX_topLeptonPt = 250.0;
const int NBINS_topLeptonsMass = 8; const double MIN_topLeptonsMass = 0.0; const double MAX_topLeptonsMass = 350.0;
const int NBINS_fourLeptonsMass = 8; const double MIN_fourLeptonsMass = 100.0; const double MAX_fourLeptonsMass = 600.0;
const int NBINS_deltaPhiTopLeptons = 8; const double MIN_deltaPhiTopLeptons = 0.0; const double MAX_deltaPhiTopLeptons = TMath::Pi();
const int NBINS_deltaPhiTopLeptonZ = 8; const double MIN_deltaPhiTopLeptonZ = 0.0; const double MAX_deltaPhiTopLeptonZ = TMath::Pi();
const int NBINS_deltaRapTopLeptons = 8; const double BINS_deltaRapTopLeptons[NBINS_deltaRapTopLeptons+1] = { 0.0, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 2.0, 3.0 };
const int NBINS_deltaRapTopLeptonZ = 8; const double BINS_deltaRapTopLeptonZ[NBINS_deltaRapTopLeptonZ+1] = { 0.0, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 2.0, 3.0 };

MyHistogram::MyHistogram(std::string _name, int nVar, int nSel, int nSys, std::function<TH1F*(std::string)> hist) :
    name(_name),
    nVariants(nVar), nSelections(nSel), nSystematics(nSys),
    nHistograms(nVariants*nSelections*(1+nSystematics)),
    hists(new TH1F*[nHistograms])
{
    for(int i=0; i<nHistograms; i++) hists[i] = hist("hist_"+name+"_var"+std::to_string(i));
}

MyHistogram::~MyHistogram() {
    for(int i=0; i<nHistograms; i++) delete hists[i];
    delete [] hists;
}

void MyHistogram::Fill(int iVar, int iSel, int iSys, double weight) {
    if(iSel>=0) {
        hists[iVar+nVariants*iSel+nVariants*nSelections*(1+iSys)]->Fill(value, weight);
        if(iSys==0) hists[iVar+nVariants*iSel]->Fill(value, weight>0.0 ? 1.0 : -1.0);
    }
}

void MyHistogram::Write(TDirectory* directory, std::string dirname) {
    directory->mkdir(dirname.c_str())->cd();
    for(int iVar=0; iVar<nVariants; iVar++) {
        const std::string sVar = "hist_var"+std::to_string(iVar);
        for(int iSel=0; iSel<nSelections; iSel++) {
            const std::string sSel = sVar+"_sel"+std::to_string(iSel);
            hists[iVar+nVariants*iSel]->Write((sSel+"_stat").c_str());
            for(int iSys=0; iSys<nSystematics; iSys++) {
                const std::string sSys = sSel+"_sys"+std::to_string(iSys);
                hists[iVar+nVariants*iSel+nVariants*nSelections*(1+iSys)]->Write(sSys.c_str());
            }
        }
    }
}

MyTruthHistogram::MyTruthHistogram(std::string _name, int nVar, int nSel, int nSys, std::function<TH1F*(std::string)> hist, std::function<TH1F*(std::string)> genhist, std::function<TH2F*(std::string)> matrix) :
    MyHistogram(_name, nVar, nSel, nSys, hist),
    nGenHistograms(nVariants*(1+nSystematics)),
    gens(new TH1F*[nGenHistograms]),
    fails(new TH1F*[nGenHistograms]),
    responses(new TH2F*[nHistograms])
{
    for(int i=0; i<nGenHistograms; i++) {
        gens[i] = genhist("gen_"+name+"_var"+std::to_string(i));
        fails[i] = genhist("fail_"+name+"_var"+std::to_string(i));
    }
    for(int i=0; i<nHistograms; i++) responses[i] = matrix("response_"+name+"_var"+std::to_string(i));
}

MyTruthHistogram::~MyTruthHistogram() {
    for(int i=0; i<nGenHistograms; i++) {
        delete gens[i];
        delete fails[i];
    }
    delete [] gens;
    delete [] fails;
    for(int i=0; i<nHistograms; i++) delete responses[i];
    delete [] responses;
}

void MyTruthHistogram::Fill(int iVar, int iSel, int iSys, double weight, double genweight) {
    MyHistogram::Fill(iVar, iSel, iSys, weight);
    gens[iVar+nVariants*(1+iSys)]->Fill(genvalue, genweight);
    if(iSys==0) gens[iVar]->Fill(genvalue, genweight>0.0 ? 1.0 : -1.0);
    if(iSel>=0) {
        responses[iVar+nVariants*iSel+nVariants*nSelections*(1+iSys)]->Fill(genvalue, value, genweight);
        if(iSys==0) responses[iVar+nVariants*iSel]->Fill(genvalue, value, genweight>0.0 ? 1.0 : -1.0);
    } else {
        fails[iVar+nVariants*(1+iSys)]->Fill(genvalue, genweight);
        if(iSys==0) fails[iVar]->Fill(genvalue, genweight>0.0 ? 1.0 : -1.0);

    }
}

void MyTruthHistogram::Write(TDirectory* directory, std::string dirname, std::string gendirname) {
    MyHistogram::Write(directory, dirname);
    directory->mkdir(gendirname.c_str())->cd();
    for(int iVar=0; iVar<nVariants; iVar++) {
        const std::string sVar = "_var"+std::to_string(iVar);
        gens[iVar]->Write(("gen"+sVar+"_stat").c_str());
        for(int iSys=0; iSys<nSystematics; iSys++) {
            const std::string sSys = sVar+"_sys"+std::to_string(iSys);
            gens[iVar+nVariants*(1+iSys)]->Write(("gen"+sSys).c_str());
        }
        for(int iSel=0; iSel<nSelections; iSel++) {
            const std::string sSel = sVar+"_sel"+std::to_string(iSel);
            responses[iVar+nVariants*iSel]->Write(("response"+sSel+"_stat").c_str());
            for(int iSys=0; iSys<nSystematics; iSys++) {
                const std::string sSys = sSel+"_sys"+std::to_string(iSys);
                responses[iVar+nVariants*iSel+nVariants*nSelections*(1+iSys)]->Write(("response"+sSys).c_str());
            }
        }
        fails[iVar]->Write(("fail"+sVar+"_stat").c_str());
        for(int iSys=0; iSys<nSystematics; iSys++) {
            const std::string sSys = sVar+"_sys"+std::to_string(iSys);
            fails[iVar+nVariants*(1+iSys)]->Write(("fail"+sSys).c_str());
        }
    }
}

MyHistogram* ttZHistograms::makeHistogram(std::string name, int nVar, int nSel, int nSys, int nBins, double mini, double maxi) {
    std::function<TH1F*(std::string)> hist = [&](std::string n){ return new TH1F(n.c_str(), "", nBins, mini, maxi); };
    MyHistogram* myhist = new MyHistogram(name, nVar, nSel, nSys, hist);
    hists.push_back(myhist);
    return myhist;
}

MyHistogram* ttZHistograms::makeHistogram(std::string name, int nVar, int nSel, int nSys, int nBins, const double* bins) {
    std::function<TH1F*(std::string)> hist = [&](std::string n){ return new TH1F(n.c_str(), "", nBins, bins); };
    MyHistogram* myhist = new MyHistogram(name, nVar, nSel, nSys, hist);
    hists.push_back(myhist);
    return myhist;
}

MyHistogram* ttZHistograms::makeTruthHistogram(std::string name, bool is_ttz, int nVar, int nSel, int nSys, int nBins, double mini, double maxi) {
    if(!is_ttz) return makeHistogram(name, nVar, nSel, nSys, nBins, mini, maxi);
    std::function<TH1F*(std::string)> hist1 = [&](std::string n){ return new TH1F(n.c_str(), "", nBins, mini, maxi); };
    std::function<TH1F*(std::string)> hist2 = [&](std::string n){ return new TH1F(n.c_str(), "", nBins/2, mini, maxi); };
    std::function<TH2F*(std::string)> hist3 = [&](std::string n){ return new TH2F(n.c_str(), "", nBins/2, mini, maxi, nBins, mini, maxi); };
    MyTruthHistogram* myhist = new MyTruthHistogram(name, nVar, nSel, nSys, hist1, hist2, hist3);
    truthhists.push_back(myhist);
    return myhist;
}

MyHistogram* ttZHistograms::makeTruthHistogram(std::string name, bool is_ttz, int nVar, int nSel, int nSys, int nBins, const double* bins) {
    if(!is_ttz) return makeHistogram(name, nVar, nSel, nSys, nBins, bins);
    double genbins[nBins/2+1];
    for(int i=0; i<=nBins/2; i++) genbins[i] = bins[i*2];
    std::function<TH1F*(std::string)> hist1 = [&](std::string n){ return new TH1F(n.c_str(), "", nBins, bins); };
    std::function<TH1F*(std::string)> hist2 = [&](std::string n){ return new TH1F(n.c_str(), "", nBins/2, genbins); };
    std::function<TH2F*(std::string)> hist3 = [&](std::string n){ return new TH2F(n.c_str(), "", nBins/2, genbins, nBins, bins); };
    MyTruthHistogram* myhist = new MyTruthHistogram(name, nVar, nSel, nSys, hist1, hist2, hist3);
    truthhists.push_back(myhist);
    return myhist;
}

ttZHistograms::ttZHistograms(bool isTtz, int nVar, int nSel, int nSys) :
    is_ttz(isTtz)
{
    hists.reserve(100);
    truthhists.reserve(100);
    category3l = makeHistogram("category3l", nVar, nSel, nSys, NBINS_category3l, MIN_category3l, MAX_category3l);
    category4l = makeHistogram("category4l", nVar, nSel, nSys, NBINS_category4l, MIN_category4l, MAX_category4l);
    nJets = makeHistogram("nJets", nVar, nSel, nSys, NBINS_nJets, MIN_nJets, MAX_nJets);
    nBjets = makeHistogram("nBjets", nVar, nSel, nSys, NBINS_nBjets, MIN_nBjets, MAX_nBjets);
    dilepPt = makeHistogram("dilepPt", nVar, nSel, nSys, NBINS_dilepPt, MIN_dilepPt, MAX_dilepPt);
    dilepEta = makeHistogram("dilepEta", nVar, nSel, nSys, NBINS_dilepEta, MIN_dilepEta, MAX_dilepEta);
    dilepPhi = makeHistogram("dilepPhi", nVar, nSel, nSys, NBINS_dilepPhi, MIN_dilepPhi, MAX_dilepPhi);
    dilepMass3l = makeHistogram("dilepMass3l", nVar, nSel, nSys, NBINS_dilepMass3l, MIN_dilepMass3l, MAX_dilepMass3l);
    dilepMass4l = makeHistogram("dilepMass4l", nVar, nSel, nSys, NBINS_dilepMass4l, MIN_dilepMass4l, MAX_dilepMass4l);
    missingEt = makeHistogram("missingEt", nVar, nSel, nSys, NBINS_missingEt, MIN_missingEt, MAX_missingEt);
    missingPhi = makeHistogram("missingPhi", nVar, nSel, nSys, NBINS_missingPhi, MIN_missingPhi, MAX_missingPhi);
    firstLepPt = makeHistogram("firstLepPt", nVar, nSel, nSys, NBINS_firstLepPt, MIN_firstLepPt, MAX_firstLepPt);
    firstLepEta = makeHistogram("firstLepEta", nVar, nSel, nSys, NBINS_firstLepEta, MIN_firstLepEta, MAX_firstLepEta);
    firstLepPhi = makeHistogram("firstLepPhi", nVar, nSel, nSys, NBINS_firstLepPhi, MIN_firstLepPhi, MAX_firstLepPhi);
    secondLepPt = makeHistogram("secondLepPt", nVar, nSel, nSys, NBINS_secondLepPt, MIN_secondLepPt, MAX_secondLepPt);
    secondLepEta = makeHistogram("secondLepEta", nVar, nSel, nSys, NBINS_secondLepEta, MIN_secondLepEta, MAX_secondLepEta);
    secondLepPhi = makeHistogram("secondLepPhi", nVar, nSel, nSys, NBINS_secondLepPhi, MIN_secondLepPhi, MAX_secondLepPhi);
    thirdLepPt = makeHistogram("thirdLepPt", nVar, nSel, nSys, NBINS_thirdLepPt, MIN_thirdLepPt, MAX_thirdLepPt);
    thirdLepEta = makeHistogram("thirdLepEta", nVar, nSel, nSys, NBINS_thirdLepEta, MIN_thirdLepEta, MAX_thirdLepEta);
    thirdLepPhi = makeHistogram("thirdLepPhi", nVar, nSel, nSys, NBINS_thirdLepPhi, MIN_thirdLepPhi, MAX_thirdLepPhi);
    fourthLepPt = makeHistogram("fourthLepPt", nVar, nSel, nSys, NBINS_fourthLepPt, MIN_fourthLepPt, MAX_fourthLepPt);
    fourthLepEta = makeHistogram("fourthLepEta", nVar, nSel, nSys, NBINS_fourthLepEta, MIN_fourthLepEta, MAX_fourthLepEta);
    fourthLepPhi = makeHistogram("fourthLepPhi", nVar, nSel, nSys, NBINS_fourthLepPhi, MIN_fourthLepPhi, MAX_fourthLepPhi);
    firstJetPt = makeHistogram("firstJetPt", nVar, nSel, nSys, NBINS_firstJetPt, MIN_firstJetPt, MAX_firstJetPt);
    firstJetEta = makeHistogram("firstJetEta", nVar, nSel, nSys, NBINS_firstJetEta, MIN_firstJetEta, MAX_firstJetEta);
    firstJetPhi = makeHistogram("firstJetPhi", nVar, nSel, nSys, NBINS_firstJetPhi, MIN_firstJetPhi, MAX_firstJetPhi);
    secondJetPt = makeHistogram("secondJetPt", nVar, nSel, nSys, NBINS_secondJetPt, MIN_secondJetPt, MAX_secondJetPt);
    secondJetEta = makeHistogram("secondJetEta", nVar, nSel, nSys, NBINS_secondJetEta, MIN_secondJetEta, MAX_secondJetEta);
    secondJetPhi = makeHistogram("secondJetPhi", nVar, nSel, nSys, NBINS_secondJetPhi, MIN_secondJetPhi, MAX_secondJetPhi);
    thirdJetPt = makeHistogram("thirdJetPt", nVar, nSel, nSys, NBINS_thirdJetPt, MIN_thirdJetPt, MAX_thirdJetPt);
    thirdJetEta = makeHistogram("thirdJetEta", nVar, nSel, nSys, NBINS_thirdJetEta, MIN_thirdJetEta, MAX_thirdJetEta);
    thirdJetPhi = makeHistogram("thirdJetPhi", nVar, nSel, nSys, NBINS_thirdJetPhi, MIN_thirdJetPhi, MAX_thirdJetPhi);
    fourthJetPt = makeHistogram("fourthJetPt", nVar, nSel, nSys, NBINS_fourthJetPt, MIN_fourthJetPt, MAX_fourthJetPt);
    fourthJetEta = makeHistogram("fourthJetEta", nVar, nSel, nSys, NBINS_fourthJetEta, MIN_fourthJetEta, MAX_fourthJetEta);
    fourthJetPhi = makeHistogram("fourthJetPhi", nVar, nSel, nSys, NBINS_fourthJetPhi, MIN_fourthJetPhi, MAX_fourthJetPhi);
    lepTopMass = makeHistogram("lepTopMass", nVar, nSel, nSys, NBINS_lepTopMass, MIN_lepTopMass, MAX_lepTopMass);
    hadTopMass = makeHistogram("hadTopMass", nVar, nSel, nSys, NBINS_hadTopMass, MIN_hadTopMass, MAX_hadTopMass);
    zbosonPt = makeTruthHistogram("zbosonPt", is_ttz, nVar, nSel, nSys, NBINS_zbosonPt, BINS_zbosonPt);
    ttzMass = makeTruthHistogram("ttzMass", is_ttz, nVar, nSel, nSys, NBINS_ttzMass, BINS_ttzMass);
    ttbarMass = makeTruthHistogram("ttbarMass", is_ttz, nVar, nSel, nSys, NBINS_ttbarMass, BINS_ttbarMass);
    topPt = makeTruthHistogram("topPt", is_ttz, nVar, nSel, nSys, NBINS_topPt, BINS_topPt);
    deltaPhiTtbar = makeTruthHistogram("deltaPhiTtbar", is_ttz, nVar, nSel, nSys, NBINS_deltaPhiTtbar, BINS_deltaPhiTtbar);
    deltaPhiTopZ = makeTruthHistogram("deltaPhiTopZ", is_ttz, nVar, nSel, nSys, NBINS_deltaPhiTopZ, BINS_deltaPhiTopZ);
    deltaRapTtbar = makeTruthHistogram("deltaRapTtbar", is_ttz, nVar, nSel, nSys, NBINS_deltaRapTtbar, BINS_deltaRapTtbar);
    deltaRapTopZ = makeTruthHistogram("deltaRapTopZ", is_ttz, nVar, nSel, nSys, NBINS_deltaRapTopZ, BINS_deltaRapTopZ);
    topLeptonPt = makeTruthHistogram("topLeptonPt", is_ttz, nVar, nSel, nSys, NBINS_topLeptonPt, MIN_topLeptonPt, MAX_topLeptonPt);
    topLeptonsMass = makeTruthHistogram("topLeptonsMass", is_ttz, nVar, nSel, nSys, NBINS_topLeptonsMass, MIN_topLeptonsMass, MAX_topLeptonsMass);
    fourLeptonsMass = makeTruthHistogram("fourLeptonsMass", is_ttz, nVar, nSel, nSys, NBINS_fourLeptonsMass, MIN_fourLeptonsMass, MAX_fourLeptonsMass);
    deltaPhiTopLeptons = makeTruthHistogram("deltaPhiTopLeptons", is_ttz, nVar, nSel, nSys, NBINS_deltaPhiTopLeptons, MIN_deltaPhiTopLeptons, MAX_deltaPhiTopLeptons);
    deltaPhiTopLeptonZ = makeTruthHistogram("deltaPhiTopLeptonZ", is_ttz, nVar, nSel, nSys, NBINS_deltaPhiTopLeptonZ, MIN_deltaPhiTopLeptonZ, MAX_deltaPhiTopLeptonZ);
    deltaRapTopLeptons = makeTruthHistogram("deltaRapTopLeptons", is_ttz, nVar, nSel, nSys, NBINS_deltaRapTopLeptons, BINS_deltaRapTopLeptons);
    deltaRapTopLeptonZ = makeTruthHistogram("deltaRapTopLeptonZ", is_ttz, nVar, nSel, nSys, NBINS_deltaRapTopLeptonZ, BINS_deltaRapTopLeptonZ);
    hists.shrink_to_fit();
    truthhists.shrink_to_fit();
}

ttZHistograms::~ttZHistograms() {
    for(std::vector<MyHistogram*>::iterator hist=hists.begin(); hist!=hists.end(); ++hist) delete (*hist);
    for(std::vector<MyTruthHistogram*>::iterator hist=truthhists.begin(); hist!=truthhists.end(); ++hist) delete (*hist);
}

#define IN_RANGE(VAL, MIN, MAX) \
    (VAL<-900.0 ? -999.0 : TMath::Min(MAX-0.001*(MAX-MIN), TMath::Max(MIN+0.001*(MAX-MIN), VAL)))
#define PHI_IN_RANGE(VAL) \
    (VAL<-900.0 || isnan(VAL) ? -999.0 : TMath::Min(0.999*TMath::Pi(), TMath::Max(-0.999*TMath::Pi(), TVector2::Phi_mpi_pi(VAL))))

void ttZHistograms::SetValues(ttZ::leptonVariables leptonVars, ttZ::jetVariables jetVars, ttZ::reconstructedVariables reconstructedVars, ttZ::fourLeptonVariables fourLeptonVars) {
    category3l->SetValue(leptonVars.isThreeLeptons ? leptonVars.category : -1);
    category4l->SetValue(leptonVars.isThreeLeptons ? -1 : leptonVars.category);
    nJets->SetValue(jetVars.nJets>MAX_nJets ? MAX_nJets-0.5 : jetVars.nJets);
    nBjets->SetValue(jetVars.nBjets>MAX_nBjets ? MAX_nBjets-0.5 : jetVars.nBjets);
    dilepPt->SetValue(IN_RANGE(leptonVars.dilepPt, MIN_dilepPt, MAX_dilepPt));
    dilepEta->SetValue(IN_RANGE(leptonVars.dilepEta, MIN_dilepEta, MAX_dilepEta));
    dilepPhi->SetValue(PHI_IN_RANGE(leptonVars.dilepPhi));
    dilepMass3l->SetValue(IN_RANGE(leptonVars.dilepMass, MIN_dilepMass3l, MAX_dilepMass3l));
    dilepMass4l->SetValue(IN_RANGE(leptonVars.dilepMass, MIN_dilepMass4l, MAX_dilepMass4l));
    missingEt->SetValue(IN_RANGE(jetVars.missingEt, MIN_missingEt, MAX_missingEt));
    missingPhi->SetValue(PHI_IN_RANGE(jetVars.missingPhi));
    firstLepPt->SetValue(IN_RANGE(leptonVars.firstLepPt, MIN_firstLepPt, MAX_firstLepPt));
    firstLepEta->SetValue(IN_RANGE(leptonVars.firstLepEta, MIN_firstLepEta, MAX_firstLepEta));
    firstLepPhi->SetValue(PHI_IN_RANGE(leptonVars.firstLepPhi));
    secondLepPt->SetValue(IN_RANGE(leptonVars.secondLepPt, MIN_secondLepPt, MAX_secondLepPt));
    secondLepEta->SetValue(IN_RANGE(leptonVars.secondLepEta, MIN_secondLepEta, MAX_secondLepEta));
    secondLepPhi->SetValue(PHI_IN_RANGE(leptonVars.secondLepPhi));
    thirdLepPt->SetValue(IN_RANGE(leptonVars.thirdLepPt, MIN_thirdLepPt, MAX_thirdLepPt));
    thirdLepEta->SetValue(IN_RANGE(leptonVars.thirdLepEta, MIN_thirdLepEta, MAX_thirdLepEta));
    thirdLepPhi->SetValue(PHI_IN_RANGE(leptonVars.thirdLepPhi));
    fourthLepPt->SetValue(IN_RANGE(leptonVars.fourthLepPt, MIN_fourthLepPt, MAX_fourthLepPt));
    fourthLepEta->SetValue(IN_RANGE(leptonVars.fourthLepEta, MIN_fourthLepEta, MAX_fourthLepEta));
    fourthLepPhi->SetValue(PHI_IN_RANGE(leptonVars.thirdLepPhi));
    firstJetPt->SetValue(IN_RANGE(jetVars.firstJetPt, MIN_firstJetPt, MAX_firstJetPt));
    firstJetEta->SetValue(IN_RANGE(jetVars.firstJetEta, MIN_firstJetEta, MAX_firstJetEta));
    firstJetPhi->SetValue(PHI_IN_RANGE(jetVars.firstJetPhi));
    secondJetPt->SetValue(IN_RANGE(jetVars.secondJetPt, MIN_secondJetPt, MAX_secondJetPt));
    secondJetEta->SetValue(IN_RANGE(jetVars.secondJetEta, MIN_secondJetEta, MAX_secondJetEta));
    secondJetPhi->SetValue(PHI_IN_RANGE(jetVars.secondJetPhi));
    thirdJetPt->SetValue(IN_RANGE(jetVars.thirdJetPt, MIN_thirdJetPt, MAX_thirdJetPt));
    thirdJetEta->SetValue(IN_RANGE(jetVars.thirdJetEta, MIN_thirdJetEta, MAX_thirdJetEta));
    thirdJetPhi->SetValue(PHI_IN_RANGE(jetVars.thirdJetPhi));
    fourthJetPt->SetValue(IN_RANGE(jetVars.fourthJetPt, MIN_fourthJetPt, MAX_fourthJetPt));
    fourthJetEta->SetValue(IN_RANGE(jetVars.fourthJetEta, MIN_fourthJetEta, MAX_fourthJetEta));
    fourthJetPhi->SetValue(PHI_IN_RANGE(jetVars.fourthJetPhi));
    lepTopMass->SetValue(IN_RANGE(reconstructedVars.lepTopMass, MIN_lepTopMass, MAX_lepTopMass));
    hadTopMass->SetValue(IN_RANGE(reconstructedVars.hadTopMass, MIN_hadTopMass, MAX_hadTopMass));
    zbosonPt->SetValue(IN_RANGE(leptonVars.dilepPt, BINS_zbosonPt[0], BINS_zbosonPt[NBINS_zbosonPt]));
    ttzMass->SetValue(IN_RANGE(reconstructedVars.ttzMass, BINS_ttzMass[0], BINS_ttzMass[NBINS_ttzMass]));
    ttbarMass->SetValue(IN_RANGE(reconstructedVars.ttbarMass, BINS_ttbarMass[0], BINS_ttbarMass[NBINS_ttbarMass]));
    topPt->SetValue(IN_RANGE(reconstructedVars.topPt, BINS_topPt[0], BINS_topPt[NBINS_topPt]));
    deltaPhiTtbar->SetValue(IN_RANGE(reconstructedVars.deltaPhiTtbar, BINS_deltaPhiTtbar[0], BINS_deltaPhiTtbar[NBINS_deltaPhiTtbar]));
    deltaPhiTopZ->SetValue(IN_RANGE(reconstructedVars.deltaPhiTopZ, BINS_deltaPhiTopZ[0], BINS_deltaPhiTopZ[NBINS_deltaPhiTopZ]));
    deltaRapTtbar->SetValue(IN_RANGE(reconstructedVars.deltaRapTtbar, BINS_deltaRapTtbar[0], BINS_deltaRapTtbar[NBINS_deltaRapTtbar]));
    deltaRapTopZ->SetValue(IN_RANGE(reconstructedVars.deltaRapTopZ, BINS_deltaRapTopZ[0], BINS_deltaRapTopZ[NBINS_deltaRapTopZ]));
    topLeptonPt->SetValue(IN_RANGE(fourLeptonVars.topLeptonPt, MIN_topLeptonPt, MAX_topLeptonPt));
    topLeptonsMass->SetValue(IN_RANGE(fourLeptonVars.topLeptonsMass, MIN_topLeptonsMass, MAX_topLeptonsMass));
    fourLeptonsMass->SetValue(IN_RANGE(fourLeptonVars.fourLeptonsMass, MIN_fourLeptonsMass, MAX_fourLeptonsMass));
    deltaPhiTopLeptons->SetValue(IN_RANGE(fourLeptonVars.deltaPhiTopLeptons, MIN_deltaPhiTopLeptons, MAX_deltaPhiTopLeptons));
    deltaPhiTopLeptonZ->SetValue(IN_RANGE(fourLeptonVars.deltaPhiTopLeptonZ, MIN_deltaPhiTopLeptonZ, MAX_deltaPhiTopLeptonZ));
    deltaRapTopLeptons->SetValue(IN_RANGE(fourLeptonVars.deltaRapTopLeptons, BINS_deltaRapTopLeptons[0], BINS_deltaRapTopLeptons[NBINS_deltaRapTopLeptons]));
    deltaRapTopLeptonZ->SetValue(IN_RANGE(fourLeptonVars.deltaRapTopLeptonZ, BINS_deltaRapTopLeptonZ[0], BINS_deltaRapTopLeptonZ[NBINS_deltaRapTopLeptonZ]));
}

void ttZHistograms::SetTruthValues(ttZ::ttzTruth truthVars) {
    ((MyTruthHistogram*)zbosonPt)->SetTruthValue(IN_RANGE(truthVars.zbosonPt, BINS_zbosonPt[0], BINS_zbosonPt[NBINS_zbosonPt]));
    ((MyTruthHistogram*)ttzMass)->SetTruthValue(IN_RANGE(truthVars.ttzMass, BINS_ttzMass[0], BINS_ttzMass[NBINS_ttzMass]));
    ((MyTruthHistogram*)ttbarMass)->SetTruthValue(IN_RANGE(truthVars.ttbarMass, BINS_ttbarMass[0], BINS_ttbarMass[NBINS_ttbarMass]));
    ((MyTruthHistogram*)topPt)->SetTruthValue(IN_RANGE(truthVars.topPt, BINS_topPt[0], BINS_topPt[NBINS_topPt]));
    ((MyTruthHistogram*)deltaPhiTtbar)->SetTruthValue(IN_RANGE(truthVars.deltaPhiTtbar, BINS_deltaPhiTtbar[0], BINS_deltaPhiTtbar[NBINS_deltaPhiTtbar]));
    ((MyTruthHistogram*)deltaPhiTopZ)->SetTruthValue(IN_RANGE(truthVars.deltaPhiTopZ, BINS_deltaPhiTopZ[0], BINS_deltaPhiTopZ[NBINS_deltaPhiTopZ]));
    ((MyTruthHistogram*)deltaRapTtbar)->SetTruthValue(IN_RANGE(truthVars.deltaRapTtbar, BINS_deltaRapTtbar[0], BINS_deltaRapTtbar[NBINS_deltaRapTtbar]));
    ((MyTruthHistogram*)deltaRapTopZ)->SetTruthValue(IN_RANGE(truthVars.deltaRapTopZ, BINS_deltaRapTopZ[0], BINS_deltaRapTopZ[NBINS_deltaRapTopZ]));
    ((MyTruthHistogram*)topLeptonPt)->SetTruthValue(IN_RANGE(truthVars.topLeptonPt, MIN_topLeptonPt, MAX_topLeptonPt));
    ((MyTruthHistogram*)topLeptonsMass)->SetTruthValue(IN_RANGE(truthVars.topLeptonsMass, MIN_topLeptonsMass, MAX_topLeptonsMass));
    ((MyTruthHistogram*)fourLeptonsMass)->SetTruthValue(IN_RANGE(truthVars.fourLeptonsMass, MIN_fourLeptonsMass, MAX_fourLeptonsMass));
    ((MyTruthHistogram*)deltaPhiTopLeptons)->SetTruthValue(IN_RANGE(truthVars.deltaPhiTopLeptons, MIN_deltaPhiTopLeptons, MAX_deltaPhiTopLeptons));
    ((MyTruthHistogram*)deltaPhiTopLeptonZ)->SetTruthValue(IN_RANGE(truthVars.deltaPhiTopLeptonZ, MIN_deltaPhiTopLeptonZ, MAX_deltaPhiTopLeptonZ));
    ((MyTruthHistogram*)deltaRapTopLeptons)->SetTruthValue(IN_RANGE(truthVars.deltaRapTopLeptons, BINS_deltaRapTopLeptons[0], BINS_deltaRapTopLeptons[NBINS_deltaRapTopLeptons]));
    ((MyTruthHistogram*)deltaRapTopLeptonZ)->SetTruthValue(IN_RANGE(truthVars.deltaRapTopLeptonZ, BINS_deltaRapTopLeptonZ[0], BINS_deltaRapTopLeptonZ[NBINS_deltaRapTopLeptonZ]));
}

void ttZHistograms::Fill(int iVar, int iSel, int iSys, double weight, double genweight) {
    for(std::vector<MyHistogram*>::iterator hist=hists.begin(); hist!=hists.end(); ++hist) (*hist)->Fill(iVar, iSel, iSys, weight);
    for(std::vector<MyTruthHistogram*>::iterator hist=truthhists.begin(); hist!=truthhists.end(); ++hist) (*hist)->Fill(iVar, iSel, iSys, weight, genweight);
}

void ttZHistograms::Write(TDirectory* directory) {
    category3l->Write(directory, "control_category3l");
    category4l->Write(directory, "control_category4l");
    nJets->Write(directory, "control_nJets");
    nBjets->Write(directory, "control_nBjets");
    dilepPt->Write(directory, "control_dilepPt");
    dilepEta->Write(directory, "control_dilepEta");
    dilepPhi->Write(directory, "control_dilepPhi");
    dilepMass3l->Write(directory, "control_dilepMass3l");
    dilepMass4l->Write(directory, "control_dilepMass4l");
    missingEt->Write(directory, "control_missingEt");
    missingPhi->Write(directory, "control_missingPhi");
    firstLepPt->Write(directory, "control_firstLepPt");
    firstLepEta->Write(directory, "control_firstLepEta");
    firstLepPhi->Write(directory, "control_firstLepPhi");
    secondLepPt->Write(directory, "control_secondLepPt");
    secondLepEta->Write(directory, "control_secondLepEta");
    secondLepPhi->Write(directory, "control_secondLepPhi");
    thirdLepPt->Write(directory, "control_thirdLepPt");
    thirdLepEta->Write(directory, "control_thirdLepEta");
    thirdLepPhi->Write(directory, "control_thirdLepPhi");
    fourthLepPt->Write(directory, "control_fourthLepPt");
    fourthLepEta->Write(directory, "control_fourthLepEta");
    fourthLepPhi->Write(directory, "control_fourthLepPhi");
    firstJetPt->Write(directory, "control_firstJetPt");
    firstJetEta->Write(directory, "control_firstJetEta");
    firstJetPhi->Write(directory, "control_firstJetPhi");
    secondJetPt->Write(directory, "control_secondJetPt");
    secondJetEta->Write(directory, "control_secondJetEta");
    secondJetPhi->Write(directory, "control_secondJetPhi");
    thirdJetPt->Write(directory, "control_thirdJetPt");
    thirdJetEta->Write(directory, "control_thirdJetEta");
    thirdJetPhi->Write(directory, "control_thirdJetPhi");
    fourthJetPt->Write(directory, "control_fourthJetPt");
    fourthJetEta->Write(directory, "control_fourthJetEta");
    fourthJetPhi->Write(directory, "control_fourthJetPhi");
    lepTopMass->Write(directory, "reco3l_lepTopMass");
    hadTopMass->Write(directory, "reco3l_hadTopMass");
    zbosonPt->Write(directory, "reco_zbosonPt", "response_zbosonPt");
    ttzMass->Write(directory, "reco3l_ttzMass", "response3l_ttzMass");
    ttbarMass->Write(directory, "reco3l_ttbarMass", "response3l_ttbarMass");
    topPt->Write(directory, "reco3l_topPt", "response3l_topPt");
    deltaPhiTtbar->Write(directory, "reco3l_deltaPhiTtbar", "response3l_deltaPhiTtbar");
    deltaPhiTopZ->Write(directory, "reco3l_deltaPhiTopZ", "response3l_deltaPhiTopZ");
    deltaRapTtbar->Write(directory, "reco3l_deltaRapTtbar", "response3l_deltaRapTtbar");
    deltaRapTopZ->Write(directory, "reco3l_deltaRapTopZ", "response3l_deltaRapTopZ");
    topLeptonPt->Write(directory, "reco4l_topLeptonPt", "response4l_topLeptonPt");
    topLeptonsMass->Write(directory, "reco4l_topLeptonsMass", "response4l_topLeptonsMass");
    fourLeptonsMass->Write(directory, "reco4l_fourLeptonsMass", "response4l_fourLeptonsMass");
    deltaPhiTopLeptons->Write(directory, "reco4l_deltaPhiTopLeptons", "response4l_deltaPhiTopLeptons");
    deltaPhiTopLeptonZ->Write(directory, "reco4l_deltaPhiTopLeptonZ", "response4l_deltaPhiTopLeptonZ");
    deltaRapTopLeptons->Write(directory, "reco4l_deltaRapTopLeptons", "response4l_deltaRapTopLeptons");
    deltaRapTopLeptonZ->Write(directory, "reco4l_deltaRapTopLeptonZ", "response4l_deltaRapTopLeptonZ");
}
