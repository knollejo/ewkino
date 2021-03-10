#ifndef ttZHistograms_H
#define ttZHistograms_H

//include c++ library classes
#include <functional>
#include <string>

//include ROOT classes
#include "TDirectory.h"
#include "TMath.h"
#include "TH1.h"
#include "TH2.h"

//include ttZ specific code
#include "ttZTruth.h"
#include "ttZVariables.h"

class MyHistogram {
public:
    MyHistogram(std::string, int, int, int, std::function<TH1F*(std::string)>);
    virtual ~MyHistogram();

    void SetValue(double v) { value = v; }
    virtual void Fill(int, int, int, double);
    virtual void Write(TDirectory*, std::string);
    virtual void Write(TDirectory* directory, std::string dirname, std::string) { Write(directory, dirname); }

protected:
    std::string name;
    const int nVariants, nSelections, nSystematics;
    const int nHistograms;
    TH1F** hists;
    double value;
};

class MyTruthHistogram : public MyHistogram {
public:
    MyTruthHistogram(std::string, int, int, int, std::function<TH1F*(std::string)>, std::function<TH1F*(std::string)>, std::function<TH2F*(std::string)>);
    virtual ~MyTruthHistogram();

    void SetTruthValue(double g) { genvalue = g; }
    virtual void Fill(int, int, int, double, double);
    virtual void Write(TDirectory*, std::string, std::string);

protected:
    const int nGenHistograms;
    TH1F** gens, ** fails;
    TH2F** responses;
    double genvalue;
};


class ttZHistograms {
protected:
    MyHistogram* makeHistogram(std::string, int, int, int, int, double, double);
    MyHistogram* makeHistogram(std::string, int, int, int, int, const double*);
    MyHistogram* makeTruthHistogram(std::string, bool, int, int, int, int, double, double);
    MyHistogram* makeTruthHistogram(std::string, bool, int, int, int, int, const double*);

public:
    ttZHistograms(bool, int, int, int);
    virtual ~ttZHistograms();

    void SetValues(ttZ::leptonVariables, ttZ::jetVariables, ttZ::reconstructedVariables, ttZ::fourLeptonVariables);
    void SetTruthValues(ttZ::ttzTruth);
    void Fill(int, int, int, double, double);

    void Write(TDirectory*);

protected:
    const bool is_ttz;
    std::vector<MyHistogram*> hists;
    std::vector<MyTruthHistogram*> truthhists;
    MyHistogram* category3l=nullptr, * category4l=nullptr,
               * nJets=nullptr, * nBjets=nullptr,
               * dilepPt=nullptr, * dilepEta=nullptr, * dilepPhi=nullptr, * dilepMass3l=nullptr, * dilepMass4l=nullptr,
               * missingEt=nullptr, * missingPhi=nullptr,
               * firstLepPt=nullptr, * firstLepEta=nullptr, * firstLepPhi=nullptr,
               * secondLepPt=nullptr, * secondLepEta=nullptr, * secondLepPhi=nullptr,
               * thirdLepPt=nullptr, * thirdLepEta=nullptr, * thirdLepPhi=nullptr,
               * fourthLepPt=nullptr, * fourthLepEta=nullptr, * fourthLepPhi=nullptr,
               * firstJetPt=nullptr, * firstJetEta=nullptr, * firstJetPhi=nullptr,
               * secondJetPt=nullptr, * secondJetEta=nullptr, * secondJetPhi=nullptr,
               * thirdJetPt=nullptr, * thirdJetEta=nullptr, * thirdJetPhi=nullptr,
               * fourthJetPt=nullptr, * fourthJetEta=nullptr, * fourthJetPhi=nullptr,
               * lepTopMass=nullptr, * hadTopMass=nullptr,
               * zbosonPt=nullptr, * ttzMass=nullptr, * ttbarMass=nullptr, * topPt=nullptr,
               * deltaPhiTtbar=nullptr, * deltaPhiTopZ=nullptr, * deltaRapTtbar=nullptr, * deltaRapTopZ=nullptr,
               * topLeptonPt=nullptr, * topLeptonsMass=nullptr, * fourLeptonsMass=nullptr,
               * deltaPhiTopLeptons=nullptr, * deltaPhiTopLeptonZ=nullptr, * deltaRapTopLeptons=nullptr, * deltaRapTopLeptonZ=nullptr;
};

#endif
