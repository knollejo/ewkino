#ifndef GenHistogramCombiner_H
#define GenHistogramCombiner_H

#include <string>
#include <tuple>
#include <vector>

#include "TH2.h"

class GenHistogramCombiner {
public:
    GenHistogramCombiner(int, int, int, int);
    virtual ~GenHistogramCombiner();

    void SetBins(double**, double**);
    void AddSample(double**, double**, double***, int, int, double=1.0);

    void EvaluateGenOrFail(int, int, double**, bool);
    void EvaluateGen(int iDef, int iSys, double** vals) { EvaluateGenOrFail(iDef, iSys, vals, true); }
    void EvaluateFail(int iDef, int iSys, double** vals) { EvaluateGenOrFail(iDef, iSys, vals, false); }
    void EvaluateResponse(int, int, double***);

    void Evaluate();

    std::string PrintResponse(int, int);
    std::string PrintStabPurEff(int, int);

    static const int nYears = 3;
    const int nDefinitions, nSystematics, nBinsX, nBinsY;
    double** binCenters = nullptr, ** binBoundaries = nullptr;

protected:
    typedef std::tuple<double**, double**, double***, int, int, double> Sample;
    std::vector<Sample> samples;

    TH2F**** responseMatrices;
    double*** genDistributions, **** genUncertainties;
};

#endif
