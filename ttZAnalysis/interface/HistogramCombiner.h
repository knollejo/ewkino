#ifndef HistogramCombiner_H
#define HistogramCombiner_H

#include <string>
#include <tuple>
#include <vector>

class HistogramCombiner {
public:
    HistogramCombiner(int nSel, int nGrp, int nBin) : nSelections(nSel), nGroups(nGrp), nBins(nBin) {}
    virtual ~HistogramCombiner() {}

    void SetBins(double*, double*);
    void AddSample(double***, int, int, bool, double=1.0);

    void EvaluateGroup(int, int, double***);
    void EvaluateEnvelopeUncorrelated(int, int, std::vector<int>, double***, double***);
    void EvaluateEnvelopeCorrelated(int, std::vector<int>, double***, double***);

    void EvaluateCentral(int, double***);
    void EvaluateStatistical(int, double***);
    void EvaluateJetMet(int, double***, double***);
    void EvaluateBtagging(int, double***, double***);
    void EvaluateLeptons(int, double***, double***);
    void EvaluatePileup(int, double***, double***);
    void EvaluateMeScales(int, double***, double***);

    std::string PrintSingle(double*, std::vector<double*>, double*, double*);
    std::string PrintDouble(double*, std::vector<double*>, double*, double*);

    static const int nYears = 3;
    const int nSelections, nGroups, nBins;
    double* binCenters = nullptr, * binBoundaries = nullptr;

protected:
    typedef std::tuple<double***, int, int, bool, double> Sample;
    std::vector<Sample> samples;
};

#endif
