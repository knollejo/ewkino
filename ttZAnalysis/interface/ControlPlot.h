#ifndef ControlPlot_H
#define ControlPlot_H

#include "HistogramCombiner.h"

class ControlPlot : public HistogramCombiner {
public:
    ControlPlot(int, int, int);
    ~ControlPlot();

    void Evaluate();
    std::string Print(int, int);
    void Dump(int, int);

    static const int nYears = 3;
    double**** groups;
    double**** uncSq;
};

#endif
