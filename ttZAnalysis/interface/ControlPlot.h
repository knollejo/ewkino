#ifndef ControlPlot_H
#define ControlPlot_H

#include "HistogramCombiner.h"

class ControlPlot : public HistogramCombiner {
public:
    ControlPlot(int, int, int);
    ~ControlPlot();

    void Evaluate();
    std::string Print(int, int);

    static const int nYears = 3;
    static const int nUncertainties = 7;
    double**** groups;
    double**** uncSqUp, **** uncSqDown;
};

#endif
