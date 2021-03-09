#include <cmath>
#include <iostream>
#include <sstream>

#include "../interface/ControlPlot.h"

ControlPlot::ControlPlot(int nSel, int nGrp, int nBin) :
    HistogramCombiner(nSel, nGrp, nBin)
{
    groups = new double***[nGroups];
    for(int iGrp=0; iGrp<nGroups; iGrp++) {
        groups[iGrp] = new double**[nYears+1];
        for(int iYear=0; iYear<=nYears; iYear++) {
            groups[iGrp][iYear] = new double*[nSelections];
            for(int iSel=0; iSel<nSelections; iSel++) {
                groups[iGrp][iYear][iSel] = new double[nBins];
                for(int iBin=0; iBin<nBins; iBin++) {
                    groups[iGrp][iYear][iSel][iBin] = 0.0;
                }
            }
        }
    }
    uncSqUp = new double***[nUncertainties];
    uncSqDown = new double***[nUncertainties];
    for(int iUnc=0; iUnc<nUncertainties; iUnc++) {
        uncSqUp[iUnc] = new double**[nYears+1];
        uncSqDown[iUnc] = new double**[nYears+1];
        for(int iYear=0; iYear<=nYears; iYear++) {
            uncSqUp[iUnc][iYear] = new double*[nSelections];
            uncSqDown[iUnc][iYear] = new double*[nSelections];
            for(int iSel=0; iSel<nSelections; iSel++) {
                uncSqUp[iUnc][iYear][iSel] = new double[nBins];
                uncSqDown[iUnc][iYear][iSel] = new double[nBins];
                for(int iBin=0; iBin<nBins; iBin++) {
                    uncSqUp[iUnc][iYear][iSel][iBin] = 0.0;
                    uncSqDown[iUnc][iYear][iSel][iBin] = 0.0;
                }
            }
        }
    }
}

ControlPlot::~ControlPlot() {
    for(int iGrp=0; iGrp<nGroups; iGrp++) {
        for(int iYear=0; iYear<=nYears; iYear++) {
            for(int iSel=0; iSel<nSelections; iSel++) {
                delete [] groups[iGrp][iYear][iSel];
            }
            delete [] groups[iGrp][iYear];
        }
        delete [] groups[iGrp];
    }
    delete [] groups;
    for(int iUnc=0; iUnc<nUncertainties; iUnc++) {
        for(int iYear=0; iYear<=nYears; iYear++) {
            for(int iSel=0; iSel<nSelections; iSel++) {
                delete [] uncSqUp[iUnc][iYear][iSel];
                delete [] uncSqDown[iUnc][iYear][iSel];
            }
            delete [] uncSqUp[iUnc][iYear];
            delete [] uncSqDown[iUnc][iYear];
        }
        delete [] uncSqUp[iUnc];
        delete [] uncSqDown[iUnc];
    }
    delete [] uncSqUp;
    delete [] uncSqDown;
}

void ControlPlot::Evaluate() {
    std::cout << "Evaluating" << std::endl;
    // central values
    for(int iGroup=0; iGroup<nGroups; iGroup++) {
        EvaluateCentral(iGroup, groups[iGroup]);
    }
    // MC statistics (uncertainty 0)
    for(int iGroup=1; iGroup<nGroups-1; iGroup++) {
        EvaluateStatistical(iGroup, uncSqUp[0]);
        EvaluateStatistical(iGroup, uncSqDown[0]);
    }
    // nonprompt background statistics (uncertainty 1)
    EvaluateStatistical(nGroups-1, uncSqUp[1]);
    EvaluateStatistical(nGroups-1, uncSqDown[1]);
    // ME scale uncertainty (uncertainty 2)
    EvaluateMeScales(1, uncSqUp[2], uncSqDown[2]);
    // JES, JER, and MET variations (uncertainty 3)
    EvaluateJetMet(1, uncSqUp[3], uncSqDown[3]);
    // Btagging scale factor uncertainties (uncertainty 4)
    EvaluateBtagging(1, uncSqUp[4], uncSqDown[4]);
    // Pileup reweighting uncertainty (uncertainty 5)
    EvaluatePileup(1, uncSqUp[5], uncSqDown[5]);
    // Lepton scale factor uncertainties (uncertainty 6)
    EvaluateLeptons(1, uncSqUp[6], uncSqDown[6]);
}

std::string ControlPlot::Print(int iYear, int iSel) {
    std::stringstream output;

    // prepare data
    double* data = groups[0][iYear][iSel];
    std::vector<double*> predictions;
    for(int iGrp=1; iGrp<nGroups; iGrp++) predictions.push_back(groups[iGrp][iYear][iSel]);
    double theUncSqUp[nBins], theUncSqDown[nBins];
    for(int iBin=0; iBin<nBins; iBin++) {
        theUncSqUp[iBin] = 0.0;
        theUncSqDown[iBin] = 0.0;
        for(int iUnc=0; iUnc<nUncertainties; iUnc++) {
            theUncSqUp[iBin] += uncSqUp[iUnc][iYear][iSel][iBin];
            theUncSqDown[iBin] += uncSqDown[iUnc][iYear][iSel][iBin];
        }
    }

    // print data fields
    output << "$data_single << EOD" << std::endl
           << PrintSingle(data, predictions, theUncSqUp, theUncSqDown)
           << "EOD" << std::endl << std::endl;
    output << "$data_double << EOD" << std::endl
           << PrintDouble(data, predictions, theUncSqUp, theUncSqDown)
           << "EOD" << std::endl << std::endl;

    // print options
    output << "xmin = " << binBoundaries[0] << std::endl;
    output << "xmax = " << binBoundaries[nBins] << std::endl;
    double ymax = 0.0;
    for(int iBin=0; iBin<nBins; iBin++) {
        double thisy = 0.0;
        for(int iGrp=1; iGrp<nGroups; iGrp++) thisy += groups[iGrp][iYear][iSel][iBin];
        thisy += sqrt(theUncSqUp[iBin]);
        if(thisy>ymax) ymax = thisy;
        double datay = groups[0][iYear][iSel][iBin];
        datay += sqrt(datay);
        if(datay>ymax) ymax = datay;
    }
    output << "ymax = " << ymax << std::endl;

    // return result
    return output.str();
}
