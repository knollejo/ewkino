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
    uncSq = new double***[2];
    for(int iUnc=0; iUnc<2; iUnc++) {
        uncSq[iUnc] = new double**[nYears+1];
        for(int iYear=0; iYear<=nYears; iYear++) {
            uncSq[iUnc][iYear] = new double*[nSelections];
            for(int iSel=0; iSel<nSelections; iSel++) {
                uncSq[iUnc][iYear][iSel] = new double[nBins];
                for(int iBin=0; iBin<nBins; iBin++) {
                    uncSq[iUnc][iYear][iSel][iBin] = 0.0;
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
    for(int iUnc=0; iUnc<2; iUnc++) {
        for(int iYear=0; iYear<=nYears; iYear++) {
            for(int iSel=0; iSel<nSelections; iSel++) {
                delete [] uncSq[iUnc][iYear][iSel];
            }
            delete [] uncSq[iUnc][iYear];
        }
        delete [] uncSq[iUnc];
    }
    delete [] uncSq;
}

void ControlPlot::Evaluate() {
    std::cout << "Evaluating" << std::endl;
    // central values
    for(int iGroup=0; iGroup<nGroups; iGroup++) {
        EvaluateCentral(iGroup, groups[iGroup]);
    }
    // statistical uncertainties
    for(int iGroup=1; iGroup<nGroups; iGroup++) {
        EvaluateStatistical(iGroup, uncSq[0]);
        EvaluateStatistical(iGroup, uncSq[1]);
    }
    // Experimental uncertainties
    // EvaluateJetMet(1, uncSq[0], uncSq[1]);
    // EvaluateBtagging(1, uncSq[0], uncSq[1]);
    // EvaluateLeptons(1, uncSq[0], uncSq[1]);
    // EvaluatePileup(1, uncSq[0], uncSq[1]);
    // Theory uncertainties
    // EvaluateMeScales(1, uncSq[0], uncSq[1]);
}

std::string ControlPlot::Print(int iYear, int iSel) {
    std::stringstream output;

    // prepare data
    double* data = groups[0][iYear][iSel];
    std::vector<double*> predictions;
    for(int iGrp=1; iGrp<nGroups; iGrp++) predictions.push_back(groups[iGrp][iYear][iSel]);
    double* uncSqUp = uncSq[0][iYear][iSel];
    double* uncSqDown = uncSq[1][iYear][iSel];

    // print data fields
    output << "$data_single << EOD" << std::endl
           << PrintSingle(data, predictions, uncSqUp, uncSqDown)
           << "EOD" << std::endl << std::endl;
    output << "$data_double << EOD" << std::endl
           << PrintDouble(data, predictions, uncSqUp, uncSqDown)
           << "EOD" << std::endl << std::endl;

    // print options
    output << "xmin = " << binBoundaries[0] << std::endl;
    output << "xmax = " << binBoundaries[nBins] << std::endl;
    double ymax = 0.0;
    for(int iBin=0; iBin<nBins; iBin++) {
        double thisy = 0.0;
        for(int iGrp=1; iGrp<nGroups; iGrp++) thisy += groups[iGrp][iYear][iSel][iBin];
        thisy += sqrt(uncSq[0][iYear][iSel][iBin]);
        if(thisy>ymax) ymax = thisy;
        double datay = groups[0][iYear][iSel][iBin];
        datay += sqrt(datay);
        if(datay>ymax) ymax = datay;
    }
    output << "ymax = " << ymax << std::endl;

    // return result
    return output.str();
}

void ControlPlot::Dump(int iYear, int iSel) {
    // prepare data
    double* data = groups[0][iYear][iSel];
    double* signal = groups[1][iYear][iSel];
    double* ttX = groups[2][iYear][iSel];
    double* nonprompt = groups[6][iYear][iSel];
    double* uncSqUp = uncSq[0][iYear][iSel];
    double* uncSqDown = uncSq[1][iYear][iSel];

    // print data
    for(int iBin=0; iBin<nBins; iBin++) {
        std::cout << iBin << binCenters[iBin] << " " << data[iBin] << " "
                  << signal[iBin] << " " << ttX[iBin] << " "
                  << nonprompt[iBin] << std::endl;
    }
    std::cout << std::endl << std::endl;
    std::cout << PrintSingle(data, std::vector<double*>({signal, ttX, nonprompt}), uncSqUp, uncSqDown);
}
