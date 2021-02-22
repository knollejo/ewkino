#include <algorithm>
#include <cmath>
#include <sstream>

#include "../interface/HistogramCombiner.h"

HistogramCombiner::HistogramCombiner(int nSel, int nGrp, int nBin) :
    nSelections(nSel), nGroups(nGrp), nBins(nBin),
    binCenters(new double[nBins]), binBoundaries(new double[nBins+1])
{}

HistogramCombiner::~HistogramCombiner() {
    delete [] binCenters;
    delete [] binBoundaries;
}

void HistogramCombiner::SetBins(double* centers, double* boundaries) {
    for(int iBin=0; iBin<nBins; iBin++) {
        binCenters[iBin] = centers[iBin];
        binBoundaries[iBin] = boundaries[iBin];
    }
    binBoundaries[nBins] = boundaries[nBins];
}

void HistogramCombiner::AddSample(double*** values, int year, int group, bool is_data, double weight) {
    samples.push_back(std::make_tuple(values, year, group, is_data, weight));
}

void HistogramCombiner::EvaluateGroup(int iGroup, int iSys, double*** values) {
    for(auto const& sample: samples) {
        if(std::get<2>(sample)!=iGroup) continue;
        const int iYear = std::get<1>(sample);
        const bool isData = std::get<3>(sample);
        const double weight = std::get<4>(sample);
        for(int iSel=0; iSel<nSelections; iSel++) {
            double** input = std::get<0>(sample)[iSel];
            for(int iBin=0; iBin<nBins; iBin++) {
                double value;
                if(iSys==0) { // evaluate statistical uncertainty
                    const double val = weight*input[1][iBin];
                    const double stat = input[0][iBin];
                    if(stat>0.0) value = val*val/stat;
                    else value = 0.0;
                } else if(isData && iSys>1) { // no systematic variations of data
                    value = weight*input[1][iBin];
                } else { // value at specified variation
                    value = weight*input[iSys][iBin];
                }
                values[iYear][iSel][iBin] += value;
                values[nYears][iSel][iBin] += value;
            }
        }
    }
}

void HistogramCombiner::EvaluateEnvelopeUncorrelated(int iYear, int startGroup, std::vector<int> vecSys, double*** valuesUp, double*** valuesDown) {
    const int nSys = vecSys.size();
    double*** varied = new double**[nSelections];
    for(int iSel=0; iSel<nSelections; iSel++) {
        varied[iSel] = new double*[nBins];
        for(int iBin=0; iBin<nBins; iBin++) {
            varied[iSel][iBin] = new double[nSys+1];
            for(int iSys=0; iSys<=nSys; iSys++) {
                varied[iSel][iBin][iSys] = 0.0;
            }
        }
    }
    for(auto const& sample: samples) {
        if(std::get<1>(sample)!=iYear) continue;
        if(std::get<2>(sample)<startGroup) continue;
        const bool isData = std::get<3>(sample);
        const double weight = std::get<4>(sample);
        for(int iSel=0; iSel<nSelections; iSel++) {
            double** input = std::get<0>(sample)[iSel];
            for(int i_iSys=0; i_iSys<=nSys; i_iSys++) {
                const int iSys = (isData || i_iSys==0) ? 1 : vecSys.at(i_iSys-1);
                for(int iBin=0; iBin<nBins; iBin++) {
                    const double val = weight*input[iSys][iBin];
                    varied[iSel][iBin][i_iSys] += val;
                }

            }
        }
    }
    for(int iSel=0; iSel<nSelections; iSel++) {
        for(int iBin=0; iBin<nBins; iBin++) {
            const double central = varied[iSel][iBin][0];
            const double maximum = *std::max_element(varied[iSel][iBin], varied[iSel][iBin]+(nSys+1));
            const double minimum = *std::min_element(varied[iSel][iBin], varied[iSel][iBin]+(nSys+1));
            valuesUp[iYear][iSel][iBin] += (maximum-central)*(maximum-central);
            valuesDown[iYear][iSel][iBin] += (central-minimum)*(central-minimum);
            valuesUp[nYears][iSel][iBin] += (maximum-central)*(maximum-central);
            valuesDown[nYears][iSel][iBin] += (central-minimum)*(central-minimum);
            delete [] varied[iSel][iBin];
        }
        delete [] varied[iSel];
    }
    delete [] varied;
}

void HistogramCombiner::EvaluateEnvelopeCorrelated(int startGroup, std::vector<int> vecSys, double*** valuesUp, double*** valuesDown) {
    const int nSys = vecSys.size();
    double**** varied = new double***[nYears+1];
    for(int iYear=0; iYear<=nYears; iYear++) {
        varied[iYear] = new double**[nSelections];
        for(int iSel=0; iSel<nSelections; iSel++) {
            varied[iYear][iSel] = new double*[nBins];
            for(int iBin=0; iBin<nBins; iBin++) {
                varied[iYear][iSel][iBin] = new double[nSys+1];
                for(int iSys=0; iSys<=nSys; iSys++) {
                    varied[iYear][iSel][iBin][iSys] = 0.0;
                }
            }
        }
    }
    for(auto const& sample: samples) {
        if(std::get<2>(sample)<startGroup) continue;
        const int iYear = std::get<1>(sample);
        const bool isData = std::get<3>(sample);
        const double weight = std::get<4>(sample);
        for(int iSel=0; iSel<nSelections; iSel++) {
            double** input = std::get<0>(sample)[iSel];
            for(int i_iSys=0; i_iSys<=nSys; i_iSys++) {
                const int iSys = (isData || i_iSys==0) ? 1 : vecSys.at(i_iSys-1);
                for(int iBin=0; iBin<nBins; iBin++) {
                    const double val = weight*input[iSys][iBin];
                    varied[iYear][iSel][iBin][i_iSys] += val;
                    varied[nYears][iSel][iBin][i_iSys] += val;
                }

            }
        }
    }
    for(int iYear=0; iYear<=nYears; iYear++) {
        for(int iSel=0; iSel<nSelections; iSel++) {
            for(int iBin=0; iBin<nBins; iBin++) {
                const double central = varied[iYear][iSel][iBin][0];
                const double maximum = *std::max_element(varied[iYear][iSel][iBin], varied[iYear][iSel][iBin]+(nSys+1));
                const double minimum = *std::min_element(varied[iYear][iSel][iBin], varied[iYear][iSel][iBin]+(nSys+1));
                valuesUp[iYear][iSel][iBin] += (maximum-central)*(maximum-central);
                valuesDown[iYear][iSel][iBin] += (central-minimum)*(central-minimum);
                delete [] varied[iYear][iSel][iBin];
            }
            delete [] varied[iYear][iSel];
        }
        delete [] varied[iYear];
    }
    delete [] varied;
}

void HistogramCombiner::EvaluateCentral(int iGroup, double*** values) {
    EvaluateGroup(iGroup, 1, values);
}

void HistogramCombiner::EvaluateStatistical(int iGroup, double*** values) {
    EvaluateGroup(iGroup, 0, values);
}

void HistogramCombiner::EvaluateJetMet(int startGroup, double*** valuesUp, double*** valuesDown) {
    // JEC
    EvaluateEnvelopeUncorrelated(0, startGroup, std::vector<int>({2,3}), valuesUp, valuesDown);
    EvaluateEnvelopeUncorrelated(1, startGroup, std::vector<int>({2,3}), valuesUp, valuesDown);
    EvaluateEnvelopeUncorrelated(2, startGroup, std::vector<int>({2,3}), valuesUp, valuesDown);
    // JER
    EvaluateEnvelopeUncorrelated(0, startGroup, std::vector<int>({4,5}), valuesUp, valuesDown);
    EvaluateEnvelopeUncorrelated(1, startGroup, std::vector<int>({4,5}), valuesUp, valuesDown);
    EvaluateEnvelopeUncorrelated(2, startGroup, std::vector<int>({4,5}), valuesUp, valuesDown);
    // Unclustered
    EvaluateEnvelopeUncorrelated(0, startGroup, std::vector<int>({6,7}), valuesUp, valuesDown);
    EvaluateEnvelopeUncorrelated(1, startGroup, std::vector<int>({6,7}), valuesUp, valuesDown);
    EvaluateEnvelopeUncorrelated(2, startGroup, std::vector<int>({6,7}), valuesUp, valuesDown);
}

void HistogramCombiner::EvaluateBtagging(int startGroup, double*** valuesUp, double*** valuesDown) {
    // heavy flavours
    EvaluateEnvelopeUncorrelated(0, startGroup, std::vector<int>({16,17}), valuesUp, valuesDown);
    EvaluateEnvelopeUncorrelated(1, startGroup, std::vector<int>({16,17}), valuesUp, valuesDown);
    EvaluateEnvelopeUncorrelated(2, startGroup, std::vector<int>({16,17}), valuesUp, valuesDown);
    // light flavours
    EvaluateEnvelopeUncorrelated(0, startGroup, std::vector<int>({18,19}), valuesUp, valuesDown);
    EvaluateEnvelopeUncorrelated(1, startGroup, std::vector<int>({18,19}), valuesUp, valuesDown);
    EvaluateEnvelopeUncorrelated(2, startGroup, std::vector<int>({18,19}), valuesUp, valuesDown);
}

void HistogramCombiner::EvaluateLeptons(int startGroup, double*** valuesUp, double*** valuesDown) {
    // muon ID
    EvaluateEnvelopeCorrelated(startGroup, std::vector<int>({22,23}), valuesUp, valuesDown);
    // electron reco
    EvaluateEnvelopeCorrelated(startGroup, std::vector<int>({24,25}), valuesUp, valuesDown);
    // electron ID
    EvaluateEnvelopeCorrelated(startGroup, std::vector<int>({26,27}), valuesUp, valuesDown);
}

void HistogramCombiner::EvaluatePileup(int startGroup, double*** valuesUp, double*** valuesDown) {
    // pileup reweighting
    EvaluateEnvelopeCorrelated(startGroup, std::vector<int>({14,15}), valuesUp, valuesDown);
}

void HistogramCombiner::EvaluateMeScales(int startGroup, double*** valuesUp, double*** valuesDown) {
    // ME scale variations
    EvaluateEnvelopeCorrelated(startGroup, std::vector<int>({8,9,10,11,12,13}), valuesUp, valuesDown);
}

std::string HistogramCombiner::PrintSingle(double* data, std::vector<double*> predictions, double* uncSqUp, double* uncSqDown) {
    const int nPredictions = predictions.size();
    std::stringstream output;
    for(int iBin=0; iBin<nBins; iBin++) {
        output << iBin << "  " << binCenters[iBin] << "  " << data[iBin];
        for(int iPred=0; iPred<nPredictions; iPred++) {
            double val = 0.0;
            for(int jPred=iPred; jPred<nPredictions; jPred++) val += predictions.at(jPred)[iBin];
            output << "  " << val;
        }
        output << "  " << sqrt(uncSqUp[iBin]) << "  " << sqrt(uncSqDown[iBin]) << std::endl;
    }
    return output.str();
}

std::string HistogramCombiner::PrintDouble(double* data, std::vector<double*> predictions, double* uncSqUp, double* uncSqDown) {
    const int nPredictions = predictions.size();
    std::stringstream output, row;
    for(int iBin=0; iBin<nBins; iBin++) {
        row.str("");
        row << data[iBin];
        for(int iPred=0; iPred<nPredictions; iPred++) {
            double val = 0.0;
            for(int jPred=iPred; jPred<nPredictions; jPred++) val += predictions.at(jPred)[iBin];
            row << "  " << val;
        }
        row << "  " << sqrt(uncSqUp[iBin]) << "  " << sqrt(uncSqDown[iBin]);
        output << iBin << "  " << binBoundaries[iBin] << "  " << row.str() << std::endl
               << iBin << "  " << binBoundaries[iBin+1] << "  " << row.str() << std::endl;
    }
    return output.str();
}
