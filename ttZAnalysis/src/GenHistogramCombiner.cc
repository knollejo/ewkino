#include <algorithm>
#include <cmath>
#include <sstream>

#include "../interface/GenHistogramCombiner.h"

GenHistogramCombiner::GenHistogramCombiner(int nDef, int nSys, int nBinX, int nBinY) :
    nDefinitions(nDef), nSystematics(nSys), nBinsX(nBinX), nBinsY(nBinY),
    responseMatrices(new TH2F***[nDef]),
    genDistributions(new double**[nDef]), genUncertainties(new double***[nDef])
{
    for(int iDef=0; iDef<nDefinitions; iDef++) {
        responseMatrices[iDef] = new TH2F**[nYears+1];
        genDistributions[iDef] = new double*[nYears+1];
        genUncertainties[iDef] = new double**[2];
        genUncertainties[iDef][0] = new double*[nYears+1];
        genUncertainties[iDef][1] = new double*[nYears+1];
        for(int iYear=0; iYear<=nYears; iYear++) {
            responseMatrices[iDef][iYear] = new TH2F*[nSystematics+1];
            for(int iSys=0; iSys<nSystematics+1; iSys++) {
                responseMatrices[iDef][iYear][iSys] = new TH2F(
                    ("resp_"+std::to_string(iDef)+"_"+std::to_string(iYear)+"_"+std::to_string(iSys)).c_str(), "",
                    nBinsX, 0.0, 1.0, nBinsY+1, 0.0, 1.0
                );
            }
            genDistributions[iDef][iYear] = new double[nBinsX];
            genUncertainties[iDef][0][iYear] = new double[nBinsX];
            genUncertainties[iDef][1][iYear] = new double[nBinsX];
        }
    }
}

GenHistogramCombiner::~GenHistogramCombiner() {
    for(int iDef=0; iDef<nDefinitions; iDef++) {
        for(int iYear=0; iYear<=nYears; iYear++) {
            for(int iSys=0; iSys<nSystematics+1; iSys++) {
                delete responseMatrices[iDef][iYear][iSys];
            }
            delete [] genDistributions[iDef][iYear];
            delete [] genUncertainties[iDef][0][iYear];
            delete [] genUncertainties[iDef][1][iYear];
        }
        delete [] responseMatrices[iDef];
        delete [] genDistributions[iDef];
        delete [] genUncertainties[iDef][0];
        delete [] genUncertainties[iDef][1];
        delete [] genUncertainties[iDef];
    }
    delete [] responseMatrices;
    delete [] genDistributions;
    delete [] genUncertainties;
}

void GenHistogramCombiner::SetBins(double** centers, double** boundaries) {
    binCenters = centers;
    binBoundaries = boundaries;
}

void GenHistogramCombiner::AddSample(double** valuesg, double** valuesf, double*** valuesr, int year, int definition, double weight) {
    samples.push_back(std::make_tuple(valuesg, valuesf, valuesr, year, definition, weight));
}

void GenHistogramCombiner::EvaluateGenOrFail(int iDef, int iSys, double** values, bool is_gen) {
    for(auto const& sample: samples) {
        if(std::get<4>(sample)!=iDef) continue;
        const int iYear = std::get<3>(sample);
        const double weight = std::get<5>(sample);
        double** input;
        if(is_gen) input = std::get<0>(sample);
        else input = std::get<1>(sample);
        for(int iBinX=0; iBinX<nBinsX; iBinX++) {
            double value;
            if(iSys==0) { // evaluate statistical uncertainty
                const double val = weight*input[1][iBinX];
                const double stat = input[0][iBinX];
                if(stat>0.0) value = val*val/stat;
                else value = 0.0;
            } else { // value at specified variation
                value = weight*input[iSys][iBinX];
            }
            values[iYear][iBinX] += value;
            values[nYears][iBinX] += value;
        }
    }
}

void GenHistogramCombiner::EvaluateResponse(int iDef, int iSys, double*** values) {
    for(auto const& sample: samples) {
        if(std::get<4>(sample)!=iDef) continue;
        const int iYear = std::get<3>(sample);
        const double weight = std::get<5>(sample);
        double*** input = std::get<2>(sample);
        for(int iBinX=0; iBinX<nBinsX; iBinX++) {
            for(int iBinY=0; iBinY<nBinsY; iBinY++) {
                double value;
                if(iSys==0) { // evaluate statistical uncertainty
                    const double val = weight*input[1][iBinX][iBinY];
                    const double stat = input[0][iBinX][iBinY];
                    if(stat>0.0) value = val*val/stat;
                    else value = 0.0;
                } else { // value at specified variation
                    value = weight*input[iSys][iBinX][iBinY];
                }
                values[iYear][iBinX][iBinY] += value;
                values[nYears][iBinX][iBinY] += value;
            }
        }
    }
}

void GenHistogramCombiner::Evaluate() {
    // evaluate response matrices
    double*** values2d = new double**[nYears+1];
    double** values1d = new double*[nYears+1];
    for(int iYear=0; iYear<=nYears; iYear++) {
        values2d[iYear] = new double*[nBinsX];
        values1d[iYear] = new double[nBinsX];
        for(int iBinX=0; iBinX<nBinsX; iBinX++) {
            values2d[iYear][iBinX] = new double[nBinsY];
        }
    }
    for(int iDef=0; iDef<nDefinitions; iDef++) {
        for(int iSys=0; iSys<nSystematics+1; iSys++) {
            // reset values to zero
            for(int iYear=0; iYear<=nYears; iYear++) {
                for(int iBinX=0; iBinX<nBinsX; iBinX++) {
                    values1d[iYear][iBinX] = 0.0;
                    for(int iBinY=0; iBinY<nBinsY; iBinY++) {
                        values2d[iYear][iBinX][iBinY] = 0.0;
                    }
                }
            }
            // evaluate and save values for given variation
            EvaluateResponse(iDef, iSys+1, values2d);
            EvaluateFail(iDef, iSys+1, values1d);
            for(int iYear=0; iYear<=nYears; iYear++) {
                for(int iBinX=0; iBinX<nBinsX; iBinX++) {
                    responseMatrices[iDef][iYear][iSys]->SetBinContent(responseMatrices[iDef][iYear][iSys]->GetBin(iBinX+1, 0), values1d[iYear][iBinX]);
                    for(int iBinY=0; iBinY<nBinsY; iBinY++) {
                        responseMatrices[iDef][iYear][iSys]->SetBinContent(responseMatrices[iDef][iYear][iSys]->GetBin(iBinX+1, iBinY+1), values2d[iYear][iBinX][iBinY]);
                    }
                }
            }
            if(iSys==0) { // for central variation: additionally evaluate statistical uncertainty
                for(int iYear=0; iYear<=nYears; iYear++) {
                    for(int iBinX=0; iBinX<nBinsX; iBinX++) {
                        values1d[iYear][iBinX] = 0.0;
                        for(int iBinY=0; iBinY<nBinsY; iBinY++) {
                            values2d[iYear][iBinX][iBinY] = 0.0;
                        }
                    }
                }
                EvaluateResponse(iDef, 0, values2d);
                EvaluateFail(iDef, 0, values1d);
                for(int iYear=0; iYear<=nYears; iYear++) {
                    for(int iBinX=0; iBinX<nBinsX; iBinX++) {
                        responseMatrices[iDef][iYear][iSys]->SetBinError(responseMatrices[iDef][iYear][iSys]->GetBin(iBinX+1, 0), sqrt(values1d[iYear][iBinX]));
                        for(int iBinY=0; iBinY<nBinsY; iBinY++) {
                            responseMatrices[iDef][iYear][iSys]->SetBinError(responseMatrices[iDef][iYear][iSys]->GetBin(iBinX+1, iBinY+1), sqrt(values2d[iYear][iBinX][iBinY]));
                        }
                    }
                }
            } else { // for other variations: use statistical uncertainty of central variation
                for(int iYear=0; iYear<=nYears; iYear++) {
                    for(int iBinX=0; iBinX<nBinsX; iBinX++) {
                        const int failbin = responseMatrices[iDef][iYear][iSys]->GetBin(iBinX+1, 0);
                        responseMatrices[iDef][iYear][iSys]->SetBinError(failbin, responseMatrices[iDef][iYear][0]->GetBinError(failbin));
                        for(int iBinY=0; iBinY<nBinsY; iBinY++) {
                            const int thebin = responseMatrices[iDef][iYear][iSys]->GetBin(iBinX+1, iBinY+1);
                            responseMatrices[iDef][iYear][iSys]->SetBinError(thebin, responseMatrices[iDef][iYear][0]->GetBinError(failbin));
                        }
                    }
                }
            }
        }
    }
    for(int iYear=0; iYear<=nYears; iYear++) {
        for(int iBinX=0; iBinX<nBinsX; iBinX++) {
            delete [] values2d[iYear][iBinX];
        }
        delete [] values2d[iYear];
        delete [] values1d[iYear];
    }
    delete [] values2d;
    delete [] values1d;
    // evaluate generated distributions
    for(int iDef=0; iDef<nDefinitions; iDef++) {
        // reset values to zero
        for(int iYear=0; iYear<=nYears; iYear++) {
            for(int iBinX=0; iBinX<nBinsX; iBinX++) {
                genDistributions[iDef][iYear][iBinX] = 0.0;
                genUncertainties[iDef][0][iYear][iBinX] = 0.0;
                genUncertainties[iDef][1][iYear][iBinX] = 0.0;
            }
        }
        EvaluateGen(iDef, 1, genDistributions[iDef]);
        EvaluateGen(iDef, 0, genUncertainties[iDef][0]);
        EvaluateGen(iDef, 0, genUncertainties[iDef][1]);
    }
}

std::string GenHistogramCombiner::PrintResponse(int iYear, int iDef) {
    std::stringstream output;

    // prepare data
    TH2F* response = responseMatrices[iDef][iYear][0];
    double data[nBinsX][nBinsY];
    double uncertainty[nBinsX][nBinsY];
    for(int iBinX=0; iBinX<nBinsX; iBinX++) {
        double sum = 0.0;
        for(int iBinY=0; iBinY<=nBinsY; iBinY++) {
            sum += response->GetBinContent(response->GetBin(iBinX+1, iBinY));
        }
        for(int iBinY=0; iBinY<nBinsY; iBinY++) {
            const int iBin = response->GetBin(iBinX+1, iBinY+1);
            data[iBinX][iBinY] = sum>0.0 ? response->GetBinContent(iBin)/sum : 0.0;
            uncertainty[iBinX][iBinY] = sum>0.0 ? response->GetBinError(iBin)/sum : 0.0;
        }
    }

    // print data fields
    output << "$matrix << EOD" << std::endl;
    for(int iBinY=0; iBinY<nBinsY; iBinY++) {
        for(int iBinX=0; iBinX<nBinsX; iBinX++) {
            output << "  " << data[iBinX][iBinY];
        }
        output << std::endl;
    }
    output << "EOD" << std::endl << std::endl;
    output << "$values << EOD" << std::endl;
    for(int iBinX=0; iBinX<nBinsX; iBinX++) {
        for(int iBinY=0; iBinY<nBinsY; iBinY++) {
            output << "  " << iBinX << "  " << iBinY
                   << "  " << data[iBinX][iBinY] << "  " << uncertainty[iBinX][iBinY]
                   << std::endl;
        }
    }
    output << "EOD" << std::endl << std::endl;

    // print options
    output << "nBinsRec = " << nBinsY << std::endl;
    output << "nBinsGen = " << nBinsX << std::endl;
    output << "binsRec = '(";
    for(int iBinY=0; iBinY<nBinsY; iBinY+=2) {
        output << '"' << binBoundaries[1][iBinY] << "\" " << iBinY-0.5 << ", ";
    }
    output << '"' << binBoundaries[1][nBinsY] << "\" " << nBinsY-0.5 << ")'" << std::endl;
    output << "binsGen = '(";
    for(int iBinX=0; iBinX<nBinsX; iBinX++) {
        output << '"' << binBoundaries[0][iBinX] << "\" " << iBinX-0.5 << ", ";
    }
    output << '"' << binBoundaries[0][nBinsX] << "\" " << nBinsX-0.5 << ")'" << std::endl;

    // return result
    return output.str();
}

std::string GenHistogramCombiner::PrintStabPurEff(int iYear, int iDef) {
    std::stringstream output;

    // prepare data
    TH2F* response = responseMatrices[iDef][iYear][0];
    const int nRebin = nBinsY/nBinsX;
    // calculate efficiency
    double eff[nBinsX];
    for(int iBinX=0; iBinX<nBinsX; iBinX++) {
        double sum = 0.0;
        for(int iBinY=0; iBinY<nBinsY; iBinY++) {
            sum += response->GetBinContent(response->GetBin(iBinX+1, iBinY+1));
        }
        eff[iBinX] = sum>0.0 ? sum/(sum+response->GetBinContent(response->GetBin(iBinX+1, 0))) : 0.0;
    }
    // calculate stability
    double stab[nBinsX];
    for(int iBinX=0; iBinX<nBinsX; iBinX++) {
        double sum = 0.0;
        for(int iBinY=0; iBinY<nBinsY; iBinY++) {
            sum += response->GetBinContent(response->GetBin(iBinX+1, iBinY+1));
        }
        double selected = 0.0;
        for(int iBinY=iBinX*nRebin; iBinY<(iBinX+1)*nRebin; iBinY++) {
            selected += response->GetBinContent(response->GetBin(iBinX+1, iBinY+1));
        }
        stab[iBinX] = sum>0.0 ? selected/sum : 0.0;
    }
    // calculate purity
    double pur[nBinsY];
    for(int iBinY=0; iBinY<nBinsY; iBinY++) {
        double sum = 0.0;
        for(int iBinX=0; iBinX<nBinsX; iBinX++) {
            sum += response->GetBinContent(response->GetBin(iBinX+1, iBinY+1));
        }
        pur[iBinY] = sum>0.0 ? response->GetBinContent(response->GetBin(iBinY/nRebin+1, iBinY+1))/sum : 0.0;
    }

    // print data fields
    output << "$stabeff_single << EOD" << std::endl;
    for(int iBinX=0; iBinX<nBinsX; iBinX++) {
        output << "  " << binCenters[0][iBinX]
               << "  " << stab[iBinX] << "  " << eff[iBinX]
               << std::endl;
    }
    output << "EOD" << std::endl << std::endl;
    output << "$stabeff_double << EOD" << std::endl;
    for(int iBinX=0; iBinX<nBinsX; iBinX++) {
        for(int iDouble=0; iDouble<2; iDouble++) {
            output << "  " << binBoundaries[0][iBinX+iDouble]
                   << "  " << stab[iBinX] << "  " << eff[iBinX]
                   << std::endl;
        }
    }
    output << "EOD" << std::endl << std::endl;
    output << "$purity_single << EOD" << std::endl;
    for(int iBinY=0; iBinY<nBinsY; iBinY++) {
        output << "  " << binCenters[1][iBinY]
               << "  " << pur[iBinY]
               << std::endl;
    }
    output << "EOD" << std::endl << std::endl;
    output << "$purity_double << EOD" << std::endl;
    for(int iBinY=0; iBinY<nBinsY; iBinY++) {
        for(int iDouble=0; iDouble<2; iDouble++) {
            output << "  " << binBoundaries[1][iBinY+iDouble]
                   << "  " << pur[iBinY]
                   << std::endl;
        }
    }
    output << "EOD" << std::endl << std::endl;

    // print options
    output << "xmin = " << binBoundaries[0][0] << std::endl;
    output << "xmax = " << binBoundaries[0][nBinsX] << std::endl;

    // return result
    return output.str();
}
