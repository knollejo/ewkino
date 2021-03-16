#include <iostream>

#include <TFile.h>
#include <TH1.h>
#include <TH2.h>

#include "../interface/GenHistogramFile.h"

GenHistogramFile::GenHistogramFile(int nDef, int nSys, int nBinX, int nBinY) :
    nDefinitions(nDef), nSystematics(nSys), nBinsX(nBinX), nBinsY(nBinY),
    values_gen(new double**[nDefinitions]), values_fail(new double**[nDefinitions]),
    values_response(new double***[nDefinitions]),
    binCenters(new double*[2]), binBoundaries(new double*[2])
{
    for(int iDef=0; iDef<nDefinitions; iDef++) {
        values_gen[iDef] = new double*[nSystematics+2];
        values_fail[iDef] = new double*[nSystematics+2];
        values_response[iDef] = new double**[nSystematics+2];
        for(int iSys=0; iSys<nSystematics+2; iSys++) {
            values_gen[iDef][iSys] = new double[nBinsX];
            values_fail[iDef][iSys] = new double[nBinsX];
            values_response[iDef][iSys] = new double*[nBinsX];
            for(int iBinX=0; iBinX<nBinsX; iBinX++) {
                values_gen[iDef][iSys][iBinX] = 0.0;
                values_fail[iDef][iSys][iBinX] = 0.0;
                values_response[iDef][iSys][iBinX] = new double[nBinsY];
                for(int iBinY=0; iBinY<nBinsY; iBinY++) {
                    values_response[iDef][iSys][iBinX][iBinY] = 0.0;
                }
            }
        }
    }
    for(int iDim=0; iDim<2; iDim++) {
        binCenters[iDim] = new double[iDim==0 ? nBinsX : nBinsY];
        binBoundaries[iDim] = new double[(iDim==0 ? nBinsX : nBinsY)+1];
    }
}

void GenHistogramFile::readDirectory(TDirectory* dir, const GenHistogramFile::Definitions definitionMap, double xsec) {
    TH2* histr = nullptr;
    for(int iDef=0; iDef<nDefinitions; iDef++) {
        Definition def = definitionMap.at(iDef);
        for(const auto iVar: def.first) {
            for(int iSys=0; iSys<nSystematics+2; iSys++) {
                const double factor = iSys==0 ? 1.0 : xsec;
                TH1* histg = (TH1*)dir->Get((
                    "gen_var"+std::to_string(iVar)
                    +(iSys==0 ? std::string("_stat") : "_sys"+std::to_string(iSys-1))
                ).c_str());
                TH1* histf = (TH1*)dir->Get((
                    "fail_var"+std::to_string(iVar)
                    +(iSys==0 ? std::string("_stat") : "_sys"+std::to_string(iSys-1))
                ).c_str());
                for(int iBinX=0; iBinX<nBinsX; iBinX++) {
                    values_gen[iDef][iSys][iBinX] += factor*histg->GetBinContent(iBinX+1);
                    values_fail[iDef][iSys][iBinX] += factor*histf->GetBinContent(iBinX+1);
                }
                for(const auto iSel: def.second) {
                    histr = (TH2*)dir->Get((
                        "response_var"+std::to_string(iVar)
                        +"_sel"+std::to_string(iSel.first)
                        +(iSys==0 ? std::string("_stat") : "_sys"+std::to_string(iSys-1))
                    ).c_str());
                    for(int iBinX=0; iBinX<nBinsX; iBinX++) {
                        for(int iBinY=0; iBinY<nBinsY; iBinY++) {
                            const double value = factor*histr->GetBinContent(histr->GetBin(iBinX+1, iBinY+1));
                            if(iSel.second) values_response[iDef][iSys][iBinX][iBinY] += value;
                            else values_fail[iDef][iSys][iBinX] += value;
                        }
                    }
                }
            }
        }
    }
    for(int iDim=0; iDim<2; iDim++) {
        const int nBins = iDim==0 ? nBinsX : nBinsY;
        const auto axis = iDim==0 ? histr->GetXaxis() : histr->GetYaxis();
        for(int iBin=0; iBin<nBins; iBin++) {
            binCenters[iDim][iBin] = axis->GetBinCenter(iBin+1);
            binBoundaries[iDim][iBin] = axis->GetBinLowEdge(iBin+1);
        }
        binBoundaries[iDim][nBins] = axis->GetBinUpEdge(nBins);
    }
}

void GenHistogramFile::readFile(std::string filename, std::string obsname, const GenHistogramFile::Definitions definitionMap, double xsec) {
    std::cout << "Reading " << filename << std::endl;
    TFile* tfile = TFile::Open(filename.c_str());
    readDirectory((TDirectory*)tfile->Get(obsname.c_str()), definitionMap, xsec);
    tfile->Close();
}

GenHistogramFile::GenHistogramFile(TDirectory* dir, int nDef, const GenHistogramFile::Definitions definitionMap, int nSys, int nBinX, int nBinY, double xsec) :
    GenHistogramFile(nDef, nSys, nBinX, nBinY)
{
    readDirectory(dir, definitionMap, xsec);
}

GenHistogramFile::GenHistogramFile(std::string filename, std::string obsname, int nDef, const GenHistogramFile::Definitions definitionMap, int nSys, int nBinX, int nBinY, double xsec) :
    GenHistogramFile(nDef, nSys, nBinX, nBinY)
{
    readFile(filename, obsname, definitionMap, xsec);
}

GenHistogramFile::~GenHistogramFile() {
    for(int iDef=0; iDef<nDefinitions; iDef++) {
        for(int iSys=0; iSys<nSystematics+2; iSys++) {
            delete [] values_gen[iDef][iSys];
            delete [] values_fail[iDef][iSys];
            for(int iBinX=0; iBinX<nBinsX; iBinX++) {
                delete [] values_response[iDef][iSys][iBinX];
            }
            delete [] values_response[iDef][iSys];
        }
        delete [] values_gen[iDef];
        delete [] values_fail[iDef];
        delete [] values_response[iDef];
    }
    delete [] values_gen;
    delete [] values_fail;
    delete [] values_response;
    for(int iDim=0; iDim<2; iDim++) {
        delete [] binCenters[iDim];
        delete [] binBoundaries[iDim];
    }
    delete [] binCenters;
    delete [] binBoundaries;
}
