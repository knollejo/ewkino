#include <iostream>

#include <TFile.h>
#include <TH1.h>

#include "../interface/HistogramFile.h"

HistogramFile::HistogramFile(int nVar, int nSel, int nSys, int nBin) :
    nVariants(nVar), nSelections(nSel), nSystematics(nSys), nBins(nBin),
    values(new double***[nVariants]),
    binCenters(new double[nBins]), binBoundaries(new double[nBins+1])
{
    for(int iVar=0; iVar<nVariants; iVar++) {
        values[iVar] = new double**[nSelections];
        for(int iSel=0; iSel<nSelections; iSel++) {
            values[iVar][iSel] = new double*[nSystematics+2];
            for(int iSys=0; iSys<nSystematics+2; iSys++) {
                values[iVar][iSel][iSys] = new double[nBins];
            }
        }
    }
}

void HistogramFile::readDirectory(TDirectory* dir, std::string obsname, double xsec) {
    TH1* hist = nullptr;
    for(int iVar=0; iVar<nVariants; iVar++) {
        for(int iSel=0; iSel<nSelections; iSel++) {
            for(int iSys=0; iSys<nSystematics+2; iSys++) {
                hist = (TH1*)dir->Get((
                    obsname
                    +"_var"+std::to_string(iVar)
                    +"_sel"+std::to_string(iSel)
                    +(iSys==0 ? std::string("_stat") : "_sys"+std::to_string(iSys-1))
                ).c_str());
                const double factor = iSys==0 ? 1.0 : xsec;
                for(int iBin=0; iBin<nBins; iBin++) {
                    values[iVar][iSel][iSys][iBin] = factor*hist->GetBinContent(iBin+1);
                }
            }
        }
    }
    for(int iBin=0; iBin<nBins; iBin++) {
        binCenters[iBin] = hist->GetBinCenter(iBin+1);
        binBoundaries[iBin] = hist->GetBinLowEdge(iBin+1);
    }
    binBoundaries[nBins] = hist->GetXaxis()->GetBinUpEdge(nBins);
}

HistogramFile::HistogramFile(TDirectory* dir, std::string obsname, int nVar, int nSel, int nSys, int nBin, double xsec) :
    HistogramFile(nVar, nSel, nSys, nBin)
{
    readDirectory(dir, obsname, xsec);
}

HistogramFile::HistogramFile(std::string filename, std::string obsname, int nVar, int nSel, int nSys, int nBin, double xsec) :
    HistogramFile(nVar, nSel, nSys, nBin)
{
    std::cout << "Reading " << filename << std::endl;
    TFile* tfile = TFile::Open(filename.c_str());
    readDirectory((TDirectory*)tfile->Get(obsname.c_str()), obsname, xsec);
    tfile->Close();
}

HistogramFile::~HistogramFile() {
    for(int iVar=0; iVar<nVariants; iVar++) {
        for(int iSel=0; iSel<nSelections; iSel++) {
            for(int iSys=0; iSys<nSystematics+2; iSys++) {
                delete [] values[iVar][iSel][iSys];
            }
            delete [] values[iVar][iSel];
        }
        delete [] values[iVar];
    }
    delete [] values;
    delete [] binCenters;
    delete [] binBoundaries;
}
