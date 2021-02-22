#ifndef HistogramFile_H
#define HistogramFile_H

#include <string>

#include <TDirectory.h>

class HistogramFile {
protected:
    HistogramFile(int, int, int, int);
    void readDirectory(TDirectory*, std::string, double);

public:
    HistogramFile(TDirectory*, std::string, int, int, int, int, double);
    HistogramFile(std::string, std::string, int, int, int, int, double);
    ~HistogramFile();

    const int nVariants, nSelections, nSystematics, nBins;
    double**** values;
    double* binCenters, * binBoundaries;
};

#endif
