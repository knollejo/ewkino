#ifndef HistogramFile_H
#define HistogramFile_H

#include <string>
#include <utility>
#include <vector>

#include <TDirectory.h>

class HistogramFile {
protected:
    HistogramFile(int, int, int, int);

public:
    void readDirectory(TDirectory*, std::vector<std::pair<int, int>>, double);
    void readFile(std::string, std::string, std::vector<std::pair<int, int>>, double);
    HistogramFile(TDirectory*, int, int, std::vector<std::pair<int, int>>, int, int, double);
    HistogramFile(std::string, std::string, int, int, std::vector<std::pair<int, int>>, int, int, double);
    ~HistogramFile();

    const int nVariants, nSelections, nSystematics, nBins;
    double**** values;
    double* binCenters, * binBoundaries;
};

#endif
