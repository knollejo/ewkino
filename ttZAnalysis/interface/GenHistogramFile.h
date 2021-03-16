#ifndef GenHistogramFile_H
#define GenHistogramFile_H

#include <string>
#include <utility>
#include <vector>

#include <TDirectory.h>

class GenHistogramFile {
protected:
    GenHistogramFile(int, int, int, int);

public:
    typedef std::pair<std::vector<int>, std::vector<std::pair<int, bool>>> Definition;
    typedef std::vector<Definition> Definitions;

    void readDirectory(TDirectory*, const Definitions, double);
    void readFile(std::string, std::string, const Definitions, double);
    GenHistogramFile(TDirectory*, int, const Definitions, int, int, int, double);
    GenHistogramFile(std::string, std::string, int, const Definitions, int, int, int, double);
    ~GenHistogramFile();

    const int nDefinitions, nSystematics, nBinsX, nBinsY;
    double*** values_gen, *** values_fail, **** values_response;
    double** binCenters, ** binBoundaries;
};

#endif
