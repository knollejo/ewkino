#include <fstream>
#include <iostream>
#include <string>
#include <utility>
#include <vector>

#include <TFile.h>

#include "interface/GenHistogramCombiner.h"
#include "interface/CrossSections.h"
#include "interface/GenHistogramFile.h"
#include "interface/PlotAxes.h"

const int nYears = 3;
const int nDefinitions = 3;
std::vector<std::pair<int, bool>> selections3l3j({
    std::make_pair(0, false), // 4l
    std::make_pair(1, true), // 3l3j
    std::make_pair(2, false), // 3l4j
});
std::vector<std::pair<int, bool>> selections3l4j({
    std::make_pair(0, false), // 4l
    std::make_pair(1, false), // 3l3j
    std::make_pair(2, true), // 3l4j
});
std::vector<std::pair<int, bool>> selections3l({
    std::make_pair(0, false), // 4l
    std::make_pair(1, true), // 3l3j
    std::make_pair(2, true), // 3l4j
});
std::vector<std::pair<int, bool>> selections4l({
    std::make_pair(0, true), // 4l
    std::make_pair(1, false), // 3l3j
    std::make_pair(2, false), // 3l4j
});
std::vector<std::pair<int, bool>> selectionsAll({
    std::make_pair(0, true), // 4l
    std::make_pair(1, true), // 3l3j
    std::make_pair(2, true), // 3l4j
});
const GenHistogramFile::Definitions definitions3l({
    std::make_pair(std::vector<int>({0}), selections3l), // 4l -> 3l
    std::make_pair(std::vector<int>({3}), selections3l), // 3l -> 3l
    std::make_pair(std::vector<int>({3}), selections3l3j), // 3l -> 3l3j
    std::make_pair(std::vector<int>({3}), selections3l4j), // 3l -> 3l4j
    std::make_pair(std::vector<int>({9}), selections3l), // Ztau -> 3l
    std::make_pair(std::vector<int>({0, 3, 6, 9, 12}), selections3l), // on-shell -> 3l
});
const std::vector<std::string> sDefinitions3l({"4lto3l", "3lto3l", "3lto3l3j", "3lto3l4j", "tauto3l", "allto3l"});
const GenHistogramFile::Definitions definitions4l({
    std::make_pair(std::vector<int>({0}), selections4l), // 4l -> 4l
    std::make_pair(std::vector<int>({9}), selections4l), // Ztau -> 4l
    std::make_pair(std::vector<int>({0, 9, 12}), selections4l), // 4l+taus -> 4l
});
const std::vector<std::string> sDefinitions4l({"4lto4l", "tauto4l", "allto4l"});
const GenHistogramFile::Definitions definitionsAll({
    std::make_pair(std::vector<int>({0, 3}), selectionsAll), // 3l4l -> 3l4l
    std::make_pair(std::vector<int>({9}), selectionsAll), // Ztau -> 3l4l
    std::make_pair(std::vector<int>({0, 3, 6, 9, 12}), selectionsAll), // on-shell -> 3l4l
});
const std::vector<std::string> sDefinitionsAll({"3l4lto3l4l", "tauto3l4l", "allto4l"});
const int nSys = 126;

const std::string path = "output/210315_095935/";

void makeResponseMatrices(std::string obsname, int nBinsX, int nBinsY, std::string outputdir) {
    const GenHistogramFile::Definitions* definitions = obsname.rfind("response_", 0)==0
                                                     ? &definitionsAll
                                                     : obsname.rfind("response3l_", 0)==0
                                                     ? &definitions3l
                                                     : &definitions4l;
    const std::vector<std::string>* sDefinitions = obsname.rfind("response_", 0)==0
                                                ? &sDefinitionsAll
                                                : obsname.rfind("response3l_", 0)==0
                                                ? &sDefinitions3l
                                                : &sDefinitions4l;
    const int nDefinitions = definitions->size();
    GenHistogramCombiner plot(nDefinitions, nSys, nBinsX, nBinsY);

    #define HISTOGRAM_FILE(OBJNAME, FILENAME, XSEC) \
        GenHistogramFile OBJNAME(path+FILENAME, obsname, nDefinitions, *definitions, nSys, nBinsX, nBinsY, XSEC)
    #define ADD_HISTOGRAM_FILE(OBJNAME, FILENAME, XSEC) \
        OBJNAME.readFile(path+FILENAME, obsname, *definitions, XSEC)

    HISTOGRAM_FILE(ttZ_2016, "ttZ1_2016.root", xsec_ttZ*0.3334);
    // ADD_HISTOGRAM_FILE(ttZ_2016, "ttZ2_2016.root", xsec_ttZ);
    ADD_HISTOGRAM_FILE(ttZ_2016, "ttZ3_2016.root", xsec_ttZ*0.3333);
    ADD_HISTOGRAM_FILE(ttZ_2016, "ttZ4_2016.root", xsec_ttZ*0.3333);
    // ADD_HISTOGRAM_FILE(ttZ_2016, "ttZ5_2016.root", xsec_ttZ);
    HISTOGRAM_FILE(ttZ_2017, "ttZ1_2017.root", xsec_ttZ*0.2);
    ADD_HISTOGRAM_FILE(ttZ_2017, "ttZ2_2017.root", xsec_ttZ*0.2);
    ADD_HISTOGRAM_FILE(ttZ_2017, "ttZ3_2017.root", xsec_ttZ*0.2);
    ADD_HISTOGRAM_FILE(ttZ_2017, "ttZ4_2017.root", xsec_ttZ*0.2);
    ADD_HISTOGRAM_FILE(ttZ_2017, "ttZ5_2017.root", xsec_ttZ*0.2);
    HISTOGRAM_FILE(ttZ_2018, "ttZ1_2018.root", xsec_ttZ*0.2);
    ADD_HISTOGRAM_FILE(ttZ_2018, "ttZ2_2018.root", xsec_ttZ*0.2);
    ADD_HISTOGRAM_FILE(ttZ_2018, "ttZ3_2018.root", xsec_ttZ*0.2);
    ADD_HISTOGRAM_FILE(ttZ_2018, "ttZ4_2018.root", xsec_ttZ*0.2);
    ADD_HISTOGRAM_FILE(ttZ_2018, "ttZ5_2018.root", xsec_ttZ*0.2);
    plot.SetBins(ttZ_2016.binCenters, ttZ_2016.binBoundaries);
    for(int iDef=0; iDef<nDefinitions; iDef++) {
        plot.AddSample(ttZ_2016.values_gen[iDef], ttZ_2016.values_fail[iDef], ttZ_2016.values_response[iDef], 0, iDef);
        plot.AddSample(ttZ_2017.values_gen[iDef], ttZ_2017.values_fail[iDef], ttZ_2017.values_response[iDef], 1, iDef);
        plot.AddSample(ttZ_2018.values_gen[iDef], ttZ_2018.values_fail[iDef], ttZ_2018.values_response[iDef], 2, iDef);
    }

    plot.Evaluate();

    for(int iYear=0; iYear<=nYears; iYear++) {
        const std::string sYear = (iYear==0) ? "2016" : (iYear==1) ? "2017" : (iYear==2) ? "2018" : "run2";
        const std::string sLumi = (iYear==0) ? "2016" : (iYear==1) ? "2017" : (iYear==2) ? "2018" : "Run 2";
        for(int iDef=0; iDef<nDefinitions; iDef++) {
            const std::string sDef = sDefinitions->at(iDef);
            const int iTo = sDef.find("to");
            const std::string sTitle = sDef.substr(0, iTo)+" → "+sDef.substr(iTo+2);
            const std::string filenamepart = "_"+sYear+"_"+sDef+"_"+obsname;
            const std::string data_response = plot.PrintResponse(iYear, iDef);
            const std::string fullpath_response = outputdir+"/response"+filenamepart+".plot";
            std::cout << "Writing " << fullpath_response << std::endl;
            std::ofstream out_response(fullpath_response);
            out_response << data_response << std::endl
                         << "filename = " << '"' << "response" << filenamepart << '"' << std::endl
                         << "my_reclabel = " << '"' << "reconstructed " << get_xaxis_label(obsname) << '"' << std::endl
                         << "my_genlabel = " << '"' << "true " << get_xaxis_label(obsname) << '"' << std::endl
                         << "lumi = " << '"' << sLumi << '"' << std::endl
                         << "title = " << '"' << sTitle << '"' << std::endl
                         << "load " << '"' << "responsematrix.plot" << '"' << std::endl;
            out_response.close();
            const std::string data_stabpureff = plot.PrintStabPurEff(iYear, iDef);
            const std::string fullpath_stabpureff = outputdir+"/stabpureff"+filenamepart+".plot";
            std::cout << "Writing " << fullpath_stabpureff << std::endl;
            std::ofstream out_stabpureff(fullpath_stabpureff);
            out_stabpureff << data_stabpureff << std::endl
                         << "filename = " << '"' << "stabpureff" << filenamepart << '"' << std::endl
                         << "my_xlabel = " << '"' << get_xaxis_label(obsname) << '"' << std::endl
                         << "lumi = " << '"' << sLumi << '"' << std::endl
                         << "title = " << '"' << sTitle << '"' << std::endl
                         << "load " << '"' << "stabpureff.plot" << '"' << std::endl;
            out_stabpureff.close();
        }
    }
}

int main(int argc, char* argv[]){
    std::vector<std::string> argvStr(&argv[0], &argv[0]+argc);
    makeResponseMatrices(argvStr.at(1), stoi(argvStr.at(2)), stoi(argvStr.at(3)), argvStr.at(4));
}
