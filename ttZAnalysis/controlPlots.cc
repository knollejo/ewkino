#include <fstream>
#include <iostream>
#include <string>
#include <utility>
#include <vector>

#include <TFile.h>

#include "interface/ControlPlot.h"
#include "interface/CrossSections.h"
#include "interface/HistogramFile.h"
#include "interface/PlotAxes.h"

const int nYears = 3;
const int nSamples = 7;
const int nVarData = 2;
const int nVarBkgr = 3;
const int nVarSgnl = 18;
const int nSelections = 3;
const std::vector<std::pair<int, int>> selectionMap({
    std::make_pair(0, 0), // 4l -> 4l
    std::make_pair(1, 1), // 3l3j -> 3l
    std::make_pair(2, 1), // 3l4j -> 3l
    std::make_pair(0, 2), // 4l -> all
    std::make_pair(1, 2), // 3l3j -> all
    std::make_pair(2, 2), // 3l4j -> all
});
const int nSysData = 0;
const int nSysBkgr = 126;
const int nSysSgnl = 126;

const int i_data = 0;
const int i_ttZ3l = 1;
const int i_ttZ4l = 2;
const int i_ttZ2l = 3;
const int i_ttZtaus = 4;
const int i_ttZoffshell = 5;
const int i_ttZneutrinos = 6;
const int i_tZq = 6;
const int i_tWZ = 6;
const int i_ttX = 6;
const int i_VV = 6;
const int i_nonprompt = 6;

const std::string path = "output/210308_142313/";

void makeControlPlots(std::string obsname, int nBins, std::string outputdir) {

    ControlPlot plot(nSelections, nSamples, nBins);

    #define HISTOGRAM_FILE(OBJNAME, FILENAME, NVAR, NSYS, XSEC) \
        HistogramFile OBJNAME(path+FILENAME, obsname, NVAR, nSelections, selectionMap, NSYS, nBins, XSEC)
    #define ADD_HISTOGRAM_FILE(OBJNAME, FILENAME, XSEC) \
        OBJNAME.readFile(path+FILENAME, obsname, selectionMap, XSEC)

    HISTOGRAM_FILE(data_2016, "data_2016.root", nVarData, nSysData, 1.0);
    HISTOGRAM_FILE(data_2017, "data_2017.root", nVarData, nSysData, 1.0);
    HISTOGRAM_FILE(data_2018, "data_2018.root", nVarData, nSysData, 1.0);
    plot.SetBins(data_2016.binCenters, data_2017.binBoundaries);
    plot.AddSample(data_2016.values[0], 0, i_data, true);
    plot.AddSample(data_2017.values[0], 1, i_data, true);
    plot.AddSample(data_2018.values[0], 2, i_data, true);
    plot.AddSample(data_2016.values[1], 0, i_nonprompt, true);
    plot.AddSample(data_2017.values[1], 1, i_nonprompt, true);
    plot.AddSample(data_2018.values[1], 2, i_nonprompt, true);

    HISTOGRAM_FILE(ttZ_2016, "ttZ1_2016.root", nVarSgnl, nSysSgnl, xsec_ttZ*0.3334);
    // ADD_HISTOGRAM_FILE(ttZ_2016, "ttZ2_2016.root", xsec_ttZ);
    ADD_HISTOGRAM_FILE(ttZ_2016, "ttZ3_2016.root", xsec_ttZ*0.3333);
    ADD_HISTOGRAM_FILE(ttZ_2016, "ttZ4_2016.root", xsec_ttZ*0.3333);
    // ADD_HISTOGRAM_FILE(ttZ_2016, "ttZ5_2016.root", xsec_ttZ);
    HISTOGRAM_FILE(ttZ_2017, "ttZ1_2017.root", nVarSgnl, nSysSgnl, xsec_ttZ*0.2);
    ADD_HISTOGRAM_FILE(ttZ_2017, "ttZ2_2017.root", xsec_ttZ*0.2);
    ADD_HISTOGRAM_FILE(ttZ_2017, "ttZ3_2017.root", xsec_ttZ*0.2);
    ADD_HISTOGRAM_FILE(ttZ_2017, "ttZ4_2017.root", xsec_ttZ*0.2);
    ADD_HISTOGRAM_FILE(ttZ_2017, "ttZ5_2017.root", xsec_ttZ*0.2);
    HISTOGRAM_FILE(ttZ_2018, "ttZ1_2018.root", nVarSgnl, nSysSgnl, xsec_ttZ*0.2);
    ADD_HISTOGRAM_FILE(ttZ_2018, "ttZ2_2018.root", xsec_ttZ*0.2);
    ADD_HISTOGRAM_FILE(ttZ_2018, "ttZ3_2018.root", xsec_ttZ*0.2);
    ADD_HISTOGRAM_FILE(ttZ_2018, "ttZ4_2018.root", xsec_ttZ*0.2);
    ADD_HISTOGRAM_FILE(ttZ_2018, "ttZ5_2018.root", xsec_ttZ*0.2);
    plot.AddSample(ttZ_2016.values[0], 0, i_ttZ4l, false);
    plot.AddSample(ttZ_2017.values[0], 1, i_ttZ4l, false);
    plot.AddSample(ttZ_2018.values[0], 2, i_ttZ4l, false);
    plot.AddSample(ttZ_2016.values[3], 0, i_ttZ3l, false);
    plot.AddSample(ttZ_2017.values[3], 1, i_ttZ3l, false);
    plot.AddSample(ttZ_2018.values[3], 2, i_ttZ3l, false);
    plot.AddSample(ttZ_2016.values[6], 0, i_ttZ2l, false);
    plot.AddSample(ttZ_2017.values[6], 1, i_ttZ2l, false);
    plot.AddSample(ttZ_2018.values[6], 2, i_ttZ2l, false);
    plot.AddSample(ttZ_2016.values[9], 0, i_ttZtaus, false);
    plot.AddSample(ttZ_2017.values[9], 1, i_ttZtaus, false);
    plot.AddSample(ttZ_2018.values[9], 2, i_ttZtaus, false);
    plot.AddSample(ttZ_2016.values[12], 0, i_ttZoffshell, false);
    plot.AddSample(ttZ_2017.values[12], 1, i_ttZoffshell, false);
    plot.AddSample(ttZ_2018.values[12], 2, i_ttZoffshell, false);
    plot.AddSample(ttZ_2016.values[15], 0, i_ttZneutrinos, false);
    plot.AddSample(ttZ_2017.values[15], 1, i_ttZneutrinos, false);
    plot.AddSample(ttZ_2018.values[15], 2, i_ttZneutrinos, false);
    for(int iCat=0; iCat<6; iCat++) {
        plot.AddSample(ttZ_2016.values[iCat*3+2], 0, i_nonprompt, false, -1.0);
        plot.AddSample(ttZ_2017.values[iCat*3+2], 1, i_nonprompt, false, -1.0);
        plot.AddSample(ttZ_2018.values[iCat*3+2], 2, i_nonprompt, false, -1.0);
    }

    HISTOGRAM_FILE(ttH_2016, "ttH_2016.root", nVarBkgr, nSysBkgr, xsec_ttH);
    // HISTOGRAM_FILE(ttH_2017, "ttH_2017.root", nVarBkgr, nSysBkgr, xsec_ttH);
    HISTOGRAM_FILE(ttH_2018, "ttH_2018.root", nVarBkgr, nSysBkgr, xsec_ttH);
    plot.AddSample(ttH_2016.values[0], 0, i_ttX, false);
    // plot.AddSample(ttH_2017.values[0], 1, i_ttX, false);
    plot.AddSample(ttH_2018.values[0], 2, i_ttX, false);
    plot.AddSample(ttH_2016.values[2], 0, i_nonprompt, false, -1.0);
    // plot.AddSample(ttH_2017.values[2], 1, i_nonprompt, false, -1.0);
    plot.AddSample(ttH_2018.values[2], 2, i_nonprompt, false, -1.0);

    HISTOGRAM_FILE(tZq_2016, "tZq_2016.root", nVarBkgr, nSysBkgr, xsec_tZq);
    HISTOGRAM_FILE(tZq_2017, "tZq_2017.root", nVarBkgr, nSysBkgr, xsec_tZq);
    HISTOGRAM_FILE(tZq_2018, "tZq_2018.root", nVarBkgr, nSysBkgr, xsec_tZq);
    plot.AddSample(tZq_2016.values[0], 0, i_tZq, false);
    plot.AddSample(tZq_2017.values[0], 1, i_tZq, false);
    plot.AddSample(tZq_2018.values[0], 2, i_tZq, false);
    plot.AddSample(tZq_2016.values[2], 0, i_nonprompt, false, -1.0);
    plot.AddSample(tZq_2017.values[2], 1, i_nonprompt, false, -1.0);
    plot.AddSample(tZq_2018.values[2], 2, i_nonprompt, false, -1.0);

    HISTOGRAM_FILE(tWZ_2016, "tWZ_2016.root", nVarBkgr, nSysBkgr, xsec_tWZ);
    HISTOGRAM_FILE(tWZ_2017, "tWZ_2017.root", nVarBkgr, nSysBkgr, xsec_tWZ);
    HISTOGRAM_FILE(tWZ_2018, "tWZ_2018.root", nVarBkgr, nSysBkgr, xsec_tWZ);
    plot.AddSample(tWZ_2016.values[0], 0, i_tWZ, false);
    plot.AddSample(tWZ_2017.values[0], 1, i_tWZ, false);
    plot.AddSample(tWZ_2018.values[0], 2, i_tWZ, false);
    plot.AddSample(tWZ_2016.values[2], 0, i_nonprompt, false, -1.0);
    plot.AddSample(tWZ_2017.values[2], 1, i_nonprompt, false, -1.0);
    plot.AddSample(tWZ_2018.values[2], 2, i_nonprompt, false, -1.0);

    HISTOGRAM_FILE(ttW_2016, "ttW_2016.root", nVarBkgr, nSysBkgr, xsec_ttW);
    HISTOGRAM_FILE(ttW_2017, "ttW_2017.root", nVarBkgr, nSysBkgr, xsec_ttW);
    HISTOGRAM_FILE(ttW_2018, "ttW_2018.root", nVarBkgr, nSysBkgr, xsec_ttW);
    plot.AddSample(ttW_2016.values[0], 0, i_ttX, false);
    plot.AddSample(ttW_2017.values[0], 1, i_ttX, false);
    plot.AddSample(ttW_2018.values[0], 2, i_ttX, false);
    plot.AddSample(ttW_2016.values[2], 0, i_nonprompt, false, -1.0);
    plot.AddSample(ttW_2017.values[2], 1, i_nonprompt, false, -1.0);
    plot.AddSample(ttW_2018.values[2], 2, i_nonprompt, false, -1.0);

    HISTOGRAM_FILE(tHQ_2016, "tHQ_2016.root", nVarBkgr, nSysBkgr, xsec_tHQ);
    HISTOGRAM_FILE(tHQ_2017, "tHQ_2017.root", nVarBkgr, nSysBkgr, xsec_tHQ);
    HISTOGRAM_FILE(tHQ_2018, "tHQ_2018.root", nVarBkgr, nSysBkgr, xsec_tHQ);
    plot.AddSample(tHQ_2016.values[0], 0, i_ttX, false);
    plot.AddSample(tHQ_2017.values[0], 1, i_ttX, false);
    plot.AddSample(tHQ_2018.values[0], 2, i_ttX, false);
    plot.AddSample(tHQ_2016.values[2], 0, i_nonprompt, false, -1.0);
    plot.AddSample(tHQ_2017.values[2], 1, i_nonprompt, false, -1.0);
    plot.AddSample(tHQ_2018.values[2], 2, i_nonprompt, false, -1.0);

    HISTOGRAM_FILE(tHW_2016, "tHW_2016.root", nVarBkgr, nSysBkgr, xsec_tHW);
    HISTOGRAM_FILE(tHW_2017, "tHW_2017.root", nVarBkgr, nSysBkgr, xsec_tHW);
    HISTOGRAM_FILE(tHW_2018, "tHW_2018.root", nVarBkgr, nSysBkgr, xsec_tHW);
    plot.AddSample(tHW_2016.values[0], 0, i_ttX, false);
    plot.AddSample(tHW_2017.values[0], 1, i_ttX, false);
    plot.AddSample(tHW_2018.values[0], 2, i_ttX, false);
    plot.AddSample(tHW_2016.values[2], 0, i_nonprompt, false, -1.0);
    plot.AddSample(tHW_2017.values[2], 1, i_nonprompt, false, -1.0);
    plot.AddSample(tHW_2018.values[2], 2, i_nonprompt, false, -1.0);

    HISTOGRAM_FILE(ttZlight_2016, "ttZlight_2016.root", nVarBkgr, nSysBkgr, xsec_ttZlight);
    HISTOGRAM_FILE(ttZlight_2017, "ttZlight_2017.root", nVarBkgr, nSysBkgr, xsec_ttZlight);
    HISTOGRAM_FILE(ttZlight_2018, "ttZlight_2018.root", nVarBkgr, nSysBkgr, xsec_ttZlight);
    plot.AddSample(ttZlight_2016.values[0], 0, i_ttX, false);
    plot.AddSample(ttZlight_2017.values[0], 1, i_ttX, false);
    plot.AddSample(ttZlight_2018.values[0], 2, i_ttX, false);
    plot.AddSample(ttZlight_2016.values[2], 0, i_nonprompt, false, -1.0);
    plot.AddSample(ttZlight_2017.values[2], 1, i_nonprompt, false, -1.0);
    plot.AddSample(ttZlight_2018.values[2], 2, i_nonprompt, false, -1.0);

    HISTOGRAM_FILE(tttt_2016, "tttt_2016.root", nVarBkgr, nSysBkgr, xsec_tttt);
    HISTOGRAM_FILE(tttt_2017, "tttt_2017.root", nVarBkgr, nSysBkgr, xsec_tttt);
    HISTOGRAM_FILE(tttt_2018, "tttt_2018.root", nVarBkgr, nSysBkgr, xsec_tttt);
    plot.AddSample(tttt_2016.values[0], 0, i_ttX, false);
    plot.AddSample(tttt_2017.values[0], 1, i_ttX, false);
    plot.AddSample(tttt_2018.values[0], 2, i_ttX, false);
    plot.AddSample(tttt_2016.values[2], 0, i_nonprompt, false, -1.0);
    plot.AddSample(tttt_2017.values[2], 1, i_nonprompt, false, -1.0);
    plot.AddSample(tttt_2018.values[2], 2, i_nonprompt, false, -1.0);

    HISTOGRAM_FILE(ttWW_2016, "ttWW_2016.root", nVarBkgr, nSysBkgr, xsec_ttWW);
    HISTOGRAM_FILE(ttWW_2017, "ttWW_2017.root", nVarBkgr, nSysBkgr, xsec_ttWW);
    HISTOGRAM_FILE(ttWW_2018, "ttWW_2018.root", nVarBkgr, nSysBkgr, xsec_ttWW);
    plot.AddSample(ttWW_2016.values[0], 0, i_ttX, false);
    plot.AddSample(ttWW_2017.values[0], 1, i_ttX, false);
    plot.AddSample(ttWW_2018.values[0], 2, i_ttX, false);
    plot.AddSample(ttWW_2016.values[2], 0, i_nonprompt, false, -1.0);
    plot.AddSample(ttWW_2017.values[2], 1, i_nonprompt, false, -1.0);
    plot.AddSample(ttWW_2018.values[2], 2, i_nonprompt, false, -1.0);

    HISTOGRAM_FILE(ttWZ_2016, "ttWZ_2016.root", nVarBkgr, nSysBkgr, xsec_ttWZ);
    HISTOGRAM_FILE(ttWZ_2017, "ttWZ_2017.root", nVarBkgr, nSysBkgr, xsec_ttWZ);
    HISTOGRAM_FILE(ttWZ_2018, "ttWZ_2018.root", nVarBkgr, nSysBkgr, xsec_ttWZ);
    plot.AddSample(ttWZ_2016.values[0], 0, i_ttX, false);
    plot.AddSample(ttWZ_2017.values[0], 1, i_ttX, false);
    plot.AddSample(ttWZ_2018.values[0], 2, i_ttX, false);
    plot.AddSample(ttWZ_2016.values[2], 0, i_nonprompt, false, -1.0);
    plot.AddSample(ttWZ_2017.values[2], 1, i_nonprompt, false, -1.0);
    plot.AddSample(ttWZ_2018.values[2], 2, i_nonprompt, false, -1.0);

    HISTOGRAM_FILE(ttZZ_2016, "ttZZ_2016.root", nVarBkgr, nSysBkgr, xsec_ttZZ);
    HISTOGRAM_FILE(ttZZ_2017, "ttZZ_2017.root", nVarBkgr, nSysBkgr, xsec_ttZZ);
    HISTOGRAM_FILE(ttZZ_2018, "ttZZ_2018.root", nVarBkgr, nSysBkgr, xsec_ttZZ);
    plot.AddSample(ttZZ_2016.values[0], 0, i_ttX, false);
    plot.AddSample(ttZZ_2017.values[0], 1, i_ttX, false);
    plot.AddSample(ttZZ_2018.values[0], 2, i_ttX, false);
    plot.AddSample(ttZZ_2016.values[2], 0, i_nonprompt, false, -1.0);
    plot.AddSample(ttZZ_2017.values[2], 1, i_nonprompt, false, -1.0);
    plot.AddSample(ttZZ_2018.values[2], 2, i_nonprompt, false, -1.0);

    HISTOGRAM_FILE(WZ3L_2016, "WZ3L_2016.root", nVarBkgr, nSysBkgr, xsec_WZ3L);
    HISTOGRAM_FILE(WZ3L_2017, "WZ3L_2017.root", nVarBkgr, nSysBkgr, xsec_WZ3L);
    HISTOGRAM_FILE(WZ3L_2018, "WZ3L_2018.root", nVarBkgr, nSysBkgr, xsec_WZ3L);
    plot.AddSample(WZ3L_2016.values[0], 0, i_VV, false);
    plot.AddSample(WZ3L_2017.values[0], 1, i_VV, false);
    plot.AddSample(WZ3L_2018.values[0], 2, i_VV, false);
    plot.AddSample(WZ3L_2016.values[2], 0, i_nonprompt, false, -1.0);
    plot.AddSample(WZ3L_2017.values[2], 1, i_nonprompt, false, -1.0);
    plot.AddSample(WZ3L_2018.values[2], 2, i_nonprompt, false, -1.0);

    // HISTOGRAM_FILE(WZ2L_2016, "WZ2L_2016.root", nVarBkgr, nSysBkgr, xsec_WZ2L);
    // HISTOGRAM_FILE(WZ2L_2017, "WZ2L_2017.root", nVarBkgr, nSysBkgr, xsec_WZ2L);
    // HISTOGRAM_FILE(WZ2L_2018, "WZ2L_2018.root", nVarBkgr, nSysBkgr, xsec_WZ2L);
    // plot.AddSample(WZ2L_2016.values[0], 0, i_VV, false);
    // plot.AddSample(WZ2L_2017.values[0], 1, i_VV, false);
    // plot.AddSample(WZ2L_2018.values[0], 2, i_VV, false);
    // plot.AddSample(WZ2L_2016.values[2], 0, i_nonprompt, false, -1.0);
    // plot.AddSample(WZ2L_2017.values[2], 1, i_nonprompt, false, -1.0);
    // plot.AddSample(WZ2L_2018.values[2], 2, i_nonprompt, false, -1.0);

    HISTOGRAM_FILE(DY_2016, "DY_2016.root", nVarBkgr, nSysBkgr, xsec_DY);
    HISTOGRAM_FILE(DY_2017, "DY_2017.root", nVarBkgr, nSysBkgr, xsec_DY);
    HISTOGRAM_FILE(DY_2018, "DY_2018.root", nVarBkgr, nSysBkgr, xsec_DY);
    plot.AddSample(DY_2016.values[0], 0, i_VV, false);
    plot.AddSample(DY_2017.values[0], 1, i_VV, false);
    plot.AddSample(DY_2018.values[0], 2, i_VV, false);
    plot.AddSample(DY_2016.values[2], 0, i_nonprompt, false, -1.0);
    plot.AddSample(DY_2017.values[2], 1, i_nonprompt, false, -1.0);
    plot.AddSample(DY_2018.values[2], 2, i_nonprompt, false, -1.0);

    HISTOGRAM_FILE(tG_2016, "tG_2016.root", nVarBkgr, nSysBkgr, xsec_tG);
    HISTOGRAM_FILE(tG_2017, "tG_2017.root", nVarBkgr, nSysBkgr, xsec_tG);
    HISTOGRAM_FILE(tG_2018, "tG_2018.root", nVarBkgr, nSysBkgr, xsec_tG);
    plot.AddSample(tG_2016.values[0], 0, i_ttX, false);
    plot.AddSample(tG_2017.values[0], 1, i_ttX, false);
    plot.AddSample(tG_2018.values[0], 2, i_ttX, false);
    plot.AddSample(tG_2016.values[2], 0, i_nonprompt, false, -1.0);
    plot.AddSample(tG_2017.values[2], 1, i_nonprompt, false, -1.0);
    plot.AddSample(tG_2018.values[2], 2, i_nonprompt, false, -1.0);

    HISTOGRAM_FILE(ttG_2016, "ttG_2016.root", nVarBkgr, nSysBkgr, xsec_ttG);
    HISTOGRAM_FILE(ttG_2017, "ttG_2017.root", nVarBkgr, nSysBkgr, xsec_ttG);
    HISTOGRAM_FILE(ttG_2018, "ttG_2018.root", nVarBkgr, nSysBkgr, xsec_ttG);
    plot.AddSample(ttG_2016.values[0], 0, i_ttX, false);
    plot.AddSample(ttG_2017.values[0], 1, i_ttX, false);
    plot.AddSample(ttG_2018.values[0], 2, i_ttX, false);
    plot.AddSample(ttG_2016.values[2], 0, i_nonprompt, false, -1.0);
    plot.AddSample(ttG_2017.values[2], 1, i_nonprompt, false, -1.0);
    plot.AddSample(ttG_2018.values[2], 2, i_nonprompt, false, -1.0);

    HISTOGRAM_FILE(WG_2016, "WG_2016.root", nVarBkgr, nSysBkgr, xsec_WG);
    HISTOGRAM_FILE(WG_2017, "WG_2017.root", nVarBkgr, nSysBkgr, xsec_WG);
    HISTOGRAM_FILE(WG_2018, "WG_2018.root", nVarBkgr, nSysBkgr, xsec_WG);
    plot.AddSample(WG_2016.values[0], 0, i_VV, false);
    plot.AddSample(WG_2017.values[0], 1, i_VV, false);
    plot.AddSample(WG_2018.values[0], 2, i_VV, false);
    plot.AddSample(WG_2016.values[2], 0, i_nonprompt, false, -1.0);
    plot.AddSample(WG_2017.values[2], 1, i_nonprompt, false, -1.0);
    plot.AddSample(WG_2018.values[2], 2, i_nonprompt, false, -1.0);

    // HISTOGRAM_FILE(tt_2016, "tt_2016.root", nVarBkgr, nSysBkgr, xsec_tt);
    // HISTOGRAM_FILE(tt_2017, "tt_2017.root", nVarBkgr, nSysBkgr, xsec_tt);
    // HISTOGRAM_FILE(tt_2018, "tt_2018.root", nVarBkgr, nSysBkgr, xsec_tt);
    // plot.AddSample(tt_2016.values[0], 0, i_ttX, false);
    // plot.AddSample(tt_2017.values[0], 1, i_ttX, false);
    // plot.AddSample(tt_2018.values[0], 2, i_ttX, false);
    // plot.AddSample(tt_2016.values[2], 0, i_nonprompt, false, -1.0);
    // plot.AddSample(tt_2017.values[2], 1, i_nonprompt, false, -1.0);
    // plot.AddSample(tt_2018.values[2], 2, i_nonprompt, false, -1.0);

    HISTOGRAM_FILE(ZZ4L_2016, "ZZ4L_2016.root", nVarBkgr, nSysBkgr, xsec_ZZ4L);
    HISTOGRAM_FILE(ZZ4L_2017, "ZZ4L_2017.root", nVarBkgr, nSysBkgr, xsec_ZZ4L);
    HISTOGRAM_FILE(ZZ4L_2018, "ZZ4L_2018.root", nVarBkgr, nSysBkgr, xsec_ZZ4L);
    plot.AddSample(ZZ4L_2016.values[0], 0, i_VV, false);
    plot.AddSample(ZZ4L_2017.values[0], 1, i_VV, false);
    plot.AddSample(ZZ4L_2018.values[0], 2, i_VV, false);
    plot.AddSample(ZZ4L_2016.values[2], 0, i_nonprompt, false, -1.0);
    plot.AddSample(ZZ4L_2017.values[2], 1, i_nonprompt, false, -1.0);
    plot.AddSample(ZZ4L_2018.values[2], 2, i_nonprompt, false, -1.0);

    HISTOGRAM_FILE(ZZ2E2M_2016, "ZZ2E2M_2016.root", nVarBkgr, nSysBkgr, xsec_ZZ2E2M);
    HISTOGRAM_FILE(ZZ2E2M_2017, "ZZ2E2M_2017.root", nVarBkgr, nSysBkgr, xsec_ZZ2E2M);
    HISTOGRAM_FILE(ZZ2E2M_2018, "ZZ2E2M_2018.root", nVarBkgr, nSysBkgr, xsec_ZZ2E2M);
    plot.AddSample(ZZ2E2M_2016.values[0], 0, i_VV, false);
    plot.AddSample(ZZ2E2M_2017.values[0], 1, i_VV, false);
    plot.AddSample(ZZ2E2M_2018.values[0], 2, i_VV, false);
    plot.AddSample(ZZ2E2M_2016.values[2], 0, i_nonprompt, false, -1.0);
    plot.AddSample(ZZ2E2M_2017.values[2], 1, i_nonprompt, false, -1.0);
    plot.AddSample(ZZ2E2M_2018.values[2], 2, i_nonprompt, false, -1.0);

    HISTOGRAM_FILE(ZZ2E2T_2016, "ZZ2E2T_2016.root", nVarBkgr, nSysBkgr, xsec_ZZ2E2T);
    HISTOGRAM_FILE(ZZ2E2T_2017, "ZZ2E2T_2017.root", nVarBkgr, nSysBkgr, xsec_ZZ2E2T);
    HISTOGRAM_FILE(ZZ2E2T_2018, "ZZ2E2T_2018.root", nVarBkgr, nSysBkgr, xsec_ZZ2E2T);
    plot.AddSample(ZZ2E2T_2016.values[0], 0, i_VV, false);
    plot.AddSample(ZZ2E2T_2017.values[0], 1, i_VV, false);
    plot.AddSample(ZZ2E2T_2018.values[0], 2, i_VV, false);
    plot.AddSample(ZZ2E2T_2016.values[2], 0, i_nonprompt, false, -1.0);
    plot.AddSample(ZZ2E2T_2017.values[2], 1, i_nonprompt, false, -1.0);
    plot.AddSample(ZZ2E2T_2018.values[2], 2, i_nonprompt, false, -1.0);

    HISTOGRAM_FILE(ZZ2M2T_2016, "ZZ2M2T_2016.root", nVarBkgr, nSysBkgr, xsec_ZZ2M2T);
    HISTOGRAM_FILE(ZZ2M2T_2017, "ZZ2M2T_2017.root", nVarBkgr, nSysBkgr, xsec_ZZ2M2T);
    HISTOGRAM_FILE(ZZ2M2T_2018, "ZZ2M2T_2018.root", nVarBkgr, nSysBkgr, xsec_ZZ2M2T);
    plot.AddSample(ZZ2M2T_2016.values[0], 0, i_VV, false);
    plot.AddSample(ZZ2M2T_2017.values[0], 1, i_VV, false);
    plot.AddSample(ZZ2M2T_2018.values[0], 2, i_VV, false);
    plot.AddSample(ZZ2M2T_2016.values[2], 0, i_nonprompt, false, -1.0);
    plot.AddSample(ZZ2M2T_2017.values[2], 1, i_nonprompt, false, -1.0);
    plot.AddSample(ZZ2M2T_2018.values[2], 2, i_nonprompt, false, -1.0);

    HISTOGRAM_FILE(ZZ4E_2016, "ZZ4E_2016.root", nVarBkgr, nSysBkgr, xsec_ZZ4E);
    HISTOGRAM_FILE(ZZ4E_2017, "ZZ4E_2017.root", nVarBkgr, nSysBkgr, xsec_ZZ4E);
    HISTOGRAM_FILE(ZZ4E_2018, "ZZ4E_2018.root", nVarBkgr, nSysBkgr, xsec_ZZ4E);
    plot.AddSample(ZZ4E_2016.values[0], 0, i_VV, false);
    plot.AddSample(ZZ4E_2017.values[0], 1, i_VV, false);
    plot.AddSample(ZZ4E_2018.values[0], 2, i_VV, false);
    plot.AddSample(ZZ4E_2016.values[2], 0, i_nonprompt, false, -1.0);
    plot.AddSample(ZZ4E_2017.values[2], 1, i_nonprompt, false, -1.0);
    plot.AddSample(ZZ4E_2018.values[2], 2, i_nonprompt, false, -1.0);

    HISTOGRAM_FILE(ZZ4M_2016, "ZZ4M_2016.root", nVarBkgr, nSysBkgr, xsec_ZZ4M);
    HISTOGRAM_FILE(ZZ4M_2017, "ZZ4M_2017.root", nVarBkgr, nSysBkgr, xsec_ZZ4M);
    HISTOGRAM_FILE(ZZ4M_2018, "ZZ4M_2018.root", nVarBkgr, nSysBkgr, xsec_ZZ4M);
    plot.AddSample(ZZ4M_2016.values[0], 0, i_VV, false);
    plot.AddSample(ZZ4M_2017.values[0], 1, i_VV, false);
    plot.AddSample(ZZ4M_2018.values[0], 2, i_VV, false);
    plot.AddSample(ZZ4M_2016.values[2], 0, i_nonprompt, false, -1.0);
    plot.AddSample(ZZ4M_2017.values[2], 1, i_nonprompt, false, -1.0);
    plot.AddSample(ZZ4M_2018.values[2], 2, i_nonprompt, false, -1.0);

    HISTOGRAM_FILE(ggHZZ_2016, "ggHZZ_2016.root", nVarBkgr, nSysBkgr, xsec_ggHZZ);
    HISTOGRAM_FILE(ggHZZ_2017, "ggHZZ_2017.root", nVarBkgr, nSysBkgr, xsec_ggHZZ);
    HISTOGRAM_FILE(ggHZZ_2018, "ggHZZ_2018.root", nVarBkgr, nSysBkgr, xsec_ggHZZ);
    plot.AddSample(ggHZZ_2016.values[0], 0, i_VV, false);
    plot.AddSample(ggHZZ_2017.values[0], 1, i_VV, false);
    plot.AddSample(ggHZZ_2018.values[0], 2, i_VV, false);
    plot.AddSample(ggHZZ_2016.values[2], 0, i_nonprompt, false, -1.0);
    plot.AddSample(ggHZZ_2017.values[2], 1, i_nonprompt, false, -1.0);
    plot.AddSample(ggHZZ_2018.values[2], 2, i_nonprompt, false, -1.0);

    HISTOGRAM_FILE(VBFHZZ_2016, "VBFHZZ_2016.root", nVarBkgr, nSysBkgr, xsec_VBFHZZ);
    HISTOGRAM_FILE(VBFHZZ_2017, "VBFHZZ_2017.root", nVarBkgr, nSysBkgr, xsec_VBFHZZ);
    HISTOGRAM_FILE(VBFHZZ_2018, "VBFHZZ_2018.root", nVarBkgr, nSysBkgr, xsec_VBFHZZ);
    plot.AddSample(VBFHZZ_2016.values[0], 0, i_VV, false);
    plot.AddSample(VBFHZZ_2017.values[0], 1, i_VV, false);
    plot.AddSample(VBFHZZ_2018.values[0], 2, i_VV, false);
    plot.AddSample(VBFHZZ_2016.values[2], 0, i_nonprompt, false, -1.0);
    plot.AddSample(VBFHZZ_2017.values[2], 1, i_nonprompt, false, -1.0);
    plot.AddSample(VBFHZZ_2018.values[2], 2, i_nonprompt, false, -1.0);

    // HISTOGRAM_FILE(WpHZZ_2016, "WpHZZ_2016.root", nVarBkgr, nSysBkgr, xsec_WpHZZ);
    // HISTOGRAM_FILE(WpHZZ_2017, "WpHZZ_2017.root", nVarBkgr, nSysBkgr, xsec_WpHZZ);
    // HISTOGRAM_FILE(WpHZZ_2018, "WpHZZ_2018.root", nVarBkgr, nSysBkgr, xsec_WpHZZ);
    // plot.AddSample(WpHZZ_2016.values[0], 0, i_VV, false);
    // plot.AddSample(WpHZZ_2017.values[0], 1, i_VV, false);
    // plot.AddSample(WpHZZ_2018.values[0], 2, i_VV, false);
    // plot.AddSample(WpHZZ_2016.values[2], 0, i_nonprompt, false, -1.0);
    // plot.AddSample(WpHZZ_2017.values[2], 1, i_nonprompt, false, -1.0);
    // plot.AddSample(WpHZZ_2018.values[2], 2, i_nonprompt, false, -1.0);

    HISTOGRAM_FILE(WmHZZ_2016, "WmHZZ_2016.root", nVarBkgr, nSysBkgr, xsec_WmHZZ);
    HISTOGRAM_FILE(WmHZZ_2017, "WmHZZ_2017.root", nVarBkgr, nSysBkgr, xsec_WmHZZ);
    HISTOGRAM_FILE(WmHZZ_2018, "WmHZZ_2018.root", nVarBkgr, nSysBkgr, xsec_WmHZZ);
    plot.AddSample(WmHZZ_2016.values[0], 0, i_VV, false);
    plot.AddSample(WmHZZ_2017.values[0], 1, i_VV, false);
    plot.AddSample(WmHZZ_2018.values[0], 2, i_VV, false);
    plot.AddSample(WmHZZ_2016.values[2], 0, i_nonprompt, false, -1.0);
    plot.AddSample(WmHZZ_2017.values[2], 1, i_nonprompt, false, -1.0);
    plot.AddSample(WmHZZ_2018.values[2], 2, i_nonprompt, false, -1.0);

    HISTOGRAM_FILE(ZHZZ_2016, "ZHZZ_2016.root", nVarBkgr, nSysBkgr, xsec_ZHZZ);
    HISTOGRAM_FILE(ZHZZ_2017, "ZHZZ_2017.root", nVarBkgr, nSysBkgr, xsec_ZHZZ);
    HISTOGRAM_FILE(ZHZZ_2018, "ZHZZ_2018.root", nVarBkgr, nSysBkgr, xsec_ZHZZ);
    plot.AddSample(ZHZZ_2016.values[0], 0, i_VV, false);
    plot.AddSample(ZHZZ_2017.values[0], 1, i_VV, false);
    plot.AddSample(ZHZZ_2018.values[0], 2, i_VV, false);
    plot.AddSample(ZHZZ_2016.values[2], 0, i_nonprompt, false, -1.0);
    plot.AddSample(ZHZZ_2017.values[2], 1, i_nonprompt, false, -1.0);
    plot.AddSample(ZHZZ_2018.values[2], 2, i_nonprompt, false, -1.0);

    HISTOGRAM_FILE(WZG_2016, "WZG_2016.root", nVarBkgr, nSysBkgr, xsec_WZG);
    HISTOGRAM_FILE(WZG_2017, "WZG_2017.root", nVarBkgr, nSysBkgr, xsec_WZG);
    HISTOGRAM_FILE(WZG_2018, "WZG_2018.root", nVarBkgr, nSysBkgr, xsec_WZG);
    plot.AddSample(WZG_2016.values[0], 0, i_VV, false);
    plot.AddSample(WZG_2017.values[0], 1, i_VV, false);
    plot.AddSample(WZG_2018.values[0], 2, i_VV, false);
    plot.AddSample(WZG_2016.values[2], 0, i_nonprompt, false, -1.0);
    plot.AddSample(WZG_2017.values[2], 1, i_nonprompt, false, -1.0);
    plot.AddSample(WZG_2018.values[2], 2, i_nonprompt, false, -1.0);

    HISTOGRAM_FILE(ZZZ_2016, "ZZZ_2016.root", nVarBkgr, nSysBkgr, xsec_ZZZ);
    HISTOGRAM_FILE(ZZZ_2017, "ZZZ_2017.root", nVarBkgr, nSysBkgr, xsec_ZZZ);
    HISTOGRAM_FILE(ZZZ_2018, "ZZZ_2018.root", nVarBkgr, nSysBkgr, xsec_ZZZ);
    plot.AddSample(ZZZ_2016.values[0], 0, i_VV, false);
    plot.AddSample(ZZZ_2017.values[0], 1, i_VV, false);
    plot.AddSample(ZZZ_2018.values[0], 2, i_VV, false);
    plot.AddSample(ZZZ_2016.values[2], 0, i_nonprompt, false, -1.0);
    plot.AddSample(ZZZ_2017.values[2], 1, i_nonprompt, false, -1.0);
    plot.AddSample(ZZZ_2018.values[2], 2, i_nonprompt, false, -1.0);

    HISTOGRAM_FILE(WZZ_2016, "WZZ_2016.root", nVarBkgr, nSysBkgr, xsec_WZZ);
    HISTOGRAM_FILE(WZZ_2017, "WZZ_2017.root", nVarBkgr, nSysBkgr, xsec_WZZ);
    HISTOGRAM_FILE(WZZ_2018, "WZZ_2018.root", nVarBkgr, nSysBkgr, xsec_WZZ);
    plot.AddSample(WZZ_2016.values[0], 0, i_VV, false);
    plot.AddSample(WZZ_2017.values[0], 1, i_VV, false);
    plot.AddSample(WZZ_2018.values[0], 2, i_VV, false);
    plot.AddSample(WZZ_2016.values[2], 0, i_nonprompt, false, -1.0);
    plot.AddSample(WZZ_2017.values[2], 1, i_nonprompt, false, -1.0);
    plot.AddSample(WZZ_2018.values[2], 2, i_nonprompt, false, -1.0);

    HISTOGRAM_FILE(WWZ_2016, "WWZ_2016.root", nVarBkgr, nSysBkgr, xsec_WWZ);
    HISTOGRAM_FILE(WWZ_2017, "WWZ_2017.root", nVarBkgr, nSysBkgr, xsec_WWZ);
    HISTOGRAM_FILE(WWZ_2018, "WWZ_2018.root", nVarBkgr, nSysBkgr, xsec_WWZ);
    plot.AddSample(WWZ_2016.values[0], 0, i_VV, false);
    plot.AddSample(WWZ_2017.values[0], 1, i_VV, false);
    plot.AddSample(WWZ_2018.values[0], 2, i_VV, false);
    plot.AddSample(WWZ_2016.values[2], 0, i_nonprompt, false, -1.0);
    plot.AddSample(WWZ_2017.values[2], 1, i_nonprompt, false, -1.0);
    plot.AddSample(WWZ_2018.values[2], 2, i_nonprompt, false, -1.0);

    HISTOGRAM_FILE(WWW_2016, "WWW_2016.root", nVarBkgr, nSysBkgr, xsec_WWW);
    HISTOGRAM_FILE(WWW_2017, "WWW_2017.root", nVarBkgr, nSysBkgr, xsec_WWW);
    HISTOGRAM_FILE(WWW_2018, "WWW_2018.root", nVarBkgr, nSysBkgr, xsec_WWW);
    plot.AddSample(WWW_2016.values[0], 0, i_VV, false);
    plot.AddSample(WWW_2017.values[0], 1, i_VV, false);
    plot.AddSample(WWW_2018.values[0], 2, i_VV, false);
    plot.AddSample(WWW_2016.values[2], 0, i_nonprompt, false, -1.0);
    plot.AddSample(WWW_2017.values[2], 1, i_nonprompt, false, -1.0);
    plot.AddSample(WWW_2018.values[2], 2, i_nonprompt, false, -1.0);

    HISTOGRAM_FILE(WWDS_2016, "WWDS_2016.root", nVarBkgr, nSysBkgr, xsec_WWDS);
    HISTOGRAM_FILE(WWDS_2017, "WWDS_2017.root", nVarBkgr, nSysBkgr, xsec_WWDS);
    HISTOGRAM_FILE(WWDS_2018, "WWDS_2018.root", nVarBkgr, nSysBkgr, xsec_WWDS);
    plot.AddSample(WWDS_2016.values[0], 0, i_VV, false);
    plot.AddSample(WWDS_2017.values[0], 1, i_VV, false);
    plot.AddSample(WWDS_2018.values[0], 2, i_VV, false);
    plot.AddSample(WWDS_2016.values[2], 0, i_nonprompt, false, -1.0);
    plot.AddSample(WWDS_2017.values[2], 1, i_nonprompt, false, -1.0);
    plot.AddSample(WWDS_2018.values[2], 2, i_nonprompt, false, -1.0);

    HISTOGRAM_FILE(WW_2016, "WW_2016.root", nVarBkgr, nSysBkgr, xsec_WW);
    HISTOGRAM_FILE(WW_2017, "WW_2017.root", nVarBkgr, nSysBkgr, xsec_WW);
    HISTOGRAM_FILE(WW_2018, "WW_2018.root", nVarBkgr, nSysBkgr, xsec_WW);
    plot.AddSample(WW_2016.values[0], 0, i_VV, false);
    plot.AddSample(WW_2017.values[0], 1, i_VV, false);
    plot.AddSample(WW_2018.values[0], 2, i_VV, false);
    plot.AddSample(WW_2016.values[2], 0, i_nonprompt, false, -1.0);
    plot.AddSample(WW_2017.values[2], 1, i_nonprompt, false, -1.0);
    plot.AddSample(WW_2018.values[2], 2, i_nonprompt, false, -1.0);

    plot.Evaluate();

    for(int iYear=0; iYear<=nYears; iYear++) {
        const std::string sYear = (iYear==0) ? "2016" : (iYear==1) ? "2017" : (iYear==2) ? "2018" : "run2";
        const std::string sLumi = (iYear==0) ? "35.9" : (iYear==1) ? "41.5" : (iYear==2) ? "59.7" : "137";
        for(int iSel=0; iSel<nSelections; iSel++) {
            const std::string sLep = (iSel==0) ? "=4ℓ" : (iSel==1) ? "=3ℓ" : "3ℓ/4ℓ";
            const std::string sJet = (iSel==0) ? "≥2j" : (iSel==1) ? "≥3j" : "≥3j/2j";
            const std::string sTitle = sLep+", "+sJet+", ≥1b";
            const std::string sSelection = "sel"+std::to_string(iSel);
            const std::string data = plot.Print(iYear, iSel);
            const std::string filename = "plot_"+sYear+"_"+sSelection+"_"+obsname;
            const std::string fullpath = outputdir+"/"+filename+".plot";
            std::cout << "Writing " << fullpath << std::endl;
            std::ofstream out(fullpath);
            out << data << std::endl
                << "filename = " << '"' << filename << '"' << std::endl
                << "my_xlabel = " << '"' << get_xaxis_label(obsname) << '"' << std::endl
                << "add_first_plot = " << '"' << get_xaxis_add_first(obsname) << '"' << std::endl
                << "add_second_plot = " << '"' << get_xaxis_add_second(obsname) << '"' << std::endl
                << "lumi = " << '"' << sLumi << '"' << std::endl
                << "title = " << '"' << sTitle << '"' << std::endl
                << "load " << '"' << "controlplotswithratio.plot" << '"' << std::endl;
            out.close();
        }
    }

}

int main(int argc, char* argv[]){
    std::vector<std::string> argvStr(&argv[0], &argv[0]+argc);
    makeControlPlots(argvStr.at(1), stoi(argvStr.at(2)), argvStr.at(3));
}
