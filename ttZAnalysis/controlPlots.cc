#include <fstream>
#include <iostream>
#include <string>

#include <TFile.h>

#include "interface/ControlPlot.h"
#include "interface/CrossSections.h"
#include "interface/HistogramFile.h"

const int nYears = 3;
const int nSamples = 7;
const int nVarData = 2;
const int nVarBkgr = 3;
const int nVarSgnl = 3;
const int nSelections = 3;
const int nSysData = 0;
const int nSysBkgr = 126;
const int nSysSgnl = 126;

const std::string path = "output/210222_102836/";

void makeControlPlots(std::string obsname, int nBins, std::string outputdir) {

    ControlPlot plot(nSelections, nSamples, nBins);

    HistogramFile data_2016(path+"data_2016.root", obsname, nVarData, nSelections, nSysData, nBins, 1.0);
    HistogramFile data_2017(path+"data_2017.root", obsname, nVarData, nSelections, nSysData, nBins, 1.0);
    HistogramFile data_2018(path+"data_2018.root", obsname, nVarData, nSelections, nSysData, nBins, 1.0);
    plot.SetBins(data_2016.binCenters, data_2017.binBoundaries);
    plot.AddSample(data_2016.values[0], 0, 0, true);
    plot.AddSample(data_2017.values[0], 1, 0, true);
    plot.AddSample(data_2018.values[0], 2, 0, true);
    plot.AddSample(data_2016.values[1], 0, 6, true);
    plot.AddSample(data_2017.values[1], 1, 6, true);
    plot.AddSample(data_2018.values[1], 2, 6, true);

    HistogramFile ttZ_2016(path+"ttZ_2016.root", obsname, nVarSgnl, nSelections, nSysSgnl, nBins, xsec_ttZ);
    HistogramFile ttZ_2017(path+"ttZ_2017.root", obsname, nVarSgnl, nSelections, nSysSgnl, nBins, xsec_ttZ);
    HistogramFile ttZ_2018(path+"ttZ_2018.root", obsname, nVarSgnl, nSelections, nSysSgnl, nBins, xsec_ttZ);
    plot.AddSample(ttZ_2016.values[0], 0, 1, false);
    plot.AddSample(ttZ_2017.values[0], 1, 1, false);
    plot.AddSample(ttZ_2018.values[0], 2, 1, false);
    plot.AddSample(ttZ_2016.values[1], 0, 6, false, -1.0);
    plot.AddSample(ttZ_2017.values[1], 1, 6, false, -1.0);
    plot.AddSample(ttZ_2018.values[1], 2, 6, false, -1.0);

    HistogramFile ttH_2016(path+"ttH_2016.root", obsname, nVarBkgr, nSelections, nSysBkgr, nBins, xsec_ttH);
    // HistogramFile ttH_2017(path+"ttH_2017.root", obsname, nVarBkgr, nSelections, nSysBkgr, nBins, xsec_ttH);
    HistogramFile ttH_2018(path+"ttH_2018.root", obsname, nVarBkgr, nSelections, nSysBkgr, nBins, xsec_ttH);
    plot.AddSample(ttH_2016.values[0], 0, 4, false);
    // plot.AddSample(ttH_2017.values[0], 1, 4, false);
    plot.AddSample(ttH_2018.values[0], 2, 4, false);
    plot.AddSample(ttH_2016.values[2], 0, 6, false, -1.0);
    // plot.AddSample(ttH_2017.values[2], 1, 6, false, -1.0);
    plot.AddSample(ttH_2018.values[2], 2, 6, false, -1.0);

    HistogramFile tZq_2016(path+"tZq_2016.root", obsname, nVarBkgr, nSelections, nSysBkgr, nBins, xsec_tZq);
    HistogramFile tZq_2017(path+"tZq_2017.root", obsname, nVarBkgr, nSelections, nSysBkgr, nBins, xsec_tZq);
    HistogramFile tZq_2018(path+"tZq_2018.root", obsname, nVarBkgr, nSelections, nSysBkgr, nBins, xsec_tZq);
    plot.AddSample(tZq_2016.values[0], 0, 2, false);
    plot.AddSample(tZq_2017.values[0], 1, 2, false);
    plot.AddSample(tZq_2018.values[0], 2, 2, false);
    plot.AddSample(tZq_2016.values[2], 0, 6, false, -1.0);
    plot.AddSample(tZq_2017.values[2], 1, 6, false, -1.0);
    plot.AddSample(tZq_2018.values[2], 2, 6, false, -1.0);

    HistogramFile tWZ_2016(path+"tWZ_2016.root", obsname, nVarBkgr, nSelections, nSysBkgr, nBins, xsec_tWZ);
    HistogramFile tWZ_2017(path+"tWZ_2017.root", obsname, nVarBkgr, nSelections, nSysBkgr, nBins, xsec_tWZ);
    HistogramFile tWZ_2018(path+"tWZ_2018.root", obsname, nVarBkgr, nSelections, nSysBkgr, nBins, xsec_tWZ);
    plot.AddSample(tWZ_2016.values[0], 0, 3, false);
    plot.AddSample(tWZ_2017.values[0], 1, 3, false);
    plot.AddSample(tWZ_2018.values[0], 2, 3, false);
    plot.AddSample(tWZ_2016.values[2], 0, 6, false, -1.0);
    plot.AddSample(tWZ_2017.values[2], 1, 6, false, -1.0);
    plot.AddSample(tWZ_2018.values[2], 2, 6, false, -1.0);

    HistogramFile ttW_2016(path+"ttW_2016.root", obsname, nVarBkgr, nSelections, nSysBkgr, nBins, xsec_ttW);
    HistogramFile ttW_2017(path+"ttW_2017.root", obsname, nVarBkgr, nSelections, nSysBkgr, nBins, xsec_ttW);
    HistogramFile ttW_2018(path+"ttW_2018.root", obsname, nVarBkgr, nSelections, nSysBkgr, nBins, xsec_ttW);
    plot.AddSample(ttW_2016.values[0], 0, 4, false);
    plot.AddSample(ttW_2017.values[0], 1, 4, false);
    plot.AddSample(ttW_2018.values[0], 2, 4, false);
    plot.AddSample(ttW_2016.values[2], 0, 6, false, -1.0);
    plot.AddSample(ttW_2017.values[2], 1, 6, false, -1.0);
    plot.AddSample(ttW_2018.values[2], 2, 6, false, -1.0);

    // HistogramFile tHQ_2016(path+"tHQ_2016.root", obsname, nVarBkgr, nSelections, nSysBkgr, nBins, xsec_tHQ);
    // HistogramFile tHQ_2017(path+"tHQ_2017.root", obsname, nVarBkgr, nSelections, nSysBkgr, nBins, xsec_tHQ);
    HistogramFile tHQ_2018(path+"tHQ_2018.root", obsname, nVarBkgr, nSelections, nSysBkgr, nBins, xsec_tHQ);
    // plot.AddSample(tHQ_2016.values[0], 0, 4, false);
    // plot.AddSample(tHQ_2017.values[0], 1, 4, false);
    plot.AddSample(tHQ_2018.values[0], 2, 4, false);
    // plot.AddSample(tHQ_2016.values[2], 0, 6, false, -1.0);
    // plot.AddSample(tHQ_2017.values[2], 1, 6, false, -1.0);
    plot.AddSample(tHQ_2018.values[2], 2, 6, false, -1.0);

    HistogramFile tHW_2016(path+"tHW_2016.root", obsname, nVarBkgr, nSelections, nSysBkgr, nBins, xsec_tHW);
    HistogramFile tHW_2017(path+"tHW_2017.root", obsname, nVarBkgr, nSelections, nSysBkgr, nBins, xsec_tHW);
    HistogramFile tHW_2018(path+"tHW_2018.root", obsname, nVarBkgr, nSelections, nSysBkgr, nBins, xsec_tHW);
    plot.AddSample(tHW_2016.values[0], 0, 4, false);
    plot.AddSample(tHW_2017.values[0], 1, 4, false);
    plot.AddSample(tHW_2018.values[0], 2, 4, false);
    plot.AddSample(tHW_2016.values[2], 0, 6, false, -1.0);
    plot.AddSample(tHW_2017.values[2], 1, 6, false, -1.0);
    plot.AddSample(tHW_2018.values[2], 2, 6, false, -1.0);

    HistogramFile ttZlight_2016(path+"ttZlight_2016.root", obsname, nVarBkgr, nSelections, nSysBkgr, nBins, xsec_ttZlight);
    HistogramFile ttZlight_2017(path+"ttZlight_2017.root", obsname, nVarBkgr, nSelections, nSysBkgr, nBins, xsec_ttZlight);
    HistogramFile ttZlight_2018(path+"ttZlight_2018.root", obsname, nVarBkgr, nSelections, nSysBkgr, nBins, xsec_ttZlight);
    plot.AddSample(ttZlight_2016.values[0], 0, 4, false);
    plot.AddSample(ttZlight_2017.values[0], 1, 4, false);
    plot.AddSample(ttZlight_2018.values[0], 2, 4, false);
    plot.AddSample(ttZlight_2016.values[2], 0, 6, false, -1.0);
    plot.AddSample(ttZlight_2017.values[2], 1, 6, false, -1.0);
    plot.AddSample(ttZlight_2018.values[2], 2, 6, false, -1.0);

    HistogramFile tttt_2016(path+"tttt_2016.root", obsname, nVarBkgr, nSelections, nSysBkgr, nBins, xsec_tttt);
    HistogramFile tttt_2017(path+"tttt_2017.root", obsname, nVarBkgr, nSelections, nSysBkgr, nBins, xsec_tttt);
    HistogramFile tttt_2018(path+"tttt_2018.root", obsname, nVarBkgr, nSelections, nSysBkgr, nBins, xsec_tttt);
    plot.AddSample(tttt_2016.values[0], 0, 4, false);
    plot.AddSample(tttt_2017.values[0], 1, 4, false);
    plot.AddSample(tttt_2018.values[0], 2, 4, false);
    plot.AddSample(tttt_2016.values[2], 0, 6, false, -1.0);
    plot.AddSample(tttt_2017.values[2], 1, 6, false, -1.0);
    plot.AddSample(tttt_2018.values[2], 2, 6, false, -1.0);

    HistogramFile ttWW_2016(path+"ttWW_2016.root", obsname, nVarBkgr, nSelections, nSysBkgr, nBins, xsec_ttWW);
    HistogramFile ttWW_2017(path+"ttWW_2017.root", obsname, nVarBkgr, nSelections, nSysBkgr, nBins, xsec_ttWW);
    HistogramFile ttWW_2018(path+"ttWW_2018.root", obsname, nVarBkgr, nSelections, nSysBkgr, nBins, xsec_ttWW);
    plot.AddSample(ttWW_2016.values[0], 0, 4, false);
    plot.AddSample(ttWW_2017.values[0], 1, 4, false);
    plot.AddSample(ttWW_2018.values[0], 2, 4, false);
    plot.AddSample(ttWW_2016.values[2], 0, 6, false, -1.0);
    plot.AddSample(ttWW_2017.values[2], 1, 6, false, -1.0);
    plot.AddSample(ttWW_2018.values[2], 2, 6, false, -1.0);

    HistogramFile ttWZ_2016(path+"ttWZ_2016.root", obsname, nVarBkgr, nSelections, nSysBkgr, nBins, xsec_ttWZ);
    HistogramFile ttWZ_2017(path+"ttWZ_2017.root", obsname, nVarBkgr, nSelections, nSysBkgr, nBins, xsec_ttWZ);
    HistogramFile ttWZ_2018(path+"ttWZ_2018.root", obsname, nVarBkgr, nSelections, nSysBkgr, nBins, xsec_ttWZ);
    plot.AddSample(ttWZ_2016.values[0], 0, 4, false);
    plot.AddSample(ttWZ_2017.values[0], 1, 4, false);
    plot.AddSample(ttWZ_2018.values[0], 2, 4, false);
    plot.AddSample(ttWZ_2016.values[2], 0, 6, false, -1.0);
    plot.AddSample(ttWZ_2017.values[2], 1, 6, false, -1.0);
    plot.AddSample(ttWZ_2018.values[2], 2, 6, false, -1.0);

    HistogramFile ttZZ_2016(path+"ttZZ_2016.root", obsname, nVarBkgr, nSelections, nSysBkgr, nBins, xsec_ttZZ);
    HistogramFile ttZZ_2017(path+"ttZZ_2017.root", obsname, nVarBkgr, nSelections, nSysBkgr, nBins, xsec_ttZZ);
    HistogramFile ttZZ_2018(path+"ttZZ_2018.root", obsname, nVarBkgr, nSelections, nSysBkgr, nBins, xsec_ttZZ);
    plot.AddSample(ttZZ_2016.values[0], 0, 4, false);
    plot.AddSample(ttZZ_2017.values[0], 1, 4, false);
    plot.AddSample(ttZZ_2018.values[0], 2, 4, false);
    plot.AddSample(ttZZ_2016.values[2], 0, 6, false, -1.0);
    plot.AddSample(ttZZ_2017.values[2], 1, 6, false, -1.0);
    plot.AddSample(ttZZ_2018.values[2], 2, 6, false, -1.0);

    HistogramFile WZ3L_2016(path+"WZ3L_2016.root", obsname, nVarBkgr, nSelections, nSysBkgr, nBins, xsec_WZ3L);
    HistogramFile WZ3L_2017(path+"WZ3L_2017.root", obsname, nVarBkgr, nSelections, nSysBkgr, nBins, xsec_WZ3L);
    HistogramFile WZ3L_2018(path+"WZ3L_2018.root", obsname, nVarBkgr, nSelections, nSysBkgr, nBins, xsec_WZ3L);
    plot.AddSample(WZ3L_2016.values[0], 0, 5, false);
    plot.AddSample(WZ3L_2017.values[0], 1, 5, false);
    plot.AddSample(WZ3L_2018.values[0], 2, 5, false);
    plot.AddSample(WZ3L_2016.values[2], 0, 6, false, -1.0);
    plot.AddSample(WZ3L_2017.values[2], 1, 6, false, -1.0);
    plot.AddSample(WZ3L_2018.values[2], 2, 6, false, -1.0);

    // HistogramFile WZ2L_2016(path+"WZ2L_2016.root", obsname, nVarBkgr, nSelections, nSysBkgr, nBins, xsec_WZ2L);
    // HistogramFile WZ2L_2017(path+"WZ2L_2017.root", obsname, nVarBkgr, nSelections, nSysBkgr, nBins, xsec_WZ2L);
    // HistogramFile WZ2L_2018(path+"WZ2L_2018.root", obsname, nVarBkgr, nSelections, nSysBkgr, nBins, xsec_WZ2L);
    // plot.AddSample(WZ2L_2016.values[0], 0, 5, false);
    // plot.AddSample(WZ2L_2017.values[0], 1, 5, false);
    // plot.AddSample(WZ2L_2018.values[0], 2, 5, false);
    // plot.AddSample(WZ2L_2016.values[2], 0, 6, false, -1.0);
    // plot.AddSample(WZ2L_2017.values[2], 1, 6, false, -1.0);
    // plot.AddSample(WZ2L_2018.values[2], 2, 6, false, -1.0);

    HistogramFile DY_2016(path+"DY_2016.root", obsname, nVarBkgr, nSelections, nSysBkgr, nBins, xsec_DY);
    HistogramFile DY_2017(path+"DY_2017.root", obsname, nVarBkgr, nSelections, nSysBkgr, nBins, xsec_DY);
    HistogramFile DY_2018(path+"DY_2018.root", obsname, nVarBkgr, nSelections, nSysBkgr, nBins, xsec_DY);
    plot.AddSample(DY_2016.values[0], 0, 5, false);
    plot.AddSample(DY_2017.values[0], 1, 5, false);
    plot.AddSample(DY_2018.values[0], 2, 5, false);
    plot.AddSample(DY_2016.values[2], 0, 6, false, -1.0);
    plot.AddSample(DY_2017.values[2], 1, 6, false, -1.0);
    plot.AddSample(DY_2018.values[2], 2, 6, false, -1.0);

    HistogramFile tG_2016(path+"tG_2016.root", obsname, nVarBkgr, nSelections, nSysBkgr, nBins, xsec_tG);
    HistogramFile tG_2017(path+"tG_2017.root", obsname, nVarBkgr, nSelections, nSysBkgr, nBins, xsec_tG);
    HistogramFile tG_2018(path+"tG_2018.root", obsname, nVarBkgr, nSelections, nSysBkgr, nBins, xsec_tG);
    plot.AddSample(tG_2016.values[0], 0, 4, false);
    plot.AddSample(tG_2017.values[0], 1, 4, false);
    plot.AddSample(tG_2018.values[0], 2, 4, false);
    plot.AddSample(tG_2016.values[2], 0, 6, false, -1.0);
    plot.AddSample(tG_2017.values[2], 1, 6, false, -1.0);
    plot.AddSample(tG_2018.values[2], 2, 6, false, -1.0);

    HistogramFile ttG_2016(path+"ttG_2016.root", obsname, nVarBkgr, nSelections, nSysBkgr, nBins, xsec_ttG);
    HistogramFile ttG_2017(path+"ttG_2017.root", obsname, nVarBkgr, nSelections, nSysBkgr, nBins, xsec_ttG);
    HistogramFile ttG_2018(path+"ttG_2018.root", obsname, nVarBkgr, nSelections, nSysBkgr, nBins, xsec_ttG);
    plot.AddSample(ttG_2016.values[0], 0, 4, false);
    plot.AddSample(ttG_2017.values[0], 1, 4, false);
    plot.AddSample(ttG_2018.values[0], 2, 4, false);
    plot.AddSample(ttG_2016.values[2], 0, 6, false, -1.0);
    plot.AddSample(ttG_2017.values[2], 1, 6, false, -1.0);
    plot.AddSample(ttG_2018.values[2], 2, 6, false, -1.0);

    HistogramFile WG_2016(path+"WG_2016.root", obsname, nVarBkgr, nSelections, nSysBkgr, nBins, xsec_WG);
    HistogramFile WG_2017(path+"WG_2017.root", obsname, nVarBkgr, nSelections, nSysBkgr, nBins, xsec_WG);
    HistogramFile WG_2018(path+"WG_2018.root", obsname, nVarBkgr, nSelections, nSysBkgr, nBins, xsec_WG);
    plot.AddSample(WG_2016.values[0], 0, 5, false);
    plot.AddSample(WG_2017.values[0], 1, 5, false);
    plot.AddSample(WG_2018.values[0], 2, 5, false);
    plot.AddSample(WG_2016.values[2], 0, 6, false, -1.0);
    plot.AddSample(WG_2017.values[2], 1, 6, false, -1.0);
    plot.AddSample(WG_2018.values[2], 2, 6, false, -1.0);

    // HistogramFile tt_2016(path+"tt_2016.root", obsname, nVarBkgr, nSelections, nSysBkgr, nBins, xsec_tt);
    // HistogramFile tt_2017(path+"tt_2017.root", obsname, nVarBkgr, nSelections, nSysBkgr, nBins, xsec_tt);
    // HistogramFile tt_2018(path+"tt_2018.root", obsname, nVarBkgr, nSelections, nSysBkgr, nBins, xsec_tt);
    // plot.AddSample(tt_2016.values[0], 0, 4, false);
    // plot.AddSample(tt_2017.values[0], 1, 4, false);
    // plot.AddSample(tt_2018.values[0], 2, 4, false);
    // plot.AddSample(tt_2016.values[2], 0, 6, false, -1.0);
    // plot.AddSample(tt_2017.values[2], 1, 6, false, -1.0);
    // plot.AddSample(tt_2018.values[2], 2, 6, false, -1.0);

    HistogramFile ZZ4L_2016(path+"ZZ4L_2016.root", obsname, nVarBkgr, nSelections, nSysBkgr, nBins, xsec_ZZ4L);
    HistogramFile ZZ4L_2017(path+"ZZ4L_2017.root", obsname, nVarBkgr, nSelections, nSysBkgr, nBins, xsec_ZZ4L);
    HistogramFile ZZ4L_2018(path+"ZZ4L_2018.root", obsname, nVarBkgr, nSelections, nSysBkgr, nBins, xsec_ZZ4L);
    plot.AddSample(ZZ4L_2016.values[0], 0, 5, false);
    plot.AddSample(ZZ4L_2017.values[0], 1, 5, false);
    plot.AddSample(ZZ4L_2018.values[0], 2, 5, false);
    plot.AddSample(ZZ4L_2016.values[2], 0, 6, false, -1.0);
    plot.AddSample(ZZ4L_2017.values[2], 1, 6, false, -1.0);
    plot.AddSample(ZZ4L_2018.values[2], 2, 6, false, -1.0);

    HistogramFile ZZ2E2M_2016(path+"ZZ2E2M_2016.root", obsname, nVarBkgr, nSelections, nSysBkgr, nBins, xsec_ZZ2E2M);
    HistogramFile ZZ2E2M_2017(path+"ZZ2E2M_2017.root", obsname, nVarBkgr, nSelections, nSysBkgr, nBins, xsec_ZZ2E2M);
    HistogramFile ZZ2E2M_2018(path+"ZZ2E2M_2018.root", obsname, nVarBkgr, nSelections, nSysBkgr, nBins, xsec_ZZ2E2M);
    plot.AddSample(ZZ2E2M_2016.values[0], 0, 5, false);
    plot.AddSample(ZZ2E2M_2017.values[0], 1, 5, false);
    plot.AddSample(ZZ2E2M_2018.values[0], 2, 5, false);
    plot.AddSample(ZZ2E2M_2016.values[2], 0, 6, false, -1.0);
    plot.AddSample(ZZ2E2M_2017.values[2], 1, 6, false, -1.0);
    plot.AddSample(ZZ2E2M_2018.values[2], 2, 6, false, -1.0);

    HistogramFile ZZ2E2T_2016(path+"ZZ2E2T_2016.root", obsname, nVarBkgr, nSelections, nSysBkgr, nBins, xsec_ZZ2E2T);
    HistogramFile ZZ2E2T_2017(path+"ZZ2E2T_2017.root", obsname, nVarBkgr, nSelections, nSysBkgr, nBins, xsec_ZZ2E2T);
    HistogramFile ZZ2E2T_2018(path+"ZZ2E2T_2018.root", obsname, nVarBkgr, nSelections, nSysBkgr, nBins, xsec_ZZ2E2T);
    plot.AddSample(ZZ2E2T_2016.values[0], 0, 5, false);
    plot.AddSample(ZZ2E2T_2017.values[0], 1, 5, false);
    plot.AddSample(ZZ2E2T_2018.values[0], 2, 5, false);
    plot.AddSample(ZZ2E2T_2016.values[2], 0, 6, false, -1.0);
    plot.AddSample(ZZ2E2T_2017.values[2], 1, 6, false, -1.0);
    plot.AddSample(ZZ2E2T_2018.values[2], 2, 6, false, -1.0);

    HistogramFile ZZ2M2T_2016(path+"ZZ2M2T_2016.root", obsname, nVarBkgr, nSelections, nSysBkgr, nBins, xsec_ZZ2M2T);
    HistogramFile ZZ2M2T_2017(path+"ZZ2M2T_2017.root", obsname, nVarBkgr, nSelections, nSysBkgr, nBins, xsec_ZZ2M2T);
    HistogramFile ZZ2M2T_2018(path+"ZZ2M2T_2018.root", obsname, nVarBkgr, nSelections, nSysBkgr, nBins, xsec_ZZ2M2T);
    plot.AddSample(ZZ2M2T_2016.values[0], 0, 5, false);
    plot.AddSample(ZZ2M2T_2017.values[0], 1, 5, false);
    plot.AddSample(ZZ2M2T_2018.values[0], 2, 5, false);
    plot.AddSample(ZZ2M2T_2016.values[2], 0, 6, false, -1.0);
    plot.AddSample(ZZ2M2T_2017.values[2], 1, 6, false, -1.0);
    plot.AddSample(ZZ2M2T_2018.values[2], 2, 6, false, -1.0);

    HistogramFile ZZ4E_2016(path+"ZZ4E_2016.root", obsname, nVarBkgr, nSelections, nSysBkgr, nBins, xsec_ZZ4E);
    HistogramFile ZZ4E_2017(path+"ZZ4E_2017.root", obsname, nVarBkgr, nSelections, nSysBkgr, nBins, xsec_ZZ4E);
    HistogramFile ZZ4E_2018(path+"ZZ4E_2018.root", obsname, nVarBkgr, nSelections, nSysBkgr, nBins, xsec_ZZ4E);
    plot.AddSample(ZZ4E_2016.values[0], 0, 5, false);
    plot.AddSample(ZZ4E_2017.values[0], 1, 5, false);
    plot.AddSample(ZZ4E_2018.values[0], 2, 5, false);
    plot.AddSample(ZZ4E_2016.values[2], 0, 6, false, -1.0);
    plot.AddSample(ZZ4E_2017.values[2], 1, 6, false, -1.0);
    plot.AddSample(ZZ4E_2018.values[2], 2, 6, false, -1.0);

    HistogramFile ZZ4M_2016(path+"ZZ4M_2016.root", obsname, nVarBkgr, nSelections, nSysBkgr, nBins, xsec_ZZ4M);
    HistogramFile ZZ4M_2017(path+"ZZ4M_2017.root", obsname, nVarBkgr, nSelections, nSysBkgr, nBins, xsec_ZZ4M);
    HistogramFile ZZ4M_2018(path+"ZZ4M_2018.root", obsname, nVarBkgr, nSelections, nSysBkgr, nBins, xsec_ZZ4M);
    plot.AddSample(ZZ4M_2016.values[0], 0, 5, false);
    plot.AddSample(ZZ4M_2017.values[0], 1, 5, false);
    plot.AddSample(ZZ4M_2018.values[0], 2, 5, false);
    plot.AddSample(ZZ4M_2016.values[2], 0, 6, false, -1.0);
    plot.AddSample(ZZ4M_2017.values[2], 1, 6, false, -1.0);
    plot.AddSample(ZZ4M_2018.values[2], 2, 6, false, -1.0);

    HistogramFile ggHZZ_2016(path+"ggHZZ_2016.root", obsname, nVarBkgr, nSelections, nSysBkgr, nBins, xsec_ggHZZ);
    HistogramFile ggHZZ_2017(path+"ggHZZ_2017.root", obsname, nVarBkgr, nSelections, nSysBkgr, nBins, xsec_ggHZZ);
    HistogramFile ggHZZ_2018(path+"ggHZZ_2018.root", obsname, nVarBkgr, nSelections, nSysBkgr, nBins, xsec_ggHZZ);
    plot.AddSample(ggHZZ_2016.values[0], 0, 5, false);
    plot.AddSample(ggHZZ_2017.values[0], 1, 5, false);
    plot.AddSample(ggHZZ_2018.values[0], 2, 5, false);
    plot.AddSample(ggHZZ_2016.values[2], 0, 6, false, -1.0);
    plot.AddSample(ggHZZ_2017.values[2], 1, 6, false, -1.0);
    plot.AddSample(ggHZZ_2018.values[2], 2, 6, false, -1.0);

    HistogramFile VBFHZZ_2016(path+"VBFHZZ_2016.root", obsname, nVarBkgr, nSelections, nSysBkgr, nBins, xsec_VBFHZZ);
    HistogramFile VBFHZZ_2017(path+"VBFHZZ_2017.root", obsname, nVarBkgr, nSelections, nSysBkgr, nBins, xsec_VBFHZZ);
    HistogramFile VBFHZZ_2018(path+"VBFHZZ_2018.root", obsname, nVarBkgr, nSelections, nSysBkgr, nBins, xsec_VBFHZZ);
    plot.AddSample(VBFHZZ_2016.values[0], 0, 5, false);
    plot.AddSample(VBFHZZ_2017.values[0], 1, 5, false);
    plot.AddSample(VBFHZZ_2018.values[0], 2, 5, false);
    plot.AddSample(VBFHZZ_2016.values[2], 0, 6, false, -1.0);
    plot.AddSample(VBFHZZ_2017.values[2], 1, 6, false, -1.0);
    plot.AddSample(VBFHZZ_2018.values[2], 2, 6, false, -1.0);

    // HistogramFile WpHZZ_2016(path+"WpHZZ_2016.root", obsname, nVarBkgr, nSelections, nSysBkgr, nBins, xsec_WpHZZ);
    // HistogramFile WpHZZ_2017(path+"WpHZZ_2017.root", obsname, nVarBkgr, nSelections, nSysBkgr, nBins, xsec_WpHZZ);
    // HistogramFile WpHZZ_2018(path+"WpHZZ_2018.root", obsname, nVarBkgr, nSelections, nSysBkgr, nBins, xsec_WpHZZ);
    // plot.AddSample(WpHZZ_2016.values[0], 0, 5, false);
    // plot.AddSample(WpHZZ_2017.values[0], 1, 5, false);
    // plot.AddSample(WpHZZ_2018.values[0], 2, 5, false);
    // plot.AddSample(WpHZZ_2016.values[2], 0, 6, false, -1.0);
    // plot.AddSample(WpHZZ_2017.values[2], 1, 6, false, -1.0);
    // plot.AddSample(WpHZZ_2018.values[2], 2, 6, false, -1.0);

    HistogramFile WmHZZ_2016(path+"WmHZZ_2016.root", obsname, nVarBkgr, nSelections, nSysBkgr, nBins, xsec_WmHZZ);
    HistogramFile WmHZZ_2017(path+"WmHZZ_2017.root", obsname, nVarBkgr, nSelections, nSysBkgr, nBins, xsec_WmHZZ);
    HistogramFile WmHZZ_2018(path+"WmHZZ_2018.root", obsname, nVarBkgr, nSelections, nSysBkgr, nBins, xsec_WmHZZ);
    plot.AddSample(WmHZZ_2016.values[0], 0, 5, false);
    plot.AddSample(WmHZZ_2017.values[0], 1, 5, false);
    plot.AddSample(WmHZZ_2018.values[0], 2, 5, false);
    plot.AddSample(WmHZZ_2016.values[2], 0, 6, false, -1.0);
    plot.AddSample(WmHZZ_2017.values[2], 1, 6, false, -1.0);
    plot.AddSample(WmHZZ_2018.values[2], 2, 6, false, -1.0);

    HistogramFile ZHZZ_2016(path+"ZHZZ_2016.root", obsname, nVarBkgr, nSelections, nSysBkgr, nBins, xsec_ZHZZ);
    HistogramFile ZHZZ_2017(path+"ZHZZ_2017.root", obsname, nVarBkgr, nSelections, nSysBkgr, nBins, xsec_ZHZZ);
    HistogramFile ZHZZ_2018(path+"ZHZZ_2018.root", obsname, nVarBkgr, nSelections, nSysBkgr, nBins, xsec_ZHZZ);
    plot.AddSample(ZHZZ_2016.values[0], 0, 5, false);
    plot.AddSample(ZHZZ_2017.values[0], 1, 5, false);
    plot.AddSample(ZHZZ_2018.values[0], 2, 5, false);
    plot.AddSample(ZHZZ_2016.values[2], 0, 6, false, -1.0);
    plot.AddSample(ZHZZ_2017.values[2], 1, 6, false, -1.0);
    plot.AddSample(ZHZZ_2018.values[2], 2, 6, false, -1.0);

    HistogramFile WZG_2016(path+"WZG_2016.root", obsname, nVarBkgr, nSelections, nSysBkgr, nBins, xsec_WZG);
    HistogramFile WZG_2017(path+"WZG_2017.root", obsname, nVarBkgr, nSelections, nSysBkgr, nBins, xsec_WZG);
    HistogramFile WZG_2018(path+"WZG_2018.root", obsname, nVarBkgr, nSelections, nSysBkgr, nBins, xsec_WZG);
    plot.AddSample(WZG_2016.values[0], 0, 5, false);
    plot.AddSample(WZG_2017.values[0], 1, 5, false);
    plot.AddSample(WZG_2018.values[0], 2, 5, false);
    plot.AddSample(WZG_2016.values[2], 0, 6, false, -1.0);
    plot.AddSample(WZG_2017.values[2], 1, 6, false, -1.0);
    plot.AddSample(WZG_2018.values[2], 2, 6, false, -1.0);

    HistogramFile ZZZ_2016(path+"ZZZ_2016.root", obsname, nVarBkgr, nSelections, nSysBkgr, nBins, xsec_ZZZ);
    HistogramFile ZZZ_2017(path+"ZZZ_2017.root", obsname, nVarBkgr, nSelections, nSysBkgr, nBins, xsec_ZZZ);
    HistogramFile ZZZ_2018(path+"ZZZ_2018.root", obsname, nVarBkgr, nSelections, nSysBkgr, nBins, xsec_ZZZ);
    plot.AddSample(ZZZ_2016.values[0], 0, 5, false);
    plot.AddSample(ZZZ_2017.values[0], 1, 5, false);
    plot.AddSample(ZZZ_2018.values[0], 2, 5, false);
    plot.AddSample(ZZZ_2016.values[2], 0, 6, false, -1.0);
    plot.AddSample(ZZZ_2017.values[2], 1, 6, false, -1.0);
    plot.AddSample(ZZZ_2018.values[2], 2, 6, false, -1.0);

    HistogramFile WZZ_2016(path+"WZZ_2016.root", obsname, nVarBkgr, nSelections, nSysBkgr, nBins, xsec_WZZ);
    HistogramFile WZZ_2017(path+"WZZ_2017.root", obsname, nVarBkgr, nSelections, nSysBkgr, nBins, xsec_WZZ);
    HistogramFile WZZ_2018(path+"WZZ_2018.root", obsname, nVarBkgr, nSelections, nSysBkgr, nBins, xsec_WZZ);
    plot.AddSample(WZZ_2016.values[0], 0, 5, false);
    plot.AddSample(WZZ_2017.values[0], 1, 5, false);
    plot.AddSample(WZZ_2018.values[0], 2, 5, false);
    plot.AddSample(WZZ_2016.values[2], 0, 6, false, -1.0);
    plot.AddSample(WZZ_2017.values[2], 1, 6, false, -1.0);
    plot.AddSample(WZZ_2018.values[2], 2, 6, false, -1.0);

    HistogramFile WWZ_2016(path+"WWZ_2016.root", obsname, nVarBkgr, nSelections, nSysBkgr, nBins, xsec_WWZ);
    HistogramFile WWZ_2017(path+"WWZ_2017.root", obsname, nVarBkgr, nSelections, nSysBkgr, nBins, xsec_WWZ);
    HistogramFile WWZ_2018(path+"WWZ_2018.root", obsname, nVarBkgr, nSelections, nSysBkgr, nBins, xsec_WWZ);
    plot.AddSample(WWZ_2016.values[0], 0, 5, false);
    plot.AddSample(WWZ_2017.values[0], 1, 5, false);
    plot.AddSample(WWZ_2018.values[0], 2, 5, false);
    plot.AddSample(WWZ_2016.values[2], 0, 6, false, -1.0);
    plot.AddSample(WWZ_2017.values[2], 1, 6, false, -1.0);
    plot.AddSample(WWZ_2018.values[2], 2, 6, false, -1.0);

    HistogramFile WWW_2016(path+"WWW_2016.root", obsname, nVarBkgr, nSelections, nSysBkgr, nBins, xsec_WWW);
    HistogramFile WWW_2017(path+"WWW_2017.root", obsname, nVarBkgr, nSelections, nSysBkgr, nBins, xsec_WWW);
    HistogramFile WWW_2018(path+"WWW_2018.root", obsname, nVarBkgr, nSelections, nSysBkgr, nBins, xsec_WWW);
    plot.AddSample(WWW_2016.values[0], 0, 5, false);
    plot.AddSample(WWW_2017.values[0], 1, 5, false);
    plot.AddSample(WWW_2018.values[0], 2, 5, false);
    plot.AddSample(WWW_2016.values[2], 0, 6, false, -1.0);
    plot.AddSample(WWW_2017.values[2], 1, 6, false, -1.0);
    plot.AddSample(WWW_2018.values[2], 2, 6, false, -1.0);

    HistogramFile WWDS_2016(path+"WWDS_2016.root", obsname, nVarBkgr, nSelections, nSysBkgr, nBins, xsec_WWDS);
    HistogramFile WWDS_2017(path+"WWDS_2017.root", obsname, nVarBkgr, nSelections, nSysBkgr, nBins, xsec_WWDS);
    HistogramFile WWDS_2018(path+"WWDS_2018.root", obsname, nVarBkgr, nSelections, nSysBkgr, nBins, xsec_WWDS);
    plot.AddSample(WWDS_2016.values[0], 0, 5, false);
    plot.AddSample(WWDS_2017.values[0], 1, 5, false);
    plot.AddSample(WWDS_2018.values[0], 2, 5, false);
    plot.AddSample(WWDS_2016.values[2], 0, 6, false, -1.0);
    plot.AddSample(WWDS_2017.values[2], 1, 6, false, -1.0);
    plot.AddSample(WWDS_2018.values[2], 2, 6, false, -1.0);

    HistogramFile WW_2016(path+"WW_2016.root", obsname, nVarBkgr, nSelections, nSysBkgr, nBins, xsec_WW);
    // HistogramFile WW_2017(path+"WW_2017.root", obsname, nVarBkgr, nSelections, nSysBkgr, nBins, xsec_WW);
    // HistogramFile WW_2018(path+"WW_2018.root", obsname, nVarBkgr, nSelections, nSysBkgr, nBins, xsec_WW);
    plot.AddSample(WW_2016.values[0], 0, 5, false);
    // plot.AddSample(WW_2017.values[0], 1, 5, false);
    // plot.AddSample(WW_2018.values[0], 2, 5, false);
    plot.AddSample(WW_2016.values[2], 0, 6, false, -1.0);
    // plot.AddSample(WW_2017.values[2], 1, 6, false, -1.0);
    // plot.AddSample(WW_2018.values[2], 2, 6, false, -1.0);

    plot.Evaluate();

    for(int iYear=0; iYear<=nYears; iYear++) {
        const std::string sYear = (iYear==0) ? "2016" : (iYear==1) ? "2017" : (iYear==2) ? "2018" : "run2";
        const std::string sLumi = (iYear==0) ? "35.9" : (iYear==1) ? "41.5" : (iYear==2) ? "59.7" : "137";
        for(int iSel=0; iSel<nSelections; iSel++) {
            const std::string sTitle = (iSel==0) ? "≥3 jets" : (iSel==1) ? "≥4 jets" : "=3 jets";
            const std::string sSelection = "sel"+std::to_string(iSel);
            const std::string data = plot.Print(iYear, iSel);
            const std::string filename = "plot_"+sYear+"_"+sSelection+"_"+obsname;
            const std::string fullpath = outputdir+"/"+filename+".txt";
            std::cout << "Writing " << fullpath << std::endl;
            std::ofstream out(fullpath);
            out << data << std::endl
                << "filename = " << '"' << filename << '"' << std::endl
                << "my_xlabel = " << '"' << obsname << '"' << std::endl
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
