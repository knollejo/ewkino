//
//
//include c++ library classes

//include ROOT classes
#include "TMath.h"
#include "TTree.h"

//include general parts of framework
#include "../Tools/interface/analysisTools.h"
#include "../Tools/interface/systemTools.h"
#include "../Tools/interface/stringTools.h"
#include "../weights/interface/ConcreteReweighterFactory.h"

//include ttZ specific code
#include "interface/ttZHistograms.h"
#include "interface/ttZSelection.h"
#include "interface/ttZTruth.h"
#include "interface/ttZVariables.h"

#include <functional>


void analyze( const std::string& year, const std::string& sampleDirectoryPath , const std::string& procName , const std::string& subfolder ){

    //check setup
    analysisTools::checkYearString( year );
    const bool is_ttZ = (procName=="ttZ1" || procName=="ttZ2" || procName=="ttZ3" || procName=="ttZ4" || procName=="ttZ5");
    const bool is_data = (procName == "data");
    const bool is_2016 = (year=="2016");
    const bool is_2017 = (year=="2017");

    //build TreeReader and loop over samples
    std::cout << "building treeReader" << std::endl;
    TreeReader treeReader( "sampleLists/samples_ttZ_" + year + ".txt", sampleDirectoryPath );
    std::vector< Sample > sampleVec = treeReader.sampleVector();

    //build ttZ reweighter
    std::cout << "building reweighter" << std::endl;
    std::shared_ptr< ReweighterFactory >reweighterFactory( new ttZReweighterFactory() );
    CombinedReweighter reweighter = reweighterFactory->buildReweighter( "../weights/", year, treeReader.sampleVector() );

    //read FR maps
    std::cout << "building FR maps" << std::endl;
    TFile* frFileMuons = TFile::Open( ( "frMaps/fakeRateMap_data_muon_" + year + "_mT.root" ).c_str() );
    std::shared_ptr< TH2 > frMapMuons = std::shared_ptr< TH2 >( dynamic_cast< TH2* >( frFileMuons->Get( ("fakeRate_muon_" + year ).c_str() ) ) );
    frMapMuons->SetDirectory( gROOT );
    frFileMuons->Close();
    TFile* frFileElectrons = TFile::Open( ( "frMaps/fakeRateMap_data_electron_" + year + "_mT.root" ).c_str() );
    std::shared_ptr< TH2 > frMapElectrons = std::shared_ptr< TH2 >( dynamic_cast< TH2* >( frFileElectrons->Get( ( "fakeRate_electron_" + year ).c_str() ) ) );
    frMapElectrons->SetDirectory( gROOT );
    frFileElectrons->Close();

    //initialize kinematic reconstruction
    KinFitter* fitter = new ttZKinFitter(is_2016, is_2017);

    //initialize histograms
    std::cout << "building histograms" << std::endl;
    const int nVariants = is_ttZ ? 18 : is_data ? 2 : 3;
    const int nSelections = 3;
    const int nPdfVariations = 100;
    const int nSystematics = is_data ? 1 : 1+28+nPdfVariations;
    ttZHistograms histograms(nVariants, nSelections, nSystematics);

    //event loop
    std::cout << "event loop" << std::endl;
    for( unsigned sampleIndex = 0; sampleIndex < treeReader.numberOfSamples(); ++sampleIndex ){
        treeReader.initSample();
        if ( sampleVec[ sampleIndex ].processName() != procName ) continue;

        std::cout << treeReader.currentSample().fileName() << std::endl;

        bool sampleHasPdfAndScale = true;

        for( long unsigned entry = 0; entry < treeReader.numberOfEntries(); ++entry ){
            Event event = treeReader.buildEvent( entry );

            //check if sample has pdf and scale information stored, try-catch blocks for every event are too slow
            if( entry == 0 && treeReader.isMC() ){
                try{
                    event.generatorInfo().relativeWeightPdfVar( 99 );
                }catch( std::out_of_range& ){
                    sampleHasPdfAndScale = false;
                }
            }

            //for ttZ events: analyze gen-level information
            auto truthInfo = is_ttZ ? ttZ::evaluateTruthStatus(treeReader) : ttZ::ttzTruth();
            int iVariantShift = 0;
            if(is_ttZ) {
                if(!truthInfo.hasLeptonicZ) iVariantShift = 15; // tt+nunu
                else if(!truthInfo.isOnShell) iVariantShift = 12; // tt+ll off-shell
                else if(truthInfo.hasTaus) iVariantShift = 9; // ttZ with tau in the decay
                else if(truthInfo.nLeptons==4) iVariantShift = 0; // ttZ with 4 leptons
                else if(truthInfo.nLeptons==3) iVariantShift = 3; // ttZ with 3 leptons
                else if(truthInfo.nLeptons==2) iVariantShift = 6; // ttZ with 2 leptons
                else iVariantShift = 15;
            }

            //apply baseline selection
            if( !ttZ::passBaselineSelection( event, true, true ) ) continue;

            //apply lepton pT cuts
            if( !ttZ::passPtCuts( event ) ) continue;

            //require triggers
            if( !ttZ::passTriggerSelection( event ) ) continue;
            if( !( event.passMetFilters() ) ) continue;

            //require the right number of tight and FO(loose) leptons in 3(4) lepton events
            if( !ttZ::passSelectionLNumber( event ) ) continue;

            //require MC events to only contain prompt leptons
            bool is_prompt = true;
            if( event.isMC() && !ttZ::leptonsArePrompt( event ) ){
                is_prompt = false;
            }

            //apply scale-factors and reweighting
            double weight = event.weight();
            if( event.isMC() ){
                weight *= reweighter.totalWeight( event );
            }

            //apply fake-rate weight
            bool is_tight = true;
            if( !ttZ::leptonsAreTight( event ) ){
                is_tight = false;
                weight *= ttZ::fakeRateWeight( event, frMapMuons, frMapElectrons );
            }
            if( !is_data && !is_prompt && !is_tight ) continue;

            //fill nominal histograms
            const int iVariant = is_data ? (is_tight ? 0 : 1)
                               : (is_prompt ? (is_tight ? 0 : 2) : 1)
                               + iVariantShift;
            int passed_selection = ttZ::passSelectionTTZ( event, "nominal" );
            if( passed_selection >= 0 ){

                auto lepVars = ttZ::computeLeptonVariables(event);
                auto jetVars = ttZ::computeJetVariables(event, "nominal");
                auto recVars = ttZ::performKinematicReconstruction(event, "nominal", fitter);
                auto fourVars = ttZ::performFourLeptonsComputation(event);

                histograms.SetValues(lepVars, jetVars, recVars, fourVars);
                histograms.Fill(weight, iVariant, passed_selection, 0);

                //no uncertainties for data
                if( event.isData() ) continue;

                //fill scale variation histograms
                const double weightScaleRenDown = sampleHasPdfAndScale ? event.generatorInfo().relativeWeight_MuR_0p5_MuF_1() : 1.0;
                histograms.Fill(weight*weightScaleRenDown, iVariant, passed_selection, 7);
                const double weightScaleRenUp = sampleHasPdfAndScale ? event.generatorInfo().relativeWeight_MuR_2_MuF_1() : 1.0;
                histograms.Fill(weight*weightScaleRenUp, iVariant, passed_selection, 8);
                const double weightScaleFacDown = sampleHasPdfAndScale ? event.generatorInfo().relativeWeight_MuR_1_MuF_0p5() : 1.0;
                histograms.Fill(weight*weightScaleFacDown, iVariant, passed_selection, 9);
                const double weightScaleFacUp = sampleHasPdfAndScale ? event.generatorInfo().relativeWeight_MuR_1_MuF_2() : 1.0;
                histograms.Fill(weight*weightScaleFacUp, iVariant, passed_selection, 10);
                const double weightScaleBothDown = sampleHasPdfAndScale ? event.generatorInfo().relativeWeight_MuR_0p5_MuF_0p5() : 1.0;
                histograms.Fill(weight*weightScaleBothDown, iVariant, passed_selection, 11);
                const double weightScaleBothUp = sampleHasPdfAndScale ? event.generatorInfo().relativeWeight_MuR_2_MuF_2() : 1.0;
                histograms.Fill(weight*weightScaleBothUp, iVariant, passed_selection, 12);

                //fill pdf histograms
                for( unsigned pdf_i = 0; pdf_i < nPdfVariations; ++pdf_i){
                    const double weightPdf = sampleHasPdfAndScale ? event.generatorInfo().relativeWeightPdfVar(pdf_i) : 1.;
                    histograms.Fill(weight*weightPdf, iVariant, passed_selection, 27+pdf_i);
                }

                //fill pileup histograms
                double weightPileupDown = reweighter[ "pileup" ]->weightDown( event ) / reweighter[ "pileup" ]->weight( event );
                if ( std::isnan(weightPileupDown) ) weightPileupDown = 0.;
                histograms.Fill(weight*weightPileupDown, iVariant, passed_selection, 13);
                double weightPileupUp = reweighter[ "pileup" ]->weightUp( event ) / reweighter[ "pileup" ]->weight( event );
                if ( std::isnan(weightPileupUp) ) weightPileupUp = 0.;
                histograms.Fill(weight*weightPileupUp, iVariant, passed_selection, 14);

                //fill b-tag histograms
                const double weightBTagHeavyDown = reweighter[ "bTag_heavy" ]->weightDown( event ) / reweighter[ "bTag_heavy" ]->weight( event );
                histograms.Fill(weight*weightBTagHeavyDown, iVariant, passed_selection, 15);
                const double weightBTagHeavyUp = reweighter[ "bTag_heavy" ]->weightUp( event ) / reweighter[ "bTag_heavy" ]->weight( event );
                histograms.Fill(weight*weightBTagHeavyUp, iVariant, passed_selection, 16);
                const double weightBTagLightDown = reweighter[ "bTag_light" ]->weightDown( event ) / reweighter[ "bTag_light" ]->weight( event );
                histograms.Fill(weight*weightBTagLightDown, iVariant, passed_selection, 17);
                const double weightBTagLightUp = reweighter[ "bTag_light" ]->weightUp( event ) / reweighter[ "bTag_light" ]->weight( event );
                histograms.Fill(weight*weightBTagLightUp, iVariant, passed_selection, 18);

                //fill prefiring histograms
                const double weightPrefireDown = reweighter[ "prefire" ]->weightDown( event ) / reweighter[ "prefire" ]->weight( event );
                histograms.Fill(weight*weightPrefireDown, iVariant, passed_selection, 19);
                const double weightPrefireUp = reweighter[ "prefire" ]->weightUp( event ) / reweighter[ "prefire" ]->weight( event );
                histograms.Fill(weight*weightPrefireUp, iVariant, passed_selection, 20);

                //fill lepton histograms
                const double leptonIDWeightDownMuon = reweighter[ "muonID" ]->weightDown( event ) / reweighter[ "muonID" ]->weight( event );
                histograms.Fill(weight*leptonIDWeightDownMuon, iVariant, passed_selection, 21);
                const double leptonIDWeightUpMuon = reweighter[ "muonID" ]->weightUp( event ) / reweighter[ "muonID" ]->weight( event );
                histograms.Fill(weight*leptonIDWeightUpMuon, iVariant, passed_selection, 22);
                const double recoWeightDownElec = reweighter[ "electronReco_pTBelow20" ]->weightDown( event ) * reweighter[ "electronReco_pTAbove20" ]->weightDown( event ) / ( reweighter[ "electronReco_pTBelow20" ]->weight( event ) * reweighter[ "electronReco_pTAbove20" ]->weight( event ) );
                histograms.Fill(weight*recoWeightDownElec, iVariant, passed_selection, 23);
                const double recoWeightUpElec = reweighter[ "electronReco_pTBelow20" ]->weightUp( event ) * reweighter[ "electronReco_pTAbove20" ]->weightUp( event ) / ( reweighter[ "electronReco_pTBelow20" ]->weight( event ) * reweighter[ "electronReco_pTAbove20" ]->weight( event ) );
                histograms.Fill(weight*recoWeightUpElec, iVariant, passed_selection, 24);
                const double leptonIDWeightDownElec = reweighter[ "electronID" ]->weightDown( event ) / reweighter[ "electronID" ]->weight( event );
                histograms.Fill(weight*leptonIDWeightDownElec, iVariant, passed_selection, 25);
                const double leptonIDWeightUpElec = reweighter[ "electronID" ]->weightUp( event ) / reweighter[ "electronID" ]->weight( event );
                histograms.Fill(weight*leptonIDWeightUpElec, iVariant, passed_selection, 26);

            }

            //no uncertainties for data
            if( event.isData() ) continue;

            //fill JEC down histograms
            passed_selection = ttZ::passSelectionTTZ( event, "JECDown" );
            if( passed_selection >= 0 ){

                auto lepVars = ttZ::computeLeptonVariables(event);
                auto jetVars = ttZ::computeJetVariables(event, "JECDown");
                auto recVars = ttZ::performKinematicReconstruction(event, "JECDown", fitter);
                auto fourVars = ttZ::performFourLeptonsComputation(event);

                histograms.SetValues(lepVars, jetVars, recVars, fourVars);
                histograms.Fill(weight, iVariant, passed_selection, 1);

            }

            //fill JEC up histograms
            passed_selection = ttZ::passSelectionTTZ( event, "JECUp" );
            if( passed_selection >= 0 ){

                auto lepVars = ttZ::computeLeptonVariables(event);
                auto jetVars = ttZ::computeJetVariables(event, "JECUp");
                auto recVars = ttZ::performKinematicReconstruction(event, "JECUp", fitter);
                auto fourVars = ttZ::performFourLeptonsComputation(event);

                histograms.SetValues(lepVars, jetVars, recVars, fourVars);
                histograms.Fill(weight, iVariant, passed_selection, 2);

            }

            //fill JER down histograms
            passed_selection = ttZ::passSelectionTTZ( event, "JERDown" );
            if( passed_selection >= 0 ){

                auto lepVars = ttZ::computeLeptonVariables(event);
                auto jetVars = ttZ::computeJetVariables(event, "JERDown");
                auto recVars = ttZ::performKinematicReconstruction(event, "JERDown", fitter);
                auto fourVars = ttZ::performFourLeptonsComputation(event);

                histograms.SetValues(lepVars, jetVars, recVars, fourVars);
                histograms.Fill(weight, iVariant, passed_selection, 3);

            }

            //fill JER up histograms
            passed_selection = ttZ::passSelectionTTZ( event, "JERUp" );
            if( passed_selection >= 0 ){

                auto lepVars = ttZ::computeLeptonVariables(event);
                auto jetVars = ttZ::computeJetVariables(event, "JERUp");
                auto recVars = ttZ::performKinematicReconstruction(event, "JERUp", fitter);
                auto fourVars = ttZ::performFourLeptonsComputation(event);

                histograms.SetValues(lepVars, jetVars, recVars, fourVars);
                histograms.Fill(weight, iVariant, passed_selection, 4);

            }

            //fill unclustered down histograms
            passed_selection = ttZ::passSelectionTTZ( event, "UnclDown" );
            if( passed_selection >= 0 ){

                auto lepVars = ttZ::computeLeptonVariables(event);
                auto jetVars = ttZ::computeJetVariables(event, "UnclDown");
                auto recVars = ttZ::performKinematicReconstruction(event, "UnclDown", fitter);
                auto fourVars = ttZ::performFourLeptonsComputation(event);

                histograms.SetValues(lepVars, jetVars, recVars, fourVars);
                histograms.Fill(weight, iVariant, passed_selection, 5);

            }

            //fill unclustered up histograms
            passed_selection = ttZ::passSelectionTTZ( event, "UnclUp" );
            if( passed_selection >= 0 ){

                auto lepVars = ttZ::computeLeptonVariables(event);
                auto jetVars = ttZ::computeJetVariables(event, "UnclUp");
                auto recVars = ttZ::performKinematicReconstruction(event, "UnclUp", fitter);
                auto fourVars = ttZ::performFourLeptonsComputation(event);

                histograms.SetValues(lepVars, jetVars, recVars, fourVars);
                histograms.Fill(weight, iVariant, passed_selection, 6);

            }

        }
    }

    // prepare root file for datacards / combine
    std::cout << "Write output" << std::endl;
    const std::string directoryName = stringTools::formatDirectoryName( "output/"+subfolder );
    systemTools::makeDirectory( directoryName );
    const std::string fileName = procName + "_" + year + ".root";
    TFile *file = TFile::Open( (directoryName+fileName).c_str(), "RECREATE");
    histograms.Write(file);
    file->Close();
    std::cout << "Done writing output" << std::endl;

    // free kinematic reconstruction
    std::cout << "Delete kinematic reconstruction" << std::endl;
    delete fitter;
    std::cout << "Done deleting kinematic reconstruction" << std::endl;
}



int main( int argc, char* argv[] ){
    const std::string sampleDirectoryPath = "/user/joknolle/ntuples/ntuples_ttz_feb22/";
    std::vector< std::string > argvStr( &argv[0], &argv[0] + argc );

    // create lists of allowed argument values. std::set would allow for faster search,
    // but it's irrelevant here and vectors will be more convenient to handle.
    // std::set< std::string > years; years.insert("2016"); years.insert("2017"); years.insert("2018");
    std::vector< std::string > years{"2016", "2017", "2018"};
    std::vector< std::string > processes{
        "data",
        "ttZ1", "ttZ2", "ttZ3", "ttZ4", "ttZ5",
        "ttH", "tZq", "tWZ", "ttW", "tHQ", "tHW", "ttZlight", "tttt", "ttWW", "ttWZ", "ttZZ",
        "WZ3L", "WZ2L",
        "DY", "tG", "ttG", "WG", "tt",
        "ZZ4L", "ZZ2E2M", "ZZ2E2T", "ZZ2M2T", "ZZ4E", "ZZ4M", "ggHZZ", "VBFHZZ", "WpHZZ", "WmHZZ", "ZHZZ",
        "WZG", "ZZZ", "WZZ", "WWZ", "WWW", "WWDS", "WW",
    };

    // if not enough arguments, complain. Last argument is optional
    if(argc < 3){
        std::cerr << "please specify input arguments, usage:" << std::endl;
        std::cerr << "./ttZAnalysis year process" << std::endl;
        return 1;
    }
    std::string year = argvStr[1];
    std::string procName = argvStr[2];
    std::string subfolder = argc>=4 ? argvStr[3]+"/" : "";

    if ( std::find(std::begin(years),     std::end(years),     year)     == std::end(years) ||
         std::find(std::begin(processes), std::end(processes), procName) == std::end(processes)
    ){
        std::cerr << "At least one of arguments not recognized" << std::endl;
        std::cerr << "provided arguments: " << year << " " << procName << std::endl;
        std::cerr << "allowed arguments:" << std::endl;
        std::cerr << "year: 2016, 2017, 2018" << std::endl;
        std::cerr << "process to analyze:";
        for(auto const& process: processes) std::cerr << " " << process;
        std::cerr << std::endl;
        return 1;
    }

    try {
        analyze( year, sampleDirectoryPath, procName, subfolder );
        return 0;
    } catch (const std::exception& e) {
        std::cerr << e.what() << std::endl
                  << "analyze(" << year << ", " << sampleDirectoryPath << ", " << procName << ", " << subfolder << ")"
                  << " ended with error" << std::endl;
        return 1;
    }

}
