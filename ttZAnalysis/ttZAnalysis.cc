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
#include "interface/ttZVariables.h"

#include <functional>


void analyze( const std::string& year, const std::string& sampleDirectoryPath , const std::string& procName , const std::string& subfolder ){

    //check setup
    analysisTools::checkYearString( year );
    const bool is_signal = (procName == "ttZ");
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
    const int nVariants = is_signal ? 4 : 2;
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
            if(entry>1000) break;
            Event event = treeReader.buildEvent( entry );

            //check if sample has pdf and scale information stored, try-catch blocks for every event are too slow
            if( entry == 0 && treeReader.isMC() ){
                try{
                    event.generatorInfo().relativeWeightPdfVar( 99 );
                }catch( std::out_of_range& ){
                    sampleHasPdfAndScale = false;
                }
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
            if( !ttZ::leptonsAreTight( event ) ){
                if( event.isMC() ) continue;
                is_prompt = false;
                weight *= ttZ::fakeRateWeight( event, frMapMuons, frMapElectrons );
            }

            //fill nominal histograms
            int passed_selection = ttZ::passSelectionTTZ( event, "nominal" );
            if( passed_selection ){

                auto lepVars = ttZ::computeLeptonVariables(event);
                auto jetVars = ttZ::computeJetVariables(event, "nominal");
                auto recVars = ttZ::performKinematicReconstruction(event, "nominal", fitter);

                histograms.SetValues(lepVars, jetVars, recVars);
                histograms.Fill(weight, is_prompt ? 0 : 1, 0, 0);
                histograms.Fill(weight, is_prompt ? 0 : 1, passed_selection==3 ? 1 : 2, 0);

                //no uncertainties for data
                if( event.isData() ) continue;

                //fill scale variation histograms
                const double weightScaleRenDown = sampleHasPdfAndScale ? event.generatorInfo().relativeWeight_MuR_0p5_MuF_1() : 1.0;
                histograms.Fill(weight*weightScaleRenDown, is_prompt ? 0 : 1, 0, 7);
                histograms.Fill(weight*weightScaleRenDown, is_prompt ? 0 : 1, passed_selection==3 ? 1 : 2, 7);
                const double weightScaleRenUp = sampleHasPdfAndScale ? event.generatorInfo().relativeWeight_MuR_2_MuF_1() : 1.0;
                histograms.Fill(weight*weightScaleRenUp, is_prompt ? 0 : 1, 0, 8);
                histograms.Fill(weight*weightScaleRenUp, is_prompt ? 0 : 1, passed_selection==3 ? 1 : 2, 8);
                const double weightScaleFacDown = sampleHasPdfAndScale ? event.generatorInfo().relativeWeight_MuR_1_MuF_0p5() : 1.0;
                histograms.Fill(weight*weightScaleFacDown, is_prompt ? 0 : 1, 0, 9);
                histograms.Fill(weight*weightScaleFacDown, is_prompt ? 0 : 1, passed_selection==3 ? 1 : 2, 9);
                const double weightScaleFacUp = sampleHasPdfAndScale ? event.generatorInfo().relativeWeight_MuR_1_MuF_2() : 1.0;
                histograms.Fill(weight*weightScaleFacUp, is_prompt ? 0 : 1, 0, 10);
                histograms.Fill(weight*weightScaleFacUp, is_prompt ? 0 : 1, passed_selection==3 ? 1 : 2, 10);
                const double weightScaleBothDown = sampleHasPdfAndScale ? event.generatorInfo().relativeWeight_MuR_0p5_MuF_0p5() : 1.0;
                histograms.Fill(weight*weightScaleBothDown, is_prompt ? 0 : 1, 0, 11);
                histograms.Fill(weight*weightScaleBothDown, is_prompt ? 0 : 1, passed_selection==3 ? 1 : 2, 11);
                const double weightScaleBothUp = sampleHasPdfAndScale ? event.generatorInfo().relativeWeight_MuR_2_MuF_2() : 1.0;
                histograms.Fill(weight*weightScaleBothUp, is_prompt ? 0 : 1, 0, 12);
                histograms.Fill(weight*weightScaleBothUp, is_prompt ? 0 : 1, passed_selection==3 ? 1 : 2, 12);

                //fill pdf histograms
                for( unsigned pdf_i = 0; pdf_i < nPdfVariations; ++pdf_i){
                    const double weightPdf = sampleHasPdfAndScale ? event.generatorInfo().relativeWeightPdfVar(pdf_i) : 1.;
                    histograms.Fill(weight*weightPdf, is_prompt ? 0 : 1, 0, 27+pdf_i);
                    histograms.Fill(weight*weightPdf, is_prompt ? 0 : 1, passed_selection==3 ? 1 : 2, 27+pdf_i);
                }

                //fill pileup histograms
                double weightPileupDown = reweighter[ "pileup" ]->weightDown( event ) / reweighter[ "pileup" ]->weight( event );
                if ( std::isnan(weightPileupDown) ) weightPileupDown = 0.;
                histograms.Fill(weight*weightPileupDown, is_prompt ? 0 : 1, 0, 13);
                histograms.Fill(weight*weightPileupDown, is_prompt ? 0 : 1, passed_selection==3 ? 1 : 2, 13);
                double weightPileupUp = reweighter[ "pileup" ]->weightUp( event ) / reweighter[ "pileup" ]->weight( event );
                if ( std::isnan(weightPileupUp) ) weightPileupUp = 0.;
                histograms.Fill(weight*weightPileupUp, is_prompt ? 0 : 1, 0, 14);
                histograms.Fill(weight*weightPileupUp, is_prompt ? 0 : 1, passed_selection==3 ? 1 : 2, 14);

                //fill b-tag histograms
                const double weightBTagHeavyDown = reweighter[ "bTag_heavy" ]->weightDown( event ) / reweighter[ "bTag_heavy" ]->weight( event );
                histograms.Fill(weight*weightBTagHeavyDown, is_prompt ? 0 : 1, 0, 15);
                histograms.Fill(weight*weightBTagHeavyDown, is_prompt ? 0 : 1, passed_selection==3 ? 1 : 2, 15);
                const double weightBTagHeavyUp = reweighter[ "bTag_heavy" ]->weightUp( event ) / reweighter[ "bTag_heavy" ]->weight( event );
                histograms.Fill(weight*weightBTagHeavyUp, is_prompt ? 0 : 1, 0, 16);
                histograms.Fill(weight*weightBTagHeavyUp, is_prompt ? 0 : 1, passed_selection==3 ? 1 : 2, 16);
                const double weightBTagLightDown = reweighter[ "bTag_light" ]->weightDown( event ) / reweighter[ "bTag_light" ]->weight( event );
                histograms.Fill(weight*weightBTagLightDown, is_prompt ? 0 : 1, 0, 17);
                histograms.Fill(weight*weightBTagLightDown, is_prompt ? 0 : 1, passed_selection==3 ? 1 : 2, 17);
                const double weightBTagLightUp = reweighter[ "bTag_light" ]->weightUp( event ) / reweighter[ "bTag_light" ]->weight( event );
                histograms.Fill(weight*weightBTagLightUp, is_prompt ? 0 : 1, 0, 18);
                histograms.Fill(weight*weightBTagLightUp, is_prompt ? 0 : 1, passed_selection==3 ? 1 : 2, 18);

                //fill prefiring histograms
                const double weightPrefireDown = reweighter[ "prefire" ]->weightDown( event ) / reweighter[ "prefire" ]->weight( event );
                histograms.Fill(weight*weightPrefireDown, is_prompt ? 0 : 1, 0, 19);
                histograms.Fill(weight*weightPrefireDown, is_prompt ? 0 : 1, passed_selection==3 ? 1 : 2, 19);
                const double weightPrefireUp = reweighter[ "prefire" ]->weightUp( event ) / reweighter[ "prefire" ]->weight( event );
                histograms.Fill(weight*weightPrefireUp, is_prompt ? 0 : 1, 0, 20);
                histograms.Fill(weight*weightPrefireUp, is_prompt ? 0 : 1, passed_selection==3 ? 1 : 2, 20);

                //fill lepton histograms
                const double leptonIDWeightDownMuon = reweighter[ "muonID" ]->weightDown( event ) / reweighter[ "muonID" ]->weight( event );
                histograms.Fill(weight*leptonIDWeightDownMuon, is_prompt ? 0 : 1, 0, 21);
                histograms.Fill(weight*leptonIDWeightDownMuon, is_prompt ? 0 : 1, passed_selection==3 ? 1 : 2, 21);
                const double leptonIDWeightUpMuon = reweighter[ "muonID" ]->weightUp( event ) / reweighter[ "muonID" ]->weight( event );
                histograms.Fill(weight*leptonIDWeightUpMuon, is_prompt ? 0 : 1, 0, 22);
                histograms.Fill(weight*leptonIDWeightUpMuon, is_prompt ? 0 : 1, passed_selection==3 ? 1 : 2, 22);
                const double recoWeightDownElec = reweighter[ "electronReco_pTBelow20" ]->weightDown( event ) * reweighter[ "electronReco_pTAbove20" ]->weightDown( event ) / ( reweighter[ "electronReco_pTBelow20" ]->weight( event ) * reweighter[ "electronReco_pTAbove20" ]->weight( event ) );
                histograms.Fill(weight*recoWeightDownElec, is_prompt ? 0 : 1, 0, 23);
                histograms.Fill(weight*recoWeightDownElec, is_prompt ? 0 : 1, passed_selection==3 ? 1 : 2, 23);
                const double recoWeightUpElec = reweighter[ "electronReco_pTBelow20" ]->weightUp( event ) * reweighter[ "electronReco_pTAbove20" ]->weightUp( event ) / ( reweighter[ "electronReco_pTBelow20" ]->weight( event ) * reweighter[ "electronReco_pTAbove20" ]->weight( event ) );
                histograms.Fill(weight*recoWeightUpElec, is_prompt ? 0 : 1, 0, 24);
                histograms.Fill(weight*recoWeightUpElec, is_prompt ? 0 : 1, passed_selection==3 ? 1 : 2, 24);
                const double leptonIDWeightDownElec = reweighter[ "electronID" ]->weightDown( event ) / reweighter[ "electronID" ]->weight( event );
                histograms.Fill(weight*leptonIDWeightDownElec, is_prompt ? 0 : 1, 0, 25);
                histograms.Fill(weight*leptonIDWeightDownElec, is_prompt ? 0 : 1, passed_selection==3 ? 1 : 2, 25);
                const double leptonIDWeightUpElec = reweighter[ "electronID" ]->weightUp( event ) / reweighter[ "electronID" ]->weight( event );
                histograms.Fill(weight*leptonIDWeightUpElec, is_prompt ? 0 : 1, 0, 26);
                histograms.Fill(weight*leptonIDWeightUpElec, is_prompt ? 0 : 1, passed_selection==3 ? 1 : 2, 26);

            }

            //no uncertainties for data
            if( event.isData() ) continue;

            //fill JEC down histograms
            passed_selection = ttZ::passSelectionTTZ( event, "JECDown" );
            if( passed_selection ){

                auto lepVars = ttZ::computeLeptonVariables(event);
                auto jetVars = ttZ::computeJetVariables(event, "JECDown");
                auto recVars = ttZ::performKinematicReconstruction(event, "JECDown", fitter);

                histograms.SetValues(lepVars, jetVars, recVars);
                histograms.Fill(weight, is_prompt ? 0 : 1, 0, 1);
                histograms.Fill(weight, is_prompt ? 0 : 1, passed_selection==3 ? 1 : 2, 1);

            }

            //fill JEC up histograms
            passed_selection = ttZ::passSelectionTTZ( event, "JECUp" );
            if( passed_selection ){

                auto lepVars = ttZ::computeLeptonVariables(event);
                auto jetVars = ttZ::computeJetVariables(event, "JECUp");
                auto recVars = ttZ::performKinematicReconstruction(event, "JECUp", fitter);

                histograms.SetValues(lepVars, jetVars, recVars);
                histograms.Fill(weight, is_prompt ? 0 : 1, 0, 2);
                histograms.Fill(weight, is_prompt ? 0 : 1, passed_selection==3 ? 1 : 2, 2);

            }

            //fill JER down histograms
            passed_selection = ttZ::passSelectionTTZ( event, "JERDown" );
            if( passed_selection ){

                auto lepVars = ttZ::computeLeptonVariables(event);
                auto jetVars = ttZ::computeJetVariables(event, "JERDown");
                auto recVars = ttZ::performKinematicReconstruction(event, "JERDown", fitter);

                histograms.SetValues(lepVars, jetVars, recVars);
                histograms.Fill(weight, is_prompt ? 0 : 1, 0, 3);
                histograms.Fill(weight, is_prompt ? 0 : 1, passed_selection==3 ? 1 : 2, 3);

            }

            //fill JER up histograms
            passed_selection = ttZ::passSelectionTTZ( event, "JERUp" );
            if( passed_selection ){

                auto lepVars = ttZ::computeLeptonVariables(event);
                auto jetVars = ttZ::computeJetVariables(event, "JERUp");
                auto recVars = ttZ::performKinematicReconstruction(event, "JERUp", fitter);

                histograms.SetValues(lepVars, jetVars, recVars);
                histograms.Fill(weight, is_prompt ? 0 : 1, 0, 4);
                histograms.Fill(weight, is_prompt ? 0 : 1, passed_selection==3 ? 1 : 2, 4);

            }

            //fill unclustered down histograms
            passed_selection = ttZ::passSelectionTTZ( event, "UnclDown" );
            if( passed_selection ){

                auto lepVars = ttZ::computeLeptonVariables(event);
                auto jetVars = ttZ::computeJetVariables(event, "UnclDown");
                auto recVars = ttZ::performKinematicReconstruction(event, "UnclDown", fitter);

                histograms.SetValues(lepVars, jetVars, recVars);
                histograms.Fill(weight, is_prompt ? 0 : 1, 0, 5);
                histograms.Fill(weight, is_prompt ? 0 : 1, passed_selection==3 ? 1 : 2, 5);

            }

            //fill unclustered up histograms
            passed_selection = ttZ::passSelectionTTZ( event, "UnclUp" );
            if( passed_selection ){

                auto lepVars = ttZ::computeLeptonVariables(event);
                auto jetVars = ttZ::computeJetVariables(event, "UnclUp");
                auto recVars = ttZ::performKinematicReconstruction(event, "UnclUp", fitter);

                histograms.SetValues(lepVars, jetVars, recVars);
                histograms.Fill(weight, is_prompt ? 0 : 1, 0, 6);
                histograms.Fill(weight, is_prompt ? 0 : 1, passed_selection==3 ? 1 : 2, 6);

            }

        }
    }

    // prepare root file for datacards / combine
    std::cout << "Write output" << std::endl;
    const std::string directoryName = stringTools::formatDirectoryName( "output/"+subfolder );
    systemTools::makeDirectory( directoryName );
    const std::string fileName = procName + "_" + year + ".root";
    TFile *file = TFile::Open( (directoryName+fileName).c_str(), "RECREATE");
    histograms.Write();
    file->Close();
    std::cout << "Done writing output" << std::endl;

    // free kinematic reconstruction
    std::cout << "Delete kinematic reconstruction" << std::endl;
    delete fitter;
    std::cout << "Done deleting kinematic reconstruction" << std::endl;
}



int main( int argc, char* argv[] ){
    const std::string sampleDirectoryPath = "/user/mniedzie/Work/ntuples_ttz_new/";
    std::vector< std::string > argvStr( &argv[0], &argv[0] + argc );

    // create lists of allowed argument values. std::set would allow for faster search,
    // but it's irrelevant here and vectors will be more convenient to handle.
    // std::set< std::string > years; years.insert("2016"); years.insert("2017"); years.insert("2018");
    std::vector< std::string > years{"2016", "2017", "2018"};
    std::vector< std::string > processes{
        "data",
        "ttZ",
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
    std::string subfolder = argc>=3 ? argvStr[3]+"/" : "";

    if ( std::find(std::begin(years),     std::end(years),     year)     == std::end(years) ||
         std::find(std::begin(processes), std::end(processes), procName) == std::end(processes)
    ){
        std::cerr << "At least one of arguments not recognized" << std::endl;
        std::cerr << "provided arguments: " << year << " " << procName << std::endl;
        std::cerr << "allowed arguments:" << std::endl;
        std::cerr << "year: 2016, 2017, 2018" << std::endl;
        std::cerr << "process to analyze:"
        for(auto const& process: processes) std::cerr << " " << process;
        std::cerr << std::endl;
        return 1;
    }

    analyze( year, sampleDirectoryPath, procName, subfolder );

    return 0;
}
