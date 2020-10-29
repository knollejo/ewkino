//include c++ library classes 
#include <memory>
#include <thread>

//include ROOT classes
#include "TH2D.h"

//include other parts of framework
#include "../TreeReader/interface/TreeReader.h"
#include "../Event/interface/Event.h"
#include "../Tools/interface/HistInfo.h"
#include "interface/fakeRateTools.h"
#include "../plotting/plotCode.h"
#include "../plotting/tdrStyle.h"
#include "../Tools/interface/systemTools.h"
#include "../Tools/interface/stringTools.h"
#include "../Tools/interface/analysisTools.h"


std::vector< HistInfo > makeDistributionInfo(){
    std::vector< HistInfo > histInfoVec = {
        HistInfo( "leptonPtLeading", "p_{T}^{leading lepton} (GeV)", 10, 25, 200 ),
        HistInfo( "leptonPtSubLeading", "p_{T}^{subleading lepton} (GeV)", 10, 15, 150 ),
        HistInfo( "leptonPtTrailing", "P_{T}^{trailing lepton} (GeV)", 15, 10, 150 ),

        HistInfo( "leptonEtaLeading", "|#eta|^{leading lepton}", 10, 0, 2.5 ),
        HistInfo( "leptonEtaSubLeading", "|#eta|^{subleading lepton}", 10, 0, 2.5 ),
        HistInfo( "leptonEtaTrailing", "|#eta|^{trailing lepton}", 10, 0, 2.5 ),
        
        HistInfo( "met", "E_{T}^{miss} (GeV)", 10, 0, 300 ),
        HistInfo( "mt", "M_{T}^{W} (GeV)", 10, 0, 300 ),
        HistInfo( "mll", "M_{ll} (GeV)", 10, 0, 200 ),
        HistInfo( "ltmet", "L_{T} + E_{T}^{miss} (GeV)", 10, 0, 300 ),
        HistInfo( "ht", "H_{T} (GeV)", 10, 0, 600 ),
        HistInfo( "m3l", "M_{3l} (GeV)", 10, 0, 300 ), 
        HistInfo( "mt3l", "M_{T}^{3l} (GeV)", 10, 0, 300 ),

        HistInfo( "nJets", "number of jets", 8, 0, 8 ),
        HistInfo( "nBJets", "number of b-jets (medium deep CSV)", 4, 0, 4 ),
        HistInfo( "nVertex", "number of vertices", 10, 0, 70 )
    };
    return histInfoVec;
}


std::shared_ptr< TH2D > readFRMap( const std::string& flavor, const std::string& year ){
    TFile* frFile = TFile::Open( ("fakeRateMaps/fakeRateMap_MC_" + flavor + "_" + year + ".root" ).c_str() );
    std::shared_ptr< TH2D > frMap( dynamic_cast< TH2D* >( frFile->Get( ( "fakeRate_" + flavor + "_" + year ).c_str() ) ) );
    frMap->SetDirectory( gROOT );
    frFile->Close();
    return frMap;
} 


bool passClosureTestEventSelection( Event& event, const bool requireMuon = false, const bool requireElectron = false ){
    event.removeTaus();
    event.applyLeptonConeCorrection();
    event.cleanElectronsFromLooseMuons();
    event.selectFOLeptons();
    if( event.numberOfLightLeptons() != 3 ) return false;
    event.selectGoodJets();
    event.cleanJetsFromFOLeptons();
    event.sortLeptonsByPt();
    
    if( !( event.passMetFilters() ) ) return false;
    size_t numberOfNonPromptLeptons = 0;
    size_t numberOfPromptLeptons = 0;
    for( auto& leptonPtr : event.lightLeptonCollection() ){
        //if( !leptonPtr->isPrompt() && leptonPtr->matchPdgId() != 22 ){
        if( !leptonPtr->isPrompt() ){
            if( requireMuon && !leptonPtr->isMuon() ) continue;
            if( requireElectron && !leptonPtr->isElectron() ) continue;
            ++numberOfNonPromptLeptons;
        }
        //else if( leptonPtr->isPrompt() && leptonPtr->matchPdgId() != 22 ) ++numberOfPromptLeptons;
        else if( leptonPtr->isPrompt() ) ++numberOfPromptLeptons;
    }
    if( numberOfNonPromptLeptons < 1 ) return false;
    if( ( numberOfNonPromptLeptons + numberOfPromptLeptons ) != 3 ) return false;
    return true;
}


double fakeRateWeight( const Event& event, const std::shared_ptr< TH2D >& frMap_muon,  const std::shared_ptr< TH2D >& frMap_electron ){
    double weight = -1.;
    for( const auto& leptonPtr : event.lightLeptonCollection() ){
        if( leptonPtr->isFO() && !leptonPtr->isTight() ){

            double croppedPt = std::min( leptonPtr->pt(), 99. );
            double croppedAbsEta = std::min( leptonPtr->absEta(), (leptonPtr->isMuon() ? 2.4 : 2.5) );

            double fr;
            if( leptonPtr->isMuon() ){
                fr = frMap_muon->GetBinContent( frMap_muon->FindBin( croppedPt, croppedAbsEta ) );
            } else {
                fr = frMap_electron->GetBinContent( frMap_electron->FindBin( croppedPt, croppedAbsEta ) );
            }
            weight *= ( - fr / ( 1. - fr ) );
        }
    }
    return weight;
}


void closureTest_MC( const std::string& process, const std::string& year, const std::string& sampleDirectory, const std::string& flavor = ""){

    analysisTools::checkYearString( year );

    //check process string
    if( ! (process == "TT" || process == "DY" ) ){
        throw std::invalid_argument( "Given closure test process argument is '" + process + "' while it should be DY or TT." );
    }

    bool onlyMuonFakes = ( flavor == "muon" );
    bool onlyElectronFakes = ( flavor == "electron" );

    //make collection of histograms
    std::vector< std::shared_ptr< TH1D > > observedHists; 
    std::vector< std::shared_ptr< TH1D > > predictedHists;

    std::vector< HistInfo > histInfoVec = makeDistributionInfo();

    for( const auto& histInfo : histInfoVec ){
        observedHists.push_back( histInfo.makeHist( histInfo.name() + "_observed_" + process + "_" + year ) );
        predictedHists.push_back( histInfo.makeHist( histInfo.name() + "_predicted_"  + process + "_" + year ) );
    }
    
    //read fake-rate map corresponding to this year and flavor 
    std::shared_ptr< TH2D > fakeRateMap_muon = readFRMap( "muon", year );
    std::shared_ptr< TH2D > fakeRateMap_electron = readFRMap( "electron", year );


    //loop over samples to fill histograms
    std::string sampleListFile = "sampleLists/samples_closureTest_" + process + "_" + year + ".txt";
    TreeReader treeReader( sampleListFile, sampleDirectory );
    for( unsigned i = 0; i < treeReader.numberOfSamples(); ++i ){
        treeReader.initSample();
    
        for( long unsigned entry = 0; entry < treeReader.numberOfEntries(); ++entry ){
            Event event = treeReader.buildEvent( entry );

//            if(entry > 1000) break;
            //apply event selection
            //if (flavor != "electron")  continue;
            if( !passClosureTestEventSelection( event, onlyMuonFakes, onlyElectronFakes ) ) continue;
            //light lepton collection
            LightLeptonCollection lightLeptons = event.lightLeptonCollection();

            double mll = 0, mtW = 0;
            if( event.hasOSSFLightLeptonPair() ){
                mll = event.bestZBosonCandidateMass();
                mtW = event.mtW();
            } else{
                mll = ( lightLeptons[0] + lightLeptons[1] ).mass();
                mtW = mt( lightLeptons[2], event.met() );
            }
            //compute plotting variables 
            std::vector< double > variables = { lightLeptons[0].pt(), lightLeptons[1].pt(), lightLeptons[2].pt(),
                lightLeptons[0].absEta(), lightLeptons[1].absEta(), lightLeptons[2].absEta(),
                event.metPt(),
                mtW,
                mll,
                lightLeptons.scalarPtSum() + event.metPt(),
                event.HT(),
                lightLeptons.mass(),
                mt( lightLeptons.objectSum(), event.met() ),
                static_cast< double >( event.numberOfJets() ),
                static_cast< double >( event.numberOfMediumBTaggedJets() ),
                static_cast< double >( event.numberOfVertices() )
            };
                
            //event is 'observed' if all leptons are tight 
            //bool isObserved = ( event.numberOfTightLeptons() == event.numberOfLightLeptons() );
            bool isObserved = true;
            for( const auto& leptonPtr : lightLeptons ){
                if( !leptonPtr->isTight() ){
                    isObserved = false;
                }
            }

            if( isObserved ){
//            std::cout << "passed 1, flavor: " << flavor << " weight: " << event.weight() << std::endl;
                for( std::vector< double >::size_type v = 0; v < variables.size(); ++v ){
                    observedHists[v]->Fill( std::min( variables[v],  histInfoVec[v].maxBinCenter() ), event.weight() );
                }

            } else {

                //compute event weight with fake-rate
                double weight = event.weight()*fakeRateWeight( event, fakeRateMap_muon, fakeRateMap_electron );
//            std::cout << "passed 1, flavor: " << flavor << " weight: " << weight << std::endl;
                for( std::vector< double >::size_type v = 0; v < variables.size(); ++v ){
                    predictedHists[v]->Fill( std::min( variables[v],  histInfoVec[v].maxBinCenter() ), weight );
                }
            }
        }
    }

    //make plot output directory
    std::string outputDirectory_name = "./closurePlots_MC_" + process + "_" + year; 
    if( flavor != "" ){
        outputDirectory_name += ( "_" + flavor );
    }
    systemTools::makeDirectory( outputDirectory_name );
    
    //make plots 
    for( std::vector< HistInfo >::size_type v = 0; v < histInfoVec.size(); ++v ){
        std::string names[2] = {"MC observed", "fake-rate prediction"};
        std::vector< TH1D* > predicted = { predictedHists[v].get() };
        std::string header;
        if( year == "2016" ){
            header = "35.9 fb^{-1}";
        } else if( year == "2017" ){
            header = "41.5 fb^{-1}";
        } else{
            header = "59.7 fb^{-1}";
        }
        TH1D* systUnc = dynamic_cast< TH1D* >( predictedHists[v]->Clone() );
        for( int b = 1; b < systUnc->GetNbinsX() + 1; ++b ){
            systUnc->SetBinContent( b , systUnc->GetBinContent(b)*0.3 );
        }
        if( flavor == "" ){
               plotDataVSMC( observedHists[v].get(), &predicted[0], names, 1, stringTools::formatDirectoryName( outputDirectory_name ) + histInfoVec[v].name() + "_" + process + "_" + year + ".pdf", "", false, false, header, systUnc);
        } else {
               plotDataVSMC( observedHists[v].get(), &predicted[0], names, 1, stringTools::formatDirectoryName( outputDirectory_name ) + histInfoVec[v].name() + "_" + process + "_" + year + "_" + flavor + ".pdf", "", false, false, header, systUnc);
        }
    }
}


int main( int argc, char* argv[] ){

    //const std::string ntupleDirectory = "/pnfs/iihe/cms/store/user/wverbeke/ntuples_ewkino";
    //const std::string ntupleDirectory = "~/Work/ntuples_ewkino_new/";
    const std::string ntupleDirectory = "/user/mniedzie/Work/ntuples_ttz_new/";

    setTDRStyle();

    std::vector< std::string > argvStr( &argv[0], &argv[0] + argc );

    // create lists of allowed argument values. std::set would allow for faster search, 
    // but it's irrelevant here and vectors will be more convenient to handle.
    // std::set< std::string > years; years.insert("2016"); years.insert("2017"); years.insert("2018");
    std::vector< std::string > years{"2016", "2017", "2018"};
    std::vector< std::string > processes{"DY", "TT"};

    // if not enough arguments, complain. Last argument is optional
    if(argc < 3){
        std::cerr << "please specify input arguments, usage:" << std::endl;
        std::cerr << "./closureTest_MC year region " << std::endl;
        return 1;
    }    
    std::string year = argvStr[1];
    std::string process = argvStr[2];

    if ( std::find(std::begin(years),     std::end(years),     year)          == std::end(years) || 
         std::find(std::begin(processes),  std::end(processes),  process) == std::end(processes)  
    ){
        std::cerr << "At least one of arguments not recognized" << std::endl;
        std::cerr << "allowed arguments:" << std::endl;
        std::cerr << "year:    2016, 2017, 2018"                      << std::endl;
        std::cerr << "process: DY, TT"    << std::endl;
        return 1;
    }

       //make sure ROOT behaves itself when running multithreaded
    ROOT::EnableThreadSafety();

//    closureTest_MC( process, year, ntupleDirectory, "" );
    //closureTest_MC( process, year, ntupleDirectory, "muon" );
    closureTest_MC( process, year, ntupleDirectory, "electron" );
    //make threads to run different closure tests 
    //closureTest_MC( "DY", "2016", ntupleDirectory, "" );
    //closureTest_MC( "DY", "2017", ntupleDirectory, "" );
    //closureTest_MC( "DY", "2018", ntupleDirectory, "" );
    //closureTest_MC( "TT", "2016", ntupleDirectory, "" );
    //closureTest_MC( "TT", "2017", ntupleDirectory, "" );
    //closureTest_MC( "TT", "2018", ntupleDirectory, "" );

    //closureTest_MC( "DY", "2016", ntupleDirectory, "muon" );
    //closureTest_MC( "DY", "2017", ntupleDirectory, "muon" );
    //closureTest_MC( "DY", "2018", ntupleDirectory, "muon" );
    //closureTest_MC( "TT", "2016", ntupleDirectory, "muon" );
    //closureTest_MC( "TT", "2017", ntupleDirectory, "muon" );
    //closureTest_MC( "TT", "2018", ntupleDirectory, "muon" );

    //closureTest_MC( "DY", "2016", ntupleDirectory, "electron" );
    //closureTest_MC( "DY", "2017", ntupleDirectory, "electron" );
    //closureTest_MC( "DY", "2018", ntupleDirectory, "electron" );
    //closureTest_MC( "TT", "2016", ntupleDirectory, "electron" );
    //closureTest_MC( "TT", "2017", ntupleDirectory, "electron" );
    //closureTest_MC( "TT", "2018", ntupleDirectory, "electron" );
    //std::vector< std::thread > threadVector;
    //threadVector.reserve( 6 );
    //for( const auto& process : {"TT", "DY" } ){
    //    for( const auto& year : {"2016", "2017", "2018" } ){
    //        threadVector.emplace_back( closureTest_MC, process, year, "~/Work/ntuples_ewkino_new/", "" );
    //    }
    //}

    //for( auto& t : threadVector ){
    //    t.join();
    //}

    /*
    std::vector< std::thread > threadVector_muon;
    threadVector_muon.reserve( 6 );
    for( const auto& process : {"TT", "DY" } ){
        for( const auto& year : {"2016", "2017", "2018" } ){
            threadVector_muon.emplace_back( closureTest_MC, process, year, "~/Work/ntuples_ewkino_new/", "muon" );
        }
    }

    for( auto& t : threadVector_muon ){
        t.join();
    }

    std::vector< std::thread > threadVector_electron;
    threadVector_electron.reserve( 6 );
    for( const auto& process : {"TT", "DY" } ){
        for( const auto& year : {"2016", "2017", "2018" } ){
            threadVector_electron.emplace_back( closureTest_MC, process, year, "~/Work/ntuples_ewkino_new/", "electron" );
        }
    }

    for( auto& t : threadVector_electron ){
        t.join();
    }
    */
 
    return 0;
}
