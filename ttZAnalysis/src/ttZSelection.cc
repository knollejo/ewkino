#include "../interface/ttZSelection.h"

//include c++ library classes
#include <stdexcept>

//include other parts of framework
#include "../../Tools/interface/histogramTools.h"
#include "../../Tools/interface/stringTools.h"
#include "../../constants/particleMasses.h"


void ttZ::applyBaselineObjectSelection( Event& event, const bool allowUncertainties ){
    event.selectLooseLeptons();
    event.removeTaus();
    event.cleanElectronsFromLooseMuons();
    event.cleanJetsFromFOLeptons();
    if( allowUncertainties ){
        event.jetCollection().selectGoodAnyVariationJets();
    } else {
        event.selectGoodJets();
    }
}


bool ttZ::passLowMllVeto( const Event& event, const double vetoValue ){
    for( const auto& leptonPtrPair : event.lightLeptonCollection().pairCollection() ){

        //veto OSSF pairs of low mass
        Lepton& lepton1 = *( leptonPtrPair.first );
        Lepton& lepton2 = *( leptonPtrPair.second );
        if( !oppositeSignSameFlavor( lepton1, lepton2 ) ){
            continue;
        }
        if( ( lepton1 + lepton2 ).mass() < vetoValue ){
            return false;
        }
    }
    return true;
}


bool ttZ::passBaselineSelection( Event& event, const bool allowUncertainties, const bool mllVeto ){
	// only loose leptons are selected
    applyBaselineObjectSelection( event, allowUncertainties );
	// accept event with at least 3 loose leptons
    if( event.numberOfLeptons() < 3 ) return false;
    if( mllVeto && !passLowMllVeto( event, 12 ) ) return false;
    event.selectFOLeptons();
    if( event.numberOfLeptons() < 3 ) return false;
    event.applyLeptonConeCorrection();
    event.sortLeptonsByPt();
    return true;
}


JetCollection ttZ::variedJetCollection( const Event& event, const std::string& uncertainty ){
    if( uncertainty == "nominal" ){
        return event.jetCollection().goodJetCollection();
    } else if( uncertainty == "JECDown" ){
        return event.jetCollection().JECDownCollection().goodJetCollection();
    } else if( uncertainty == "JECUp" ){
        return event.jetCollection().JECUpCollection().goodJetCollection();
    } else if( uncertainty == "JERDown" ){
        return event.jetCollection().JERDownCollection().goodJetCollection();
    } else if( uncertainty == "JERUp" ){
        return event.jetCollection().JERUpCollection().goodJetCollection();
    } else if( uncertainty == "UnclDown" ){
        return event.jetCollection().goodJetCollection();
    } else if( uncertainty == "UnclUp" ){
        return event.jetCollection().goodJetCollection();
    } else {
        throw std::invalid_argument( "Uncertainty source " + uncertainty + " is unknown." );
    }
}


JetCollection::size_type ttZ::numberOfVariedBJets( const Event& event, const std::string& uncertainty ){
    return ttZ::variedJetCollection( event, uncertainty ).numberOfMediumBTaggedJets();
}

JetCollection::size_type ttZ::numberOfVariedJets( const Event& event, const std::string& uncertainty ){
    return ttZ::variedJetCollection( event, uncertainty ).numberOfGoodJets();
}


Met ttZ::variedMet( const Event& event, const std::string& uncertainty ){
    if( uncertainty == "nominal" ){
        return event.met();
    } else if( uncertainty == "JECDown" ){
        return event.met().MetJECDown();
    } else if( uncertainty == "JECUp" ){
        return event.met().MetJECUp();
    } else if( uncertainty == "JERDown" ){
        return event.met();
    } else if( uncertainty == "JERUp" ){
        return event.met();
    } else if( uncertainty == "UnclDown" ){
        return event.met().MetUnclusteredDown();
    } else if( uncertainty == "UnclUp" ){
        return event.met().MetUnclusteredUp();
    } else {
        throw std::invalid_argument( "Uncertainty source " + uncertainty + " is unknown." );
    }
}


bool ttZ::passSelectionLNumber( Event& event ){
    if ( event.numberOfFOLeptons() == 3 || event.numberOfFOLeptons() == 4 ) return true;
    return false;
}


int ttZ::passSelectionTTZ( Event& event, const std::string& uncertainty ){
// -1: failed selection
// 0: passed 4l selection
// 1: passed 3l selection with 3 jets
// 2: passed 3l selection with 4 or more jets

    if ( numberOfVariedJets( event, uncertainty ) < 2 ) return -1;
    if ( numberOfVariedBJets( event, uncertainty ) < 1 ) return -1;
    if ( !( event.hasOSSFLightLeptonPair() ) ) return -1;
    if ( event.numberOfFOLeptons() < 3 ) return -1;
    if ( event.numberOfFOLeptons() == 3 ){
        if ( std::abs( event.bestZBosonCandidateMass() - particle::mZ ) >= 10 ) return -1;
        if ( numberOfVariedJets( event, uncertainty ) < 3 ) return -1;
        else if ( numberOfVariedJets( event, uncertainty ) == 3 ) return 1;
        else return 2;
    } else if ( event.numberOfFOLeptons() == 4 ) {
        if ( std::abs( event.bestZBosonCandidateMass() - particle::mZ ) >= 20 ) return -1;
        std::vector< LeptonCollection::size_type > ind2Zcand;
        for( LeptonCollection::size_type l = 0; l < event.numberOfLeptons(); ++l ) {
            if( ( l == event.bestZBosonCandidateIndices().first ) || ( l == event.bestZBosonCandidateIndices().second ) ) continue;
            ind2Zcand.push_back(l);
        }
        auto firstOtherLepton = &event.leptonCollection()[ind2Zcand.at(0)];
        auto secondOtherLepton = &event.leptonCollection()[ind2Zcand.at(1)];
        if( oppositeSignSameFlavor(*firstOtherLepton, *secondOtherLepton) ) {
            if ( std::abs( ( *firstOtherLepton + *secondOtherLepton ).mass() - particle::mZ ) < 20 ) return -1;
        }
        return 0;
    } else return -1;
}


bool ttZ::passTriggerSelection( const Event& event ){
    if( event.numberOfMuons() >= 1 ){
        if( event.passTriggers_m() ) return true;
    }
    if( event.numberOfMuons() >= 2 ){
        if( event.passTriggers_mm() ) return true;
    }
    if( event.numberOfMuons() >= 3 ){
        if( event.passTriggers_mmm() ) return true;
    }
    if( event.numberOfElectrons() >= 1 ){
        if( event.passTriggers_e() ) return true;
    }
    if( event.numberOfElectrons() >= 2 ){
        if( event.passTriggers_ee() ) return true;
    }
    if( event.numberOfElectrons() >= 3 ){
        if( event.passTriggers_eee() ) return true;
    }
    if( ( event.numberOfMuons() >= 1 ) && ( event.numberOfElectrons() >= 1 ) ){
        if( event.passTriggers_em() ) return true;
    }
    if( ( event.numberOfMuons() >= 2 ) && ( event.numberOfElectrons() >= 1 ) ){
        if( event.passTriggers_emm() ) return true;
    }
    if( ( event.numberOfMuons() >= 1 ) && ( event.numberOfElectrons() >= 2 ) ){
        if( event.passTriggers_eem() ) return true;
    }

    return false;
}


bool ttZ::passPtCuts( const Event& event ){

    //leading lepton
    if( event.lepton( 0 ).pt() <= 40 ) return false;

    //subleading lepton
    if( event.lepton( 1 ).pt() <= 20 ) return false;

    //trailing lepton
    if( event.lepton( 2 ).pt() <= 10 ) return false;

    return true;
}


bool ttZ::leptonsArePrompt( const Event& event ){
    for( const auto& leptonPtr : event.leptonCollection() ){
        if( leptonPtr->isFO() && !leptonPtr->isPrompt() ) return false;
    }
    return true;
}


bool ttZ::leptonsAreTight( const Event& event ){
    for( const auto& leptonPtr : event.leptonCollection() ){
        if( leptonPtr->isFO() && !leptonPtr->isTight() ) return false;
    }
    return true;
}


bool leptonFromMEExternalConversion( const Lepton& lepton ){
    if( !( lepton.matchPdgId() == 22 ) ) return false;
    if( !( lepton.isPrompt() && lepton.provenanceConversion() == 0 ) ) return false;
    return true;
}


double ttZ::fakeRateWeight( const Event& event, const std::shared_ptr< TH2 >& muonMap, const std::shared_ptr< TH2 >& electronMap ){
    static constexpr double maxPt = 44.;

    double weight = -1.;
    for( const auto& leptonPtr : event.leptonCollection() ){
        if( !leptonPtr->isFO() ) continue;
        if( leptonPtr->isTight() ) continue;
        double fr;
        double pt = std::min( leptonPtr->pt(), maxPt );
        if( leptonPtr->isMuon() ){
            fr = histogram::contentAtValues( muonMap.get(), pt, leptonPtr->absEta() );
        } else if( leptonPtr->isElectron() ){
            fr = histogram::contentAtValues( electronMap.get(), pt, leptonPtr->absEta() );
        } else {
            throw std::invalid_argument( "we are not considering taus for now" );
        }
        weight *= - fr / ( 1. - fr );
    }
    return weight;
}
