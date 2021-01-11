#include "../objects/interface/JetSelector.h"

//b-tagging working points
#include "bTagWP.h"

/*
Jet ID selection
*/

bool JetSelector::isGoodBase() const{
    if( jetPtr->pt() < 30 ) return false;
    if( fabs( jetPtr->eta() ) > 2.4 ) return false;
    if( !jetPtr->isTight() ) return false;
    return true;
}


bool JetSelector::isGood2016() const{
    return true;
}


bool JetSelector::isGood2017() const{
    return true;
}


bool JetSelector::isGood2018() const{
    return true;
}


/*
b-tagging working points
*/

bool JetSelector::inBTagAcceptance() const{
    if( !isGood() ) return false;
    if( jetPtr->pt() < 30 ) return false;
    if( fabs( jetPtr->eta() ) >= 2.4 ) return false;
    return true;
}


bool JetSelector::isBTaggedLoose2016() const{
    return ( jetPtr->deepFlavor() > bTagWP::looseDeepFlavor2016()  );
}


bool JetSelector::isBTaggedLoose2017() const{
    return ( jetPtr->deepFlavor() > bTagWP::looseDeepFlavor2017() ); 
}


bool JetSelector::isBTaggedLoose2018() const{
    return ( jetPtr->deepFlavor() > bTagWP::looseDeepFlavor2018() );
}
        

bool JetSelector::isBTaggedMedium2016() const{
    return ( jetPtr->deepFlavor() > bTagWP::mediumDeepFlavor2016() );
}


bool JetSelector::isBTaggedMedium2017() const{
    return ( jetPtr->deepFlavor() > bTagWP::mediumDeepFlavor2017() );
}


bool JetSelector::isBTaggedMedium2018() const{
    return ( jetPtr->deepFlavor() > bTagWP::mediumDeepFlavor2018() );
}


bool JetSelector::isBTaggedTight2016() const{
    return ( jetPtr->deepFlavor() > bTagWP::tightDeepFlavor2016() );
}


bool JetSelector::isBTaggedTight2017() const{
    return ( jetPtr->deepFlavor() > bTagWP::tightDeepFlavor2017() );
}


bool JetSelector::isBTaggedTight2018() const{
    return ( jetPtr->deepFlavor() > bTagWP::tightDeepFlavor2018() );
}
