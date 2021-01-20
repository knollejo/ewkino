#ifndef ttZVariables_H
#define ttZVariables_H

//include c++ library classes
#include <string>
#include <map>

//include other parts of framework
#include "../../Event/interface/Event.h"
#include "../interface/kinFitter.h"

namespace ttZ {
    std::map<std::string, double> computeLeptonVariables(Event& event);
    std::map<std::string, double> computeJetVariables(Event& event, const std::string& unc);
    std::map<std::string, double> performKinematicReconstruction(Event& event, const std::string& unc, KinFitter* fitter);
}

#endif
