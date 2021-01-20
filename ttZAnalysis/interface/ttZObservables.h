#include "TMath.h"

namespace ttZObservables {

const unsigned int nBinsRec = 8;
const double pi = TMath::Pi();
namespace binsRec {
    const double ttzMass[nBinsRec+1] = { 436.0, 580.0, 640.0, 720.0, 800.0, 890.0, 1000.0, 1150.0, 1500.0 };
    const double ttbarMass[nBinsRec+1] = { 345.0, 400.0, 440.0, 485.0, 540.0, 610.0, 700.0, 850.0, 1200.0 };
    const double topPt[nBinsRec+1] = { 0.0, 60.0, 90.0, 125.0, 160.0, 200.0, 250.0, 330.0, 500.0 };
    const double deltaPhiTtbar[nBinsRec+1] = { 0.0, pi/4, pi/2, 5*pi/8, 3*pi/4, 13*pi/16, 7*pi/8, 15*pi/16, pi };
    const double deltaPhiTopZ[nBinsRec+1] = { 0.0, pi/6, pi/3, pi/2, 2*pi/3, 3*pi/4, 5*pi/6, 11*pi/12, pi };
    const double deltaRapTtbar[nBinsRec+1] = { 0.0, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 2.1, 3.0 };
    const double deltaRapTopZ[nBinsRec+1] = { 0.0, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 2.0, 3.0 };
}

}
