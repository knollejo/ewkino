#ifndef ttZObservables_H
#define ttZObservables_H

#include "TMath.h"
#include "TVector2.h"

inline double deltaPhi(double phi1, double phi2) {
    const double dphi = phi1-phi2;
    if(isnan(dphi)) return -999.0;
    else return TMath::Abs(TVector2::Phi_mpi_pi(dphi));
}

inline double deltaRap(double rap1, double rap2) {
    return TMath::Abs(rap1-rap2);
}

#endif
