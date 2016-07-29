#ifndef jetms_H_INCLUDED
#define jetms_H_INCLUDED

#include "wf.h"
#include "cubature.h"
#include "single.h"


double d2N_MRK(double pT, double qT, double phiq, double yp, double yq, double rts);
double d2N_simple_jet(double pT, double qT, double phi, double yp, double yq, double rts);

#endif
