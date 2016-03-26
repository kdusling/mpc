/*! \file single.h
    \brief
*/
#ifndef single_H_INCLUDED
#define single_H_INCLUDED

#include "wf.h"
#include "cubature.h"

struct single_params{
   double pT;
   double yp;
   double rts;
   int integral; 
} ;

/*! \fn  dNdpTdy(double pT, double yp, double rts)
    \brief \f$\frac{1}{S_\perp}\frac{dN}{d^2p_T dy_p} [GeV^0] \f$
*/
double dNdpTdy(double pT, double yp, double rts);

/*! \fn  Ntot(double pTmin, double pTmax, double ymin, double ymax, double rts)
    \brief \f$\frac{1}{S_\perp}{N}_{tot} [GeV^2] \f$
*/
double Ntot(double pTmin, double pTmax, double ymin, double ymax, double rts);

void TabulateNtot(FILE *out, double rts);

void TabulateSingle(FILE *out, double rts);

#endif
