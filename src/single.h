/*! \file single.h
    \brief
*/
#ifndef single_H_INCLUDED
#define single_H_INCLUDED

#include "wf.h"
#include "cubature.h"

#define KMAX (100.)

struct single_params{
   double pT;
   double yp;
   double rts;
   double mass;
   int integral; 
} ;

/*! \fn  dNdpTdy(double pT, double yp, double rts)
    \brief \f$\frac{1}{S_\perp}\frac{dN}{d^2p_T dy_p} [GeV^0] \f$
*/
double dNd2pTdy(double pT, double yp, double rts);
double dNdpT(double pT, double ymin, double ymax, double rts);
double dNdy(double pTmin, double yp, double rts, double mass);

/*! \fn  ng(double pTmin, double pTmax, double ymin, double ymax, double rts)
    \brief \f$\frac{1}{S_\perp}{N}_{tot} [GeV^2] \f$
*/
double ng(double pTmin, double pTmax, double ymin, double ymax, double rts);

void TabulateNtot(FILE *out, double rts);

void TabulateSingle(char *tag, double rts);

#endif
