#include "jet.h"

struct jet_params{
   double pT;
   double qT;
   double phi;
   double x1, x2;
   double dyQMRK;
} ;

void mrk_integrand(unsigned ndim, const double *x, void *params, unsigned fdim, double *fval) {
   struct jet_params myparams;
   myparams = *(struct jet_params *) params;

   double phi1 = x[0];
   double kT = x[1];
   double pT = myparams.pT;
   double qT = myparams.qT;
   double x1 = myparams.x1;
   double x2 = myparams.x2;
   double phi = myparams.phi;
   double kpppq = sqrt( pT*pT + qT*qT + kT*kT + 2.*pT*kT*cos(phi1) + 2.*qT*kT*cos(phi-phi1) + 2.*qT*pT*cos(phi) );

	fval[0] =  kT*wf(1,x1,kT)*wf(2,x2,kpppq);
}

void qmrk_integrand(unsigned ndim, const double *x, void *params, unsigned fdim, double *fval){
   struct jet_params myparams;
   myparams = *(struct jet_params *) params;   

   double kT = x[1];
   double pT = myparams.pT;
   double qT = myparams.qT;
   double x1 = myparams.x1;
   double x2 = myparams.x2;
   double phi = x[0];
   double phiq = myparams.phi;
   double k2 = sqrt( qT*qT + pT*pT + kT*kT - 2.*qT*kT*cos(phiq+phi) - 2.*pT*kT*cos(phi) + 2.*qT*pT*  cos(phiq) );
   double dyQMRK = myparams.dyQMRK;
   double E = exp(1.);

      double A;
      #include "A.h"
      double shat = 2*(-(pT*qT*cos(phiq)) + pT*qT*cosh(dyQMRK));
      double mu2 = pow(0.140,2.);
      fval[0] = kT*wf(1,x1,kT)*wf(2,x2,k2)*A/(kT*kT*k2*k2)*shat/(mu2+shat);
}

double d2N_MRK(double pT, double qT, double phi, double yp, double yq, double rts)
{
   struct jet_params params;
   params.pT = pT;
   params.qT = qT;
   params.phi = phi;
   params.x1 = (pT*exp(-yp)+qT*exp(-yq))/rts;
   params.x2 = (pT*exp(+yp)+qT*exp(+yq))/rts;

   double result, error;
   double xmin[2] = {0.0, kT_min};
   double xmax[2] = {(2.0*pi), kT_max};
   adapt_integrate(1, mrk_integrand, &params, 2, xmin, xmax, 5000000, 0, 1.e-5, &result, &error);
   
   double fac = Nc*Nc/4./pow(pi,8.)/(Nc*Nc-1.0)*alpha(pT)*alpha(qT)/(pT*pT*qT*qT);

return fac*result;
}

double d2N_QMRK(double pT, double qT, double phi, double yp, double yq, double rts)
{
   struct jet_params params;
   params.pT = pT;
   params.qT = qT;
   params.phi = phi;
   params.x1 = (pT*exp(-yp)+qT*exp(-yq))/rts;
   params.x2 = (pT*exp(+yp)+qT*exp(+yq))/rts;

params.dyQMRK = yp-yq;
if (yp == yq)
   params.dyQMRK = 1.e-10;
if (phi == 0)
   params.phi=1.e-3;

   double result, error;
   double xmin[2] = {0.0, kT_min};
   double xmax[2] = {(2.0*pi), kT_max};
   adapt_integrate(1, qmrk_integrand, &params, 2, xmin, xmax, 5000000, 0, 1.e-5, &result, &error);
   
   double fac = Nc*Nc/4./pow(pi,8.)/(Nc*Nc-1.0)*alpha(pT)*alpha(qT);

return fac*result;
}




