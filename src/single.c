#include "single.h"

void single_integrand(unsigned ndim, const double *x, void *params, unsigned fdim, double *fval) {
   struct single_params myparams;
   myparams = *(struct single_params *) params;

   double rts = myparams.rts;
   double pT, yp, x1, x2;

   double phi = x[0];
   double kT = x[1];

   pT = myparams.pT;
   yp = myparams.yp;
      
   x1 = pT/rts*exp(-yp);
   x2 = pT/rts*exp(+yp);

   double pTpkT = 0.5*sqrt( pT*pT + kT*kT + 2.0*pT*kT*cos(phi) );
   double pTmkT = 0.5*sqrt( pT*pT + kT*kT - 2.0*pT*kT*cos(phi) );

   double norm = 8.0/Cf/pow(2.*M_PI,6.)/pow(pT,2.0)*alpha(pT)/4.;
   fval[0] = norm*kT*wf(1,x1,pTmkT)*wf(2,x2,pTpkT);

return;
}

double dNdpTdy(double pT, double yp, double rts)
{
   struct single_params params;
   params.pT = pT;
   params.yp = yp;
   params.rts = rts;

   double result, error;
   double xmin[2] = {-pi,0.};
   double xmax[2] = {+pi,50.*pT};
   
   adapt_integrate(1, single_integrand, &params, 2, xmin, xmax, 0, 0, 1.e-6, &result, &error);

return result;
}

void TabulateSingle(FILE *out, double rts)
{
   //returns dN/d^2p_T dy_p / (S_\perp) [GeV^0]
   //as a function of y and sqrt(pT)

   double yp, rtpT;
   for (yp = -3; yp <= 3; yp += 0.25)
   for (rtpT = .1; rtpT <= 5.2; rtpT += .2){
   fprintf(out,"%10.2e\t%10.2e\t%10.5e\n", yp, rtpT,\
   dNdpTdy(pow(rtpT,2.),yp,rts) );
   }
   fclose(out);
}


