#include "single.h"

void single_integrand(unsigned ndim, const double *x, void *params, unsigned fdim, double *fval) {
   struct single_params myparams;
   myparams = *(struct single_params *) params;

   int integral = myparams.integral;
   double rts = myparams.rts;
   double m = myparams.mass;
   double pT, mT, yp, x1, x2, jac;

   double phi = x[0];
   double kT = x[1];

   if (integral == 1){
      yp = myparams.yp;
      pT = myparams.pT;
      jac = 1.0;
   }
   if (integral == 2){
      pT = x[2];
      yp = myparams.yp;
      mT = sqrt(m*m + pT*pT);
      jac = 2*M_PI*pT*sqrt(1.0 - m*m/pow(mT*cosh(yp),2.0) );
   }
   if (integral == 3){
      yp = x[2];
      pT = myparams.pT;
      jac = 2*M_PI*pT;
   }
   if (integral == 4){
      pT = x[2];
      yp = x[3];
      jac = 2*M_PI*pT;
   }
      
   x1 = pT/rts*exp(+yp);
   x2 = pT/rts*exp(-yp);

   double pTpkT = 0.5*sqrt( pT*pT + kT*kT + 2.0*pT*kT*cos(phi) );
   double pTmkT = 0.5*sqrt( pT*pT + kT*kT - 2.0*pT*kT*cos(phi) );

   double norm = 8.0/Cf/pow(2.*M_PI,6.)/pow(pT,2.0)*alpha(pT)/4.;
   fval[0] = jac*norm*kT*wf(1,x1,pTmkT)*wf(2,x2,pTpkT);

return;
}

double dNd2pTdy(double pT, double yp, double rts)
{
   struct single_params params;
   params.integral = 1;
   params.pT = pT;
   params.yp = yp;
   params.rts = rts;
   
   if ( pT/rts*exp(fabs(yp)) >= 1.0 || pT/rts*exp(-fabs(yp)) <= XMIN)
           return 0.;
   
   double result, error;
   double xmin[2] = {-pi,   0.};
   double xmax[2] = {+pi, KMAX};
   
   adapt_integrate(1, single_integrand, &params, 2, xmin, xmax, 100000, 0, 1.e-4, &result, &error);

return result;
}

double dNdy(double pTmin, double yp, double rts, double mass)
{
   struct single_params params;
   params.integral = 2;
   params.yp = yp;
   params.rts = rts;
   params.mass = mass;
    
   double pTmax = rts*gsl_min(exp(yp),exp(-yp));
   if (pTmax < pTmin)
      return 0.;

   double result, error;
   double xmin[3] = {-pi,   0., pTmin};
   double xmax[3] = {+pi, KMAX, rts*gsl_min(exp(yp),exp(-yp)) };
   adapt_integrate(1, single_integrand, &params, 3, xmin, xmax, 0, 0, 1.e-4, &result, &error);

return result;
}

double dNdpT(double pT, double ymin, double ymax, double rts)
{
   struct single_params params;
   params.integral = 3;
   params.pT = pT;
   params.rts = rts;
   double result, error;
   double xmin[3] = {-pi,   0., ymin};
   double xmax[3] = {+pi, KMAX, ymax};

   adapt_integrate(1, single_integrand, &params, 3, xmin, xmax, 0, 0, 1.e-4, &result, &error);

return result;
}

double ng(double pTmin, double pTmax, double ymin, double ymax, double rts)
{
   struct single_params params;
   params.integral = 4;
   params.rts = rts;
   double result, error;
   double xmin[4] = {-pi,   0., pTmin, ymin};
   double xmax[4] = {+pi, KMAX, pTmax, ymax};
   
   adapt_integrate(1, single_integrand, &params, 4, xmin, xmax, 0, 0, 1.e-4, &result, &error);

return result;
}

void TabulateSingle(char *tag, double rts)
{
   //returns dN/d^2p_T dy_p / (S_\perp) [GeV^0]
   //as a function of y and sqrt(pT)

   char fn_filename[256]; 

   strcpy(fn_filename, tag); 
   strcat(fn_filename, "_single.txt"); 
   FILE * s1 = fopen( fn_filename, "w"); 
   if (s1 == NULL) {printf("Error: file \"%s\" could not be opened.\n\n",fn_filename); exit(1);}

   double yp, rtpT;
   for (yp = -10.0; yp <= 10.0; yp += 1.0)
   for (rtpT = .1; rtpT <= 5.11; rtpT += .1){
   fprintf(s1,"%10.2f\t%10.2f\t%10.5e\n", yp, rtpT,\
   dNd2pTdy(pow(rtpT,2.),yp,rts) );
   fflush(s1);
   }
   fclose(s1);
   
   strcpy(fn_filename, tag); 
   strcat(fn_filename, "_single_midrapidity.txt"); 
   FILE * s2 = fopen( fn_filename, "w"); 
   if (s2 == NULL) {printf("Error: file \"%s\" could not be opened.\n\n",fn_filename); exit(1);}

   double rT;
   for (rT=.28125; rT<18. ; rT+=0.0625)
   {
        fprintf(s2,"%10.5f\t%10.5e\n",rT,dNd2pTdy(rT,0.,rts));
   }
   fclose(s2);
   
   return;
}


