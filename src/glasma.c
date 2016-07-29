#include <stdio.h>
#include <math.h>
#include "wf.h"
#include "alphas.h"
#include "glasma.h"
#include <gsl/gsl_integration.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>

//workspace and global variables for integration
static double abserr = 1.e-200;
static double relerr = 1.e-20;

double phiKernel(double x, void *params);

static double result,error;
static size_t neval;
static gsl_function phiIntKernel;

double d2Nglasma0(double pT, double qT, double phipq, double yp, double yq, double rts)
{
    static struct glasma_params params;
    params.pT = pT;
    params.qT = qT;
    params.yp = yp;
    params.yq = yq;
    params.rts = rts;
    params.phipq = phipq;
    phiIntKernel.function = &phiKernel; 
    (params.Kernel).function = &doubleKernel;

    phiIntKernel.params = &params;
    gsl_integration_qng(&phiIntKernel, -M_PI, M_PI, abserr, relerr,\
            &result, &error, &neval);

    double norm = (Nc*Nc)/gsl_pow_3(Nc*Nc-1.)/(4.*pow(M_PI,10.))/gsl_pow_2(pT*
        qT)*alpha(pT)*alpha(qT);

    return norm*result;
}

double d2Nglasma1(double pT, double yp, double yq, double rts)
{
    double sinres, cosres;
    static struct glasma_params params;
    params.pT = pT;
    params.yp = yp;
    params.yq = yq;
    params.rts = rts;

    phiIntKernel.function = &phiKernel; 
    (params.Kernel).function = &cosKernel;
    phiIntKernel.params = &params;
    gsl_integration_qng(&phiIntKernel, -M_PI, M_PI, abserr, relerr, &cosres, &error, &neval);
    
    (params.Kernel).function = &sinKernel;
    phiIntKernel.params = &params;
    gsl_integration_qng(&phiIntKernel, -M_PI, M_PI, abserr, relerr, &sinres, &error, &neval);

    double norm =  (Nc*Nc)/gsl_pow_3(Nc*Nc-1.)/(4.*pow(M_PI,10.))/gsl_pow_2(pT*pT)*alpha(pT)*alpha(pT) ;

    return norm*(cosres*cosres + sinres*sinres);
}


double doubleKernel(double x, void *p)
{
   struct glasma_params params = *(struct glasma_params *)p;
  
   double rts = params.rts;
   double pT, qT, yp, yq;
   double x1p, x2p, x1q, x2q, x1max, x2max;

   double phi = params.phi;
   double phipq = params.phipq;
   double kT =  x;
   
   pT = params.pT;
   qT = params.qT;
   yp = params.yp;
   yq = params.yq;

   x1p = pT/rts*exp(+yp);
   x2p = pT/rts*exp(-yp);
   x1q = qT/rts*exp(+yq);
   x2q = qT/rts*exp(-yq);

   x1max = gsl_max(x1p,x1q);
   x2max = gsl_max(x2p,x2q);
   
   double qTmkT = sqrt( qT*qT + kT*kT - 2.0*qT*kT*cos(phipq-phi) );
   double qTpkT = sqrt( qT*qT + kT*kT + 2.0*qT*kT*cos(phipq-phi) );
   double pTmkT = sqrt( pT*pT + kT*kT - 2.0*pT*kT*cos(phi) );
   double pTpkT = sqrt( pT*pT + kT*kT + 2.0*pT*kT*cos(phi) );
   
   double pTmkTmqT = sqrt( pT*pT + kT*kT + qT*qT \
                          -2.0*pT*qT*cos(phipq) \
                          -2.0*pT*kT*cos(phi) \
                          +2.0*qT*kT*cos(phipq-phi)\
                        );

   double pTmkTpqT = sqrt( pT*pT + kT*kT + qT*qT \
                          +2.0*pT*qT*cos(phipq) \
                          -2.0*pT*kT*cos(phi) \
                          -2.0*qT*kT*cos(phipq-phi)\
                        );

   //diagram A & E
   //this factor of 0.5 is from
   //symmetrizing wrt pTmkT
   double AE = 0.5*kT*pow(wf(1,x1max,kT),2.0)\
          *( wf(2,x2p,pTmkT) + wf(2,x2p,pTpkT) )\
          *( wf(2,x2q,qTmkT) + wf(2,x2q,qTpkT) );

   //diagram B & F
   double BF = 0.5*kT*pow(wf(2,x2max,kT),2.0)\
          *( wf(1,x1p,pTmkT) + wf(1,x1p,pTpkT) )\
          *( wf(1,x1q,qTmkT) + wf(1,x1q,qTpkT) );
  
   double T1D, T2D, T1H, T2H, denD, denH;
   T1D = ( kT*pT*cos(phi)-kT*kT )\
        *( pT*pT - pT*qT*cos(phipq) - pT*kT*cos(phi) - pow(pTmkTmqT,2.0) ) \
        - (pT*kT*sin(phi))*(pT*qT*sin(phipq) + pT*kT*sin(phi) );        
   
   T1H = ( kT*pT*cos(phi)-kT*kT )\
        *( pT*pT + pT*qT*cos(phipq) - pT*kT*cos(phi) - pow(pTmkTpqT,2.0) ) \
        - (pT*kT*sin(phi))*(-pT*qT*sin(phipq) + pT*kT*sin(phi) );        
   
   T2D = ( kT*qT*cos(phipq-phi)+kT*kT )\
        *( -qT*qT + pT*qT*cos(phipq) - qT*kT*cos(phipq-phi) + pow(pTmkTmqT,2.0) ) \
        + (qT*kT*sin(phipq-phi))*(pT*qT*sin(phipq) - qT*kT*sin(phipq-phi) );        
   
   T2H = ( -kT*qT*cos(phipq-phi)+kT*kT )\
        *( -qT*qT - pT*qT*cos(phipq) + qT*kT*cos(phipq-phi) + pow(pTmkTpqT,2.0) ) \
        + (qT*kT*sin(phipq-phi))*(pT*qT*sin(phipq) - qT*kT*sin(phipq-phi) );        
   
   denD = pow(kT*pTmkTmqT*pTmkT*qTpkT,2.0);
   
   denH = pow(kT*pTmkTpqT*pTmkT*qTmkT,2.0);
   
   //this factor of 0.5 is from the different f^4\delta^4 structure
   double D = 0.5*kT*T1D*T2D/denD\
               *wf(1,x1max,kT)*wf(1,x1max,pTmkTmqT)*wf(2,x2max,pTmkT)*wf(2,x2max,qTpkT);
   
   double H = 0.5*kT*T1H*T2H/denH\
               *wf(1,x1max,kT)*wf(1,x1max,pTmkTpqT)*wf(2,x2max,pTmkT)*wf(2,x2max,qTmkT);
   
   return (AE + BF + D + H);
}


double cosKernel(double x, void *p)
{
   struct glasma_params params = *(struct glasma_params *)p;
  
   double rts = params.rts;
   double pT, yp, yq;
   double x1p, x2p, x1q, x2q, x1max, x2max;

   double phi = params.phi;
   double kT =  x;
   
   pT = params.pT;
   yp = params.yp;
   yq = params.yq;

   x1p = pT/rts*exp(+yp);
   x2p = pT/rts*exp(-yp);
   x1q = pT/rts*exp(+yq);
   x2q = pT/rts*exp(-yq);

   x1max = gsl_max(x1p,x1q);
   x2max = gsl_max(x2p,x2q);
   
   double pTmkT = sqrt( pT*pT + kT*kT - 2.0*pT*kT*cos(phi) );

   double num = pow( kT*pT*cos(phi)-kT*kT, 2.0); 	 
   double den  = pow(kT*pTmkT,2.);  

   return kT*wf(1,x1p,kT)*wf(2,x2p,pTmkT)*num/den;
}

double sinKernel(double x, void *p)
{
   struct glasma_params params = *(struct glasma_params *)p;
  
   double rts = params.rts;
   double pT, yp, yq;
   double x1p, x2p, x1q, x2q, x1max, x2max;

   double phi = params.phi;
   double kT =  x;
   
   pT = params.pT;
   yp = params.yp;
   yq = params.yq;

   x1p = pT/rts*exp(+yp);
   x2p = pT/rts*exp(-yp);
   x1q = pT/rts*exp(+yq);
   x2q = pT/rts*exp(-yq);

   x1max = gsl_max(x1p,x1q);
   x2max = gsl_max(x2p,x2q);
   
   double pTmkT = sqrt( pT*pT + kT*kT - 2.0*pT*kT*cos(phi) );

   double num = pow( kT*pT*sin(phi), 2.0); 	 
   double den = pow( kT*pTmkT ,2.0);  

   return kT*wf(1,x1max,kT)*wf(2,x2max,pTmkT)*num/den;
}

//does phi integration over d^2k_\perp
double phiKernel(double x, void *p)
{
   struct glasma_params params = *(struct glasma_params *)p;
   double result, error;
   params.phi = x;
   
   gsl_function F = params.Kernel;
   F.params = &params;

   gsl_integration_qng(&F, kT_min, 2.0*kT_max, abserr, relerr,\
           &result, &error, &neval);

return result;
}

void TabulateGlasma(FILE *out, double rts)
{
   //returns dN/d^2p_T d^2q_T dy_p dy_q / (S_\perp) [GeV^-2]
   //this is the delta function contribution
   //of the glasma graphs

   double yp, yq, rtpT;
   //double yrange = 6.0;
   double dy = 0.25;
   double yrange = dy;
   
   for (yp = -yrange; yp <= yrange; yp += dy)
   for (yq = -yrange; yq <= yrange; yq += dy)
   for (rtpT = .1; rtpT <= 5.11; rtpT += .2){
   fprintf(out,"%10.2f\t%10.2f\t%10.2f\t%10.5e\n", yp, yq, rtpT,\
   d2Nglasma1(pow(rtpT,2.),yp,yq,rts) );
   }
   fclose(out);

return;
}
