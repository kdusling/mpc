#include <stdio.h>
#include <math.h>
#include "wf.h"
#include "alphas.h"
#include "glasma.h"
#include <gsl/gsl_integration.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>

//workspace and global variables for integration
static double abserr = 1.e-20;
static double relerr = 1.e-4;
static size_t wktsize = 20;
static gsl_integration_workspace * wkt;
static struct spec_params myparams;

double ktphidouble(double x, void *params);
double phidouble(double x, void *params);

double Glasma_Kernel(double x, void *params) ;
double GlasmaI1_Kernel(double x, void *params) ;
double GlasmaI2_Kernel(double x, void *params) ;
double GlasmaI3_Kernel(double x, void *params) ;

static double result,error;
static size_t neval;
static gsl_function F;
static struct spec_params params;

void setup_double()
{
wkt = gsl_integration_workspace_alloc (wktsize);
gsl_set_error_handler_off ();
}

double d2N_Glasma(double pT, double qT, double phiq, double yp, double yq, double rts)
{
params.type = 2;
params.pT = pT; params.qT = qT; params.phiq = phiq;
	
params.swap = 0;
if (yp >= yq){
   yp *= -1.;
   yq *= -1.;
   params.swap = 1;
}

params.x1 = pT/rts*exp(-yp);
params.x2 = pT/rts*exp(+yp);
params.z1 = qT/rts*exp(-yq);
params.z2 = qT/rts*exp(+yq);
//do phi integration
F.params = &params;
F.function = &phidouble;
gsl_integration_qng(&F, 0., 2.*M_PI,abserr,relerr,&result,&error,&neval);

return result;
}

double In(double pT, double yp, double yq, double rts, int n)
{
params.pT = pT; params.qT = pT; params.phiq = 0;
double fac = sqrt( (Nc*Nc)/gsl_pow_3(Nc*Nc-1.)/(4.*pow(M_PI,10.))/gsl_pow_2(pT*pT)*alpha(pT)*alpha(pT) );

params.swap = 0;

if (yp >= yq){
	yp *= -1.;
	yq *= -1.;
   params.swap = 1;
}

params.x1 = pT/rts*exp(-yp);
params.x2 = pT/rts*exp(+yp);
params.z1 = pT/rts*exp(-yq);
params.z2 = pT/rts*exp(+yq);
       
F.function = &phidouble;
F.params = &params;

double I;

if (n == 1)
   params.type = 11;

if (n == 2)
   params.type = 12;

if (n == 3)
   params.type = 13;

gsl_integration_qng(&F, 0., 2.*M_PI,abserr,relerr,&I,&error,&neval);
     
return I/fac;
} 

double d2N_Glasma_DeltaFn(double pT, double qT, double phiq, double yp, double yq, double rts)
{
pT = (pT+qT)/2.;
qT = pT;
params.pT = pT; params.qT = qT; params.phiq = phiq;

params.swap = 0;

if (yp >= yq){
	yp *= -1.;
	yq *= -1.;
   params.swap = 1;
}

params.x1 = pT/rts*exp(-yp);
params.x2 = pT/rts*exp(+yp);
params.z1 = qT/rts*exp(-yq);
params.z2 = qT/rts*exp(+yq);
       
F.function = &phidouble;
F.params = &params;

double I1, I2, I3;

params.type = 11;
gsl_integration_qng(&F, 0., 2.*M_PI,abserr,relerr,&I1,&error,&neval);
params.type = 12;
gsl_integration_qng(&F, 0., 2.*M_PI,abserr,relerr,&I2,&error,&neval);
params.type = 13;
gsl_integration_qng(&F, 0., 2.*M_PI,abserr,relerr,&I3,&error,&neval);
      
return (I1*I1 + I2*I2 + 2.*I3*I3) ;
}

//does phi integration over d^2k_\perp
double phidouble(double x, void *params)
{
   myparams = *(struct spec_params *)params;
   myparams.phi = x;

   double result,error;
   gsl_function F;
   F.function = &ktphidouble;
   F.params = &myparams;

   gsl_integration_qagiu (&F , kT_min, abserr, relerr, wktsize, wkt, &result, &error);
return result;
}



double ktphidouble(double x, void *params)
{
myparams = *(struct spec_params *)params;
int type = myparams.type;


switch (type){
	case 11:
      return GlasmaI1_Kernel(x, params);
	break;
	case 12:
      return GlasmaI2_Kernel(x, params);
	break;
	case 13:
      return GlasmaI3_Kernel(x, params);
	break;
	

   case 2:
      return Glasma_Kernel(x, params);
	break;
}

return result;
}

double Glasma_Kernel(double x, void *params)
{
   myparams = *(struct spec_params *)params;

   double kT = x;
   double pT = myparams.pT;
   double qT = myparams.qT;
   double x1 = myparams.x1;
   double x2 = myparams.x2;
   double z1 = myparams.z1;
   double z2 = myparams.z2;
   double phi = myparams.phi;
   double phiq = myparams.phiq;
   double phi2 = phiq-phi;
   double pTpkT = sqrt(pT*pT + kT*kT + 2.*pT*kT*cos(phi ));
   double pTmkT = sqrt(pT*pT + kT*kT - 2.*pT*kT*cos(phi ));
   double qTpkT = sqrt(qT*qT + kT*kT + 2.*qT*kT*cos(phi2));
   double qTmkT = sqrt(qT*qT + kT*kT - 2.*qT*kT*cos(phi2));
   double pTmqTmkT = sqrt( pT*pT + qT*qT + kT*kT - 2.*pT*qT*cos(phiq) + 2.*qT*kT*cos(phi2) - 2.*pT*kT*cos(phi) );
   double pTpqTmkT = sqrt( pT*pT + qT*qT + kT*kT + 2.*pT*qT*cos(phiq) - 2.*qT*kT*cos(phi2) - 2.*pT*kT*cos(phi) );

	double fac = (Nc*Nc)/gsl_pow_3(Nc*Nc-1.)/(4.*pow(M_PI,10.))/gsl_pow_2(pT*qT)*alpha(pT)*alpha(qT);
	
   double res1=0; double res2=0; double res5=0; double res7=0;
   double num1, num2, den;
	   
   //res1+res2 is the F^{1,2,3,6}
   if (myparams.swap == 0){
	   //the factor of 1/2 comes from symmetrizing the pTmkT wf which is not necessary but does help out the numerics
      res1 = 0.5*fac*kT*pow(wf(1,x1,kT),2.)*( wf(2,x2,pTmkT) + wf(2,x2,pTpkT) )*( wf(2,z2,qTpkT) + wf(2,z2,qTmkT) ) ; 
	   res2 = 0.5*fac*kT*pow(wf(2,z2,kT),2.)*( wf(1,x1,pTmkT) + wf(1,x1,pTpkT) )*( wf(1,z1,qTpkT) + wf(1,z1,qTmkT) ) ;
   }
   if (myparams.swap == 1){
	   res1 = 0.5*fac*kT*pow(wf(2,x1,kT),2.)*( wf(1,x2,pTmkT) + wf(1,x2,pTpkT) )*( wf(1,z2,qTpkT) + wf(1,z2,qTmkT) ) ; 
	   res2 = 0.5*fac*kT*pow(wf(1,z2,kT),2.)*( wf(2,x1,pTmkT) + wf(2,x1,pTpkT) )*( wf(2,z1,qTpkT) + wf(2,z1,qTmkT) ) ;
   }
   //this is F^{5}
   den =  pow(kT*pTmqTmkT*pTmkT*qTpkT,2.0);
   num1 = ( kT*pT*cos(phi ) - kT*kT )*(  pT*pT - pT*qT*cos(phiq) - pT*kT*cos(phi ) - pow(pTmqTmkT,2.) ) + ( kT*pT*sin(phi ) )*( -pT*qT*sin(phiq) - pT*kT*sin(phi ) );
   num2 = ( kT*qT*cos(phi2) + kT*kT )*( -qT*qT + pT*qT*cos(phiq) - qT*kT*cos(phi2) + pow(pTmqTmkT,2.) ) + ( kT*qT*sin(phi2) )*(  pT*qT*sin(phiq) - qT*kT*sin(phi2) );
   
   //this factor of 0.5 is from the different f^4\delta^4 structure
   if (myparams.swap == 0){
      res5 = 0.5*fac*kT*wf(1,x1,kT)*wf(1,x1,pTmqTmkT)*wf(2,z2,pTmkT)*wf(2,z2,qTpkT)*( num1*num2/den );
   }   
   if (myparams.swap == 1){
      res5 = 0.5*fac*kT*wf(2,x1,kT)*wf(2,x1,pTmqTmkT)*wf(1,z2,pTmkT)*wf(1,z2,qTpkT)*( num1*num2/den );
   }
   
   //this is F^{7}
   den =  pow(kT*pTpqTmkT*pTmkT*qTmkT,2.0);
   num1 = ( kT*pT*cos(phi ) - kT*kT )*(  pT*pT + pT*qT*cos(phiq) - pT*kT*cos(phi ) - pow(pTpqTmkT,2.) ) + (  kT*pT*sin(phi ) )*(  pT*qT*sin(phiq) - pT*kT*sin(phi ) );
   num2 = ( kT*qT*cos(phi2) - kT*kT )*(  qT*qT + pT*qT*cos(phiq) - qT*kT*cos(phi2) - pow(pTpqTmkT,2.) ) + ( -kT*qT*sin(phi2) )*( -pT*qT*sin(phiq) + qT*kT*sin(phi2) );
   if (myparams.swap == 0){
      res7 = 0.5*fac*kT*wf(1,x1,kT)*wf(1,x1,pTpqTmkT)*wf(2,z2,pTmkT)*wf(2,z2,qTmkT)*( num1*num2/den );
   }   
   if (myparams.swap == 1){
      res7 = 0.5*fac*kT*wf(2,x1,kT)*wf(2,x1,pTpqTmkT)*wf(1,z2,pTmkT)*wf(1,z2,qTmkT)*( num1*num2/den );
   }   

return res1 + res2 + res5 + res7;
}	

double GlasmaI1_Kernel(double x, void *params)
{
   myparams = *(struct spec_params *)params;

   double kT = x;
   double pT = myparams.pT;
   double qT = myparams.qT;
   double x1 = myparams.x1;
   double z2 = myparams.z2;
   double phi = myparams.phi;
   double pTmkT = sqrt(pT*pT + kT*kT - 2.*pT*kT*cos(phi) );

	double fac = sqrt( (Nc*Nc)/gsl_pow_3(Nc*Nc-1.)/(4.*pow(M_PI,10.))/gsl_pow_2(pT*qT)*alpha(pT)*alpha(qT) );
   double num = pow( kT*pT*cos(phi)-kT*kT, 2.0); 	 
   double den  = pow(kT*pTmkT,2.);  
   
         if (myparams.swap == 0){
            return fac*kT*wf(1,x1,kT)*wf(2,z2,pTmkT)*num/den;
         }
         if (myparams.swap == 1){
            return fac*kT*wf(2,x1,kT)*wf(1,z2,pTmkT)*num/den;
         }

return 0;
}
		
double GlasmaI2_Kernel(double x, void *params)
{
   myparams = *(struct spec_params *)params;

   double kT = x;
   double pT = myparams.pT;
   double qT = myparams.qT;
   double x1 = myparams.x1;
   double z2 = myparams.z2;
   double phi = myparams.phi;
   double pTmkT = sqrt(pT*pT + kT*kT - 2.*pT*kT*cos(phi) );
   
   double fac = sqrt( (Nc*Nc)/gsl_pow_3(Nc*Nc-1.)/(4.*pow(M_PI,10.))/gsl_pow_2(pT*qT)*alpha(pT)*alpha(qT) );
   double num = pow( kT*pT*sin(phi), 2.0); 	 
   double den = pow(kT*pTmkT,2.);  
   
         if (myparams.swap == 0){
            return fac*kT*wf(1,x1,kT)*wf(2,z2,pTmkT)*num/den;
         } 
         if (myparams.swap == 1){
            return fac*kT*wf(2,x1,kT)*wf(1,z2,pTmkT)*num/den;
         }

return 0;
}
		
double GlasmaI3_Kernel(double x, void *params)
{
   myparams = *(struct spec_params *)params;

   double kT = x;
   double pT = myparams.pT;
   double qT = myparams.qT;
   double x1 = myparams.x1;
   double z2 = myparams.z2;
   double phi = myparams.phi;
   double pTmkT = sqrt(pT*pT + kT*kT - 2.*pT*kT*cos(phi) );

      double fac = sqrt( (Nc*Nc)/gsl_pow_3(Nc*Nc-1.)/(4.*pow(M_PI,10.))/gsl_pow_2(pT*qT)*alpha(pT)*alpha(qT) );
      double num = ( kT*pT*cos(phi)-kT*kT )*( kT*pT*sin(phi) ); 	 
      double den  = pow(kT*pTmkT,2.);  
   
         if (myparams.swap == 0){
            return fac*kT*wf(1,x1,kT)*wf(2,z2,pTmkT)*num/den;
         }
         if (myparams.swap == 1){
            return fac*kT*wf(2,x1,kT)*wf(1,z2,pTmkT)*num/den;
         }

return 0;
}
