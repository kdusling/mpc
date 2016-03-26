#include "alphas.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double alpha(double mu)
{
	return alpha_nf3_fit(mu);
}

double alpha_b0(double nf)
{
   return 11. - 2./3.*nf;
}

double alpha_b1(double nf)
{
   return 51. - 19./3.*nf;
}

double alpha_b2(double nf)
{
   return 2857. - 5033./9.*nf + 325./27.*nf*nf;
}

double nf(double mu)
{
if (fixednf == 0)
   return 3. + (mu >= 1.270 ? 1. : 0.) + (mu >= 4.190 ? 1. : 0.) + (mu > 172. ? 1. : 0.) ;
else 
   return fixednf;
}

double alpha_s0(double mu)
{

   return 4.*pi/( alpha_b0(nf(mu))*2.0*log(mu/LamQCD) );
}

double alpha_s1(double mu)
{
   
   return alpha_s0(mu)*(1. - alpha_b1(nf(mu))/pow(alpha_b0(nf(mu)),2.)*log(2.*log(mu/LamQCD))/log(mu/LamQCD) );
}

double alpha_s2(double mu)
{
   
   return alpha_s1(mu) + alpha_s0(mu)*( pow( alpha_b1(nf(mu))/pow(alpha_b0(nf(mu)),2.)/log(mu/LamQCD) ,2.)*( pow( log(2.*log(mu/LamQCD)) - 0.5, 2.) + alpha_b2(nf(mu))*alpha_b0(nf(mu))/8.0/pow(alpha_b1(nf(mu)),2.) - 5./4.) );
}

double alpha_nf3_fit(double mu)
{

if (mu < mu0*LamQCD){
    //double x = mu0*LamQCD;
    //double t2 = (a+b/x)/log(e*x) ;
    //double slope = -t2/x/log(e*x)-b/x/x/log(e*x);
    double t2 = 0.3703727875;
    double slope = -0.2644016513;
    return 2.*t2/(1 + exp(2.*slope/t2*(mu0*LamQCD - mu))) ;
}

double a = 0.665861;
double b = 0.0963647;
double e = 9.38513;

return (a+b/mu)/log(e*mu);

}

void print_alpha(const char *out)
{
   char fname[256];
   sprintf(fname,"%s_alphas.dat",out);
   printf("Writing alpha_s to file %s\n",fname);
   FILE *alpha_test = fopen(fname,"w");
   for (double lmu = -7; lmu < 7; lmu += .1)
   {
      fprintf(alpha_test, "%lf\t%lf\t%lf\t%lf\t%lf\n", exp(lmu),\
      alpha_s0(exp(lmu)),\
      alpha_s1(exp(lmu)),\
      alpha_s2(exp(lmu)),\
      alpha_nf3_fit(exp(lmu)) );
   }
   fclose(alpha_test);
}

