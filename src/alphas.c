#include "alphas.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void TabulateAlpha()
{
    int pts = 100;
    double mu[pts];
    double a[pts];

    mu[0] = 0.0;  a[0] = alpha_s2((mu0-1.0)*LamQCD);
    int i;
    for (i=1; i<pts; i++)
    {
        mu[i] = exp( log(mu0*LamQCD) + (i-1)*.1 );
        a[i] = alpha_s2(mu[i]);  
        //printf("%10.5e\t%10.5e\n",mu[i],a[i]);
    }
    accAlpha = gsl_interp_accel_alloc();
    Alpha = gsl_spline_alloc(gsl_interp_cspline, pts);
    gsl_spline_init (Alpha, &mu[0], &a[0], pts); 

return ;
}

double alpha(double mu)
{
    #ifdef g2
        return g2/(4.*pi);
    #endif 

    return gsl_spline_eval(Alpha, mu, accAlpha);
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


