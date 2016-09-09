#ifndef alphas_H_INCLUDED
#define alphas_H_INCLUDED

//for fixed alpha
//#define g2 (4.6)

#include <gsl/gsl_spline.h>

//http://pdg.lbl.gov/2010/reviews/rpp2010-rev-qcd.pdf

#define Mz 91.1876
#define LamQCD 0.217
//fixednf = 0 for nf(mu)
#define fixednf 3
//extrapolates below mu0*LamQCD to value of (mu0-1.0)*LamQCD at mu=0
#define mu0 4.
#define pi 3.14159265359
#define Nc (3.)
#define Cf ((Nc*Nc-1.)/(2.*Nc))

void TabulateAlpha();
double alpha(double mu);

double nf(double mu);
double alpha_b0(double nf);
double alpha_b1(double nf);
double alpha_b2(double nf);
double alpha_s0(double mu);
double alpha_s1(double mu);
double alpha_s2(double mu);

gsl_interp_accel *accAlpha;
gsl_spline *Alpha;

#endif
