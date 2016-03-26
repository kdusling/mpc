#ifndef alphas_H_INCLUDED
#define alphas_H_INCLUDED

//http://pdg.lbl.gov/2010/reviews/rpp2010-rev-qcd.pdf

#define Mz 91.1876
#define LamQCD 0.217
#define fixednf 3
#define mu0 4.
#define pi 3.14159265359
#define Nc (3.)
#define Cf ((Nc*Nc-1.)/(2.*Nc))

double alpha(double mu);

double alpha_b0(double nf);
double alpha_b1(double nf);
double alpha_b2(double nf);
double nf(double mu);
double alpha_s0(double mu);
double alpha_s1(double mu);
double alpha_s2(double mu);
double alpha_nf3_fit(double mu);

void print_alpha(const char *out);

#endif
