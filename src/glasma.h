#ifndef glasma_H_INCLUDED
#define glasma_H_INCLUDED

#define Cf ((Nc*Nc-1.)/(2.*Nc))

struct glasma_params {
	double pT;
    double qT;
    double yp;
    double yq;
	double phi;
    double phipq;
    double rts;
    gsl_function Kernel; 
	};

void setup_double();

double d2Nglasma0(double pT, double qT, double phiq, double yp, double yq, double rts);
double d2Nglasma1(double pT, double yp, double yq, double rts);

//does phi integration over d^2k_\perp
double phikernel(double x, void *p);

//integrand functions
double doubleKernel(double x, void *p);
double cosKernel(double x, void *p);
double sinKernel(double x, void *p);

void TabulateGlasma(FILE *out, double rts);

#endif
