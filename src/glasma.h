#ifndef glasma_H_INCLUDED
#define glasma_H_INCLUDED

#define Cf ((Nc*Nc-1.)/(2.*Nc))

struct spec_params {
	int type;
   int swap;
	double pT;
	double qT;
	double phi;
	double phiq;
	double x1, x2;
	double z1, z2;
	};

void setup_double();

double d2N_Glasma(double pT, double qT, double phiq, double yp, double yq, double rts);
double d2N_Glasma_DeltaFn(double pT, double qT, double phiq, double yp, double yq, double rts);
double In(double pT, double yp, double yq, double rts, int n);

#endif
