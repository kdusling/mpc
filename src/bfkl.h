#ifndef bfkl_H_INCLUDED
#define bfkl_H_INCLUDED

#include <stdio.h>
#include <math.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include "cubature.h"
#include "string.h"
#include "alphas.h"
#include "wf.h"

//number of tabulated points in z and x for BFKL Kernel
#define Nz 94
#define Nx 1201
#define Nhar 24 
#define ZMIN (0.15)
#define ZMAX (4.8)
#define XMAX (100.00)

void ReadInBFKL();
double BFKLfunc(double z, double x, int n);
double d2N_BFKL(double pT, double qT, double phi, double yp, double yq, double rts);

struct bfkl_params {
	double pT;
	double qT;
	double phi;
	double x1, x2;
    double dy;
	} ;

#endif
