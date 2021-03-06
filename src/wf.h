#ifndef wf_H_INCLUDED
#define wf_H_INCLUDED

//#define DEBUG

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_integration.h>
#include "alphas.h"

//*
// Global variables for UGD for nucleus 1 and 2
//*

//run at fixed x
//#define fixedx (0.001)

//number of tabulated points in Y and kT for wavefunction
#define NY 151
#define Nkt 101
#define kT_min (0.02)   // phi=0 below this min kt
#define kT_max (18.)    // phi has a power lower extrapolation 
                        // beyond this max kt

#define x0 (0.01)       // separates small/large x
#define XMIN (x0*exp(-15.)) //smallest x tabulated

// Global variables for UGD for nucleus 1 and 2
int wfTAG;
double uGD1[NY][Nkt];
double uGD2[NY][Nkt];
double Y[NY];
double kT[NY][Nkt];
gsl_interp_accel *accWF1[NY];
gsl_spline *WF1[NY];
gsl_interp_accel *accWF2[NY];
gsl_spline *WF2[NY];
gsl_interp_accel *accLargeX;
gsl_spline *LargeX;
//global variable for wfTAG=2
#define QS02 (0.1)
double QS21, QS22;

void ReadInWF(int A1, int A2, int TAG);
double wf(int nucleus, double x, double kT) ;
double wfsimple(int nucleus, double x, double kT) ;
void PrintWF(double x, const char *out) ;

#endif
