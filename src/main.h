#ifndef main_H_INCLUDED
#define main_H_INCLUDED

#include "wf.h"
#include "cubature.h"
#include "jet.h"
#include "bfkl.h"
#include "glasma.h"
#include "time.h"
#include "single.h"
//#include <gsl/gsl_errno.h>
//#include "SpecialFunc.h"

char tag[256];

//*
//*
// This section of code initializes and reads in the
// following global variables from file
//*
//*
#define X_FIELDS \
   X(int, wfTAG, "%d") \
   X(int, A1, "%d") \
   X(int, A2, "%d") \
   X(double, rts, "%lf") \

#include "ReadParameters.cxx"
//


#endif
