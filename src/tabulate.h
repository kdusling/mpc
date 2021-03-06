#ifndef tabulate_H_INCLUDED
#define tabulate_H_INCLUDED

//#define DEBUG

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include "wf.h"
#include "glasma.h"
#include "jet.h"
#include "bfkl.h"
#include <gsl/gsl_errno.h>

#define END_MSG 0
#define PIPE_MSG 1
#define RET_MSG 2

int myid;
//number of processors
int np;
//each point in phi is done on an individual slave node
//the code is setup such that the number of points in phi is the number of processors - 1 

MPI_Status Status;

int master_main(int, char**) ;
int slave_main(int, char**) ;


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
