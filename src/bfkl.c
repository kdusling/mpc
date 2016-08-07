#include "bfkl.h"

//Global variables to store BFKL table
double BFKLKernel[Nhar+1][Nz][Nx];
double Z[Nz];
double X[Nz][Nx];
gsl_interp_accel *accBFKL[Nhar+1][Nz];
gsl_spline *BFKL[Nhar+1][Nz];

static struct bfkl_params myparams;


void ReadInBFKL()
{
//reads in BFKL table from file
int iN, iZ, iX;
FILE *BFKLtable = fopen("BFKLKernel.table","r");
printf("Reading in bfkl table from file: %s\n","BFKLKernel.table");

for (iZ=0; iZ < Nz; iZ++) {
      for (iX=0; iX < Nx; iX++) {
	      if (fscanf(BFKLtable,"%lf %lf %lf", &Z[iZ], &X[iZ][iX], &BFKLKernel[0][iZ][iX]) != 3) {
	      printf("Error reading BFKLtable!\n");
	      exit(1);
	      }

         for (iN = 1; iN <= Nhar; iN++){
            if (fscanf(BFKLtable,"%lf", &BFKLKernel[iN][iZ][iX]) != 1){
	         printf("Error reading BFKLtable!\n");
	         exit(1);
            }
         }
    	}
	
      for (iN = 0; iN <= Nhar; iN++){
      	accBFKL[iN][iZ] = gsl_interp_accel_alloc();
      	BFKL[iN][iZ] = gsl_spline_alloc(gsl_interp_cspline, Nx);
      	gsl_spline_init (BFKL[iN][iZ], &X[iZ][0], &BFKLKernel[iN][iZ][0], Nx);
	}
}

return;
}

double BFKLfunc(double z, double x, int n)
{
  double sign = 1.;

  if (x < 1){
      x = 1/x;
      sign = -1.;
   }

   if (x >= XMAX || z < ZMIN || z > ZMAX || n >= Nhar)
	return 0.;


   double dZ=Z[1]-Z[0];  // width of rapidity bin
   int iZ = (int) ( (z-ZMIN)/dZ );    // rapidity bin

   double r = sqrt( x-1. );

   if (iZ+1<Nz) {
      double tmp,tmp1,tmp2;
	   tmp1 = gsl_spline_eval(BFKL[n][iZ], r, accBFKL[n][iZ]);
	   tmp2 = gsl_spline_eval(BFKL[n][iZ+1], r, accBFKL[n][iZ+1]);
	   tmp = tmp1 + (z-dZ*iZ-ZMIN)/dZ * (tmp2-tmp1);  // lin interpolation in Y

	return sign*tmp;
   } 

return 0.;
}


void bfklint(unsigned ndim, const double *x, void *params, unsigned fdim, double *fval) {
   myparams = *(struct bfkl_params *)params;
   double phi0 = x[0];
   double phi3 = x[1];
   double kT3 = x[2];
   double qratio2 = x[3];
   double kT0 = kT3*sqrt(qratio2);
   double Jac = kT3/( 2.*pow(qratio2,0.5) );
   double pT = myparams.pT;
   double qT = myparams.qT;
   double phiq = myparams.phi;
   double x1 = myparams.x1;
   double x2 = myparams.x2;
   double dy = myparams.dy;

   double pTpk0 = sqrt(pT*pT + kT0*kT0 + 2.*pT*kT0*cos(phi0));
   double qTmk3 = sqrt(qT*qT + kT3*kT3 - 2.*qT*kT3*cos(phiq-phi3));
   double alphasBar = Nc/pi*sqrt( alpha(kT0) * alpha(kT3) );
  
   //BFKLfunc( z = alphasbar * y, x = (qa/qb)^2, N)
   double BFKL = 1./(kT0*kT3)*BFKLfunc(alphasBar*fabs(dy),qratio2,0);
   
   int k; 
   for (k=1; k <= Nhar; k++)
	   BFKL += 2./(kT0*kT3)*BFKLfunc(alphasBar*fabs(dy),qratio2,k)*cos( (double)k*(phi3-phi0) );

   fval[0] = Jac*kT0*kT3*wf(1,x1,pTpk0)*wf(2,x2,qTmk3)*BFKL/(1.e-10+fabs(log(qratio2)));
   

return;
}

double d2N_BFKL(double pT, double qT, double phi, double yp, double yq, double rts)
{
   double result1, result2, error1, error2;
   struct bfkl_params params;
   params.pT = pT;
   params.qT = qT;
   params.phi = phi;
   params.dy = yp-yq;
   params.x1 = (pT*exp(-yp)+qT*exp(-yq))/rts;
   params.x2 = (pT*exp(+yp)+qT*exp(+yq))/rts;
   
   double xmin1[4] = {0., 0., kT_min, 1.0/XMAX};
   double xmax1[4] = {2.*M_PI, 2.*M_PI, 3.0*kT_max, 1.00};
   adapt_integrate(1, bfklint, &params, 4, xmin1, xmax1, 10000000, 0, 1e-4, &result1, &error1);

   double xmin2[4] = {0., 0., kT_min, 1.00};
   double xmax2[4] = {2.*M_PI, 2.*M_PI, 3.0*kT_max, XMAX};
   adapt_integrate(1, bfklint, &params, 4, xmin2, xmax2, 10000000, 0, 1e-4, &result2, &error2);
   
   //two times MRK factor 
   double fac = Nc*Nc/2./pow(pi,8.)/(Nc*Nc-1.0)*alpha(pT)*alpha(qT)/(pT*pT*qT*qT);

return fac*(result2 - result1);
}




