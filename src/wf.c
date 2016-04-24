#include "wf.h"

double wf(int nucleus, double x, double k)
{
   if ( x <= XMIN || x >= 1.0 || k < kT_min )
      return 0.0;

   if (k > kT_max)
   {
      //get index corresponding to max kT
      int i = Nkt-1;
      while (kT[0][i] > kT_max)
         i--;
      double kT1 = kT[0][i];
      double kT2 = kT[0][i-1];

      double b = log( wf(nucleus, x, kT1 )/wf(nucleus, x, kT2 ) ) / log( kT1/kT2 );
      double a = wf(nucleus, x, kT1)/pow(kT1,b);
      return a*pow(k,b);
   }
   
   double norm = Nc/4.*k*k/alpha(k) ; 
   
   double myY = gsl_max(0.0,log(x0/x));
   double dY = Y[1]-Y[0];  // width of rapidity bin
   int iY = (int) (myY/dY);    // rapidity bin
   
   if (iY+1 < NY) {
	   double tmp,tmp1,tmp2;
      if (nucleus == 1){
	      tmp1 = gsl_spline_eval(WF1[iY], k, accWF1[iY]);
	      tmp2 = gsl_spline_eval(WF1[iY+1], k, accWF1[iY+1]);
      } else if (nucleus == 2){
	      tmp1 = gsl_spline_eval(WF2[iY], k, accWF2[iY]);
	      tmp2 = gsl_spline_eval(WF2[iY+1], k, accWF2[iY+1]);
      } else {
         tmp=0; tmp1=0; tmp2=0;
         printf("Error: invalid Nucleus\n");
         return 0.;
      }
	   tmp = tmp1 + (myY-dY*iY)/dY * (tmp2-tmp1);  // lin interpolation in Y
       /* large-x correction */
	   return norm*tmp*(x > x0 ? gsl_spline_eval(LargeX, x, accLargeX) : 1.0);
	   //return norm*tmp*(x > x0 ? pow( (1. - x)/(1. - x0), 4.)  : 1.0);
   } else {
	   //printf("x too small in Phi(x)!\n");
	   return 0.;
   }

return 0.;
}

void ReadInWF(const char *wf1, const char *wf2, const char *lgx)
{

#ifdef DEBUG 
   printf("Reading in wf 1 from file: %s\n", wf1);
   printf("Reading in wf 2 from file: %s\n", wf2);
#endif

double temp;
FILE *table1 = fopen(wf1,"r");
FILE *table2 = fopen(wf2,"r");
FILE *table3 = fopen(lgx,"r");

int iY, ikT;
for (iY=0; iY < NY; iY++) {
   for (ikT=0; ikT < Nkt; ikT++) {
   	if (fscanf(table1,"%lf %lf %lf %lf", &Y[iY], &kT[iY][ikT], &temp, &uGD1[iY][ikT]) != 4) {
	      printf("Error reading table 1!\n");
	      exit(1);
	   }
   	if (fscanf(table2,"%lf %lf %lf %lf", &Y[iY], &kT[iY][ikT], &temp, &uGD2[iY][ikT]) != 4) {
	      printf("Error reading table 2!\n");
	      exit(1);
	   }
 	}

   //initialize spline at each rapidity
   accWF1[iY] = gsl_interp_accel_alloc();
   WF1[iY] = gsl_spline_alloc(gsl_interp_cspline, Nkt);
   gsl_spline_init (WF1[iY], &kT[iY][0], &uGD1[iY][0], Nkt);

   accWF2[iY] = gsl_interp_accel_alloc();
   WF2[iY] = gsl_spline_alloc(gsl_interp_cspline, Nkt);
   gsl_spline_init (WF2[iY], &kT[iY][0], &uGD2[iY][0], Nkt);
}

double xvals[100];
double LxTable[100];

int iLx;
for (iLx = 0; iLx < 100; iLx++) {
   if (fscanf(table3,"%lf %lf", &xvals[iLx], &LxTable[iLx]) != 2) {
         printf("Error reading table Large X!\n");
         exit(1);
      }
}   
   
accLargeX = gsl_interp_accel_alloc();
LargeX = gsl_spline_alloc(gsl_interp_cspline, 100);
gsl_spline_init (LargeX, &xvals[0], &LxTable[0], 100);

fclose(table1);
fclose(table2);
fclose(table3);

return;
}

void PrintWF(double x, const char *out)
{
   char fname[256];
   sprintf(fname,"%s_wf_x_%0.5f.dat",out,x);
   printf("Writing wf to file %s\n",fname);
   FILE *wfout = fopen(fname,"w");
   double k;
   for(k = kT_min; k <= 2.*kT_max ; k += 0.01)
      fprintf(wfout,"%10.3e\t%10.5e\t%10.5e\t%10.5e\n",k, alpha(k), wf(1,x,k), wf(2,x,k) );

   fclose(wfout);

}
