#include "main.h"

double d2N_zero(double pT,double qT,double phi,double yp,double yq,double rts)
{
   return 0.0;
}

int main(int argc, char **argv)
{

//holds name of inputfile
char inputFileName[80];

//get input file name
   if (argc != 2)
   {
      printf("Error: no input file passed through command line.\n");
      strcpy(inputFileName,"test");
      printf("Using default input file: \"%s.input\"\n", inputFileName);
   } else {
      strcpy(inputFileName,argv[1]);
      printf("Using input file: \"%s.input\"\n", inputFileName);
   }

//read in parameters
OpenFile(input, inputFileName, ".input", "r");
int ierr = ReadParameters (input);
if (ierr == 1)
   printf("Error in ReadParameters\n");

strcat(tag,"output/");
strcat(tag,inputFileName);

gsl_set_error_handler_off ();
//read in UGD from file
char fname1[256];
char fname2[256];
char fname3[256];

if (wfTAG == 0){
   sprintf(fname1,"../wf/ft_g1119_qs02_0168_af1_399_N%d.dat",A1);
   sprintf(fname2,"../wf/ft_g1119_qs02_0168_af1_399_N%d.dat",A2);
   sprintf(fname3,"../wf/LargeX.dat");
}
if (wfTAG == 1){
   sprintf(fname1,"../wf/ft_mv_qs02_02_af1_N%d.dat",A1);
   sprintf(fname2,"../wf/ft_mv_qs02_02_af1_N%d.dat",A2);
   sprintf(fname3,"../wf/LargeX.dat");
}

ReadInWF(fname1,fname2,fname3);

print_alpha(tag);
PrintWF(.1, tag);
PrintWF(.01, tag);
PrintWF(.0001, tag);

//OpenFile(outsingle, tag, "_single.dat", "w");
//printf("Writing single inclusive spectra to file: %s\n", fn_filename);
//TabulateSingle(outsingle,rts);

//   double ptT;
//   for (ptT = .1; ptT <= 20.001; ptT += .1){
//   fprintf(stdout,"\t%10.5e\n", dNdpTdy(ptT,0.,rts) );
//   }


return 0;

ReadInBFKL();  
//SpecialFuncTests(1);
//print_alpha();
setup_double();

char fname[256];
strcpy(fname,tag);
strcat(fname,"_jet.dat");
FILE *out = fopen(fname, "w"); 
printf("Writing jet spectra to file: %s\n", fname);
   
double pT = 16.81;
double qT = 16.81;
double yp = -1.;
double yq = -3.;

int nphi = 32;
double pq; /**< Angle between \f$pT\f$ and \f$q_T\f$ */

//double mrk_zyam = d2N_MRK(pT, qT, 0., yp, yq, rts);
//double qmrk_zyam = d2N_QMRK(pT, qT, 0., yp, yq, rts);
//double bfkl_zyam = d2N_BFKL(pT, qT, 0., yp, yq, rts);

fprintf(out, "# p_T=%10.2e\tq_T=%10.2e\ty_p=%10.2e\ty_q=%10.2e\n",pT,qT,yp,yq);  		
fprintf(out, "# x_P: %10.5e\n", (pT*exp(-yp)+qT*exp(-yq))/rts );
fprintf(out, "# x_T: %10.5e\n", (pT*exp(+yp)+qT*exp(+yq))/rts );
//fprintf(out, "# MRK zyam=%10.5e\tQMRK zyam=%10.5e\tBFKL zyam=%10.5e\n",\
              mrk_zyam, qmrk_zyam, bfkl_zyam);

//for (int n = 0; n < nphi;  n++) {
//     pq = (n != nphi - 1 ? n*(pi)/(double)(nphi-1) : 3.14);
//  		fprintf(out,"%10.2e\t%10.5e\t%10.5e\t%10.5e\n", pq,\
//      d2N_MRK(pT, qT, pq, yp, yq, rts),\
//      d2N_QMRK(pT, qT, pq, yp, yq, rts),\
//      d2N_BFKL(pT, qT, pq, yp, yq, rts) );
//      fflush(out);
//}  

nphi = 64;
for (int n = 0; n < nphi;  n++) {
      pq = (n != nphi - 1 ? n*(pi)/(double)(nphi-1) : 3.14);
  		printf("%10.2e\t%10.5e\t%10.5e\n", pq,\
      d2N_Glasma(pT, qT, pq, yp, yq, rts),\
      d2N_Glasma_DeltaFn(pT, qT, pq, yp, yq, rts) );
}  


return 0;
}
