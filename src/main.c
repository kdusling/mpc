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
setup_double();
//ReadInBFKL();  

//print_alpha(tag);
//PrintWF(.01, tag);

OpenFile(outsingle, tag, "_single.dat", "w");
OpenFile(outdouble, tag, "_double.dat", "w");


double y, x1, x2;
double kT = 2.0;

for (y = -10; y <= 10; y += 0.05)
{ 
    x1 = kT/rts*exp(+y);
    x2 = kT/rts*exp(-y);
    fprintf(outsingle,"%10.2e\t%10.2e\t%10.2e\t%10.5e\n", y, x1, x2,  dNd2pTdy(kT,y,rts) );
}
clock_t begin, end;
double time_spent;

begin = clock();

//warning -- only using AEBF
double yp,yq,dy;
double x1p, x2p, x1q, x2q;
double pT=0.8;
double qT=0.8;
for (yp = -3; yp <= 3; yp += 3)
{
    for (dy = -8; dy <= 8.01; dy += .2)
    { 
        yq = dy + yp;
        x1p = pT/rts*exp(+yp);
        x2p = pT/rts*exp(-yp);
        x1q = qT/rts*exp(+yq);
        x2q = qT/rts*exp(-yq);
        fprintf(outdouble,"%10.2e\t%10.2e\t%10.2e\t%10.2e\t%10.2e\t%10.2e\t%10.5e\t%10.5e\n", yp, dy, x1p, x2p, x1q, x2q, \
        d2N(pT,qT,0.,yp,yq,rts),dNd2pTdy(pT,yp,rts)*dNd2pTdy(qT,yq,rts) );
    }
fprintf(outdouble,"\n");
fflush(outdouble);
}

end = clock();
time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
printf("time elapsed (seconds): %10.2e\n", time_spent);

return 0;
//char fname[256];
//strcpy(fname,tag);
//strcat(fname,"_jet.dat");
//FILE *out = fopen(fname, "w"); 
//printf("Writing jet spectra to file: %s\n", fname);
   
//double pT = 16.81;
//double qT = 16.81;
//double yp = -1.;
//double yq = -3.;

//int nphi = 32;
//double pq; /**< Angle between \f$pT\f$ and \f$q_T\f$ */

//double mrk_zyam = d2N_MRK(pT, qT, 0., yp, yq, rts);
//double qmrk_zyam = d2N_QMRK(pT, qT, 0., yp, yq, rts);
//double bfkl_zyam = d2N_BFKL(pT, qT, 0., yp, yq, rts);

//fprintf(out, "# p_T=%10.2e\tq_T=%10.2e\ty_p=%10.2e\ty_q=%10.2e\n",pT,qT,yp,yq);  		
//fprintf(out, "# x_P: %10.5e\n", (pT*exp(-yp)+qT*exp(-yq))/rts );
//fprintf(out, "# x_T: %10.5e\n", (pT*exp(+yp)+qT*exp(+yq))/rts );
//fprintf(out, "# MRK zyam=%10.5e\tQMRK zyam=%10.5e\tBFKL zyam=%10.5e\n",\
//              mrk_zyam, qmrk_zyam, bfkl_zyam);

//for (int n = 0; n < nphi;  n++) {
//     pq = (n != nphi - 1 ? n*(pi)/(double)(nphi-1) : 3.14);
//  		fprintf(out,"%10.2e\t%10.5e\t%10.5e\t%10.5e\n", pq,\
//      d2N_MRK(pT, qT, pq, yp, yq, rts),\
//      d2N_QMRK(pT, qT, pq, yp, yq, rts),\
//      d2N_BFKL(pT, qT, pq, yp, yq, rts) );
//      fflush(out);
//}  


return 0;
}
