#include "main.h"

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

gsl_set_error_handler_off ();

TabulateAlpha();

ReadInWF(A1, A2, wfTAG);

strcat(tag,"output/");
strcat(tag,inputFileName);
OpenFile(outsingle, tag, "_single.dat", "w");
//TabulateSingle(outsingle,rts);

//double rT;
//for (rT=.28125; rT<18. ; rT+=0.0625)
//{
 //   printf("%10.5f\t%10.5e\n",rT,dNd2pTdy(rT,0.,rts));
//}
 
OpenFile(outglasma, tag, "_double.dat", "w");
double pT, qT, phi;
for (pT=.28125; pT<18. ; pT+=0.0625)
for (qT=.28125; qT<18. ; qT+=0.0625)
for (phi=0; phi<3.1416 ; phi+=.1963)
{
fprintf(outglasma,"%10.5f\t%10.5f\t%10.3e\t%10.5e\t%10.5e\t%10.5e\t%10.5e\t%10.5e\t%10.5e\n",pT,qT,phi,\
        d2Nglasma0(pT, qT, phi, 0., 0.,rts),\
        d2Nglasma1(pT,0.,0.,rts),\
        d2N_MRK(pT, qT, phi, 0., 0., rts),\
        d2N_simple_jet(pT, qT, phi, 0., 0., rts),\
        dNd2pTdy(pT,0.,rts),\
        dNd2pTdy(qT,0.,rts) );
}

return 0;

    clock_t begin, end;
    double time_spent;
    begin = clock();
    TabulateGlasma(outglasma,rts);
    end = clock();
    time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
    printf("time elapsed (seconds): %10.2e\n", time_spent);

return 0;
}
