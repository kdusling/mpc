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
TabulateSingle(outsingle,rts);

OpenFile(outglasma, tag, "_glasma0.dat", "w");

    clock_t begin, end;
    double time_spent;
    begin = clock();
    TabulateGlasma(outglasma,rts);
    end = clock();
    time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
    printf("time elapsed (seconds): %10.2e\n", time_spent);

return 0;
}
