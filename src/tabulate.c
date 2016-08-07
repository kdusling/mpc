#include "tabulate.h"

int main(int argc, char **argv)
{

   // Initialize the MPI environment
   MPI_Init(NULL, NULL);
   // Find out rank, size
   MPI_Comm_rank(MPI_COMM_WORLD, &myid);
   MPI_Comm_size(MPI_COMM_WORLD, &np);

   if (np < 3)
   {
    printf("Error: np must be at least 3.\n");
    return 0;
   }

   //if we are the master
   if (myid == 0) { master_main(argc, argv); }
  
   //if we are the slave 
   if (myid != 0) { slave_main(argc, argv); }

   MPI_Finalize();
   return 0;
}

int master_main(int argc, char **argv)
{
//holds name of inputfile
char inputFileName[80];

//get input file name
   if (argc != 2)
   {
      printf("Error: wrong input arguments passed to command line.\n");
      printf("Arguments should be: \n");
      printf("tag\n");
      printf("Example:\n");
      printf("./tabulate.exe 01_01\n");
      exit(2);
   } else {
      strcpy(inputFileName,argv[1]);
   }

   //read in parameters
   OpenFile(input, inputFileName, ".input", "r");
   int ierr = ReadParameters (input);
   if (ierr == 1)
      printf("Error in ReadParameters\n");

   strcat(tag,"output/");
   strcat(tag,inputFileName); 

   char filename[256];
   strcpy(filename,"_bfkl.dat");
   OpenFile(output, tag, filename, "w");

   int i;
   double *phi;
   phi = (double *)malloc( (np-1)*sizeof(double));
  
   for (i = 0; i < np-1; i++){
      phi[i] = (double)i/(double)(np-2)*3.14159;
   }

   double slres;
   double yp, yq, pT, qT;
   double vals[5]; 
   yp = -1.5; yq = +1.5;
   
   for (pT=0.96875; pT<18. ; pT+=0.25)
   for (qT=0.96875; qT<18. ; qT+=0.25)
   {   
      //this loop sends out phi values to all the processors
      for(i=1; i <= np-1; i++) {
         vals[0] = pT; vals[1] = qT;
         vals[2] = phi[i-1];
         vals[3] = yp; vals[4] = yq;

         MPI_Send(vals, 5, MPI_DOUBLE, i, PIPE_MSG, MPI_COMM_WORLD);
      }
      
      for(i=1; i <= np-1; i++) {
    	  MPI_Recv(&slres, 1, MPI_DOUBLE, i, RET_MSG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
         fprintf(output,"%10.2f\t%10.2f\t%10.2f\t%10.2f\t%10.4f\t%10.5e\n",\
         yp, yq, pT, qT, phi[i-1], slres);
      }
      //finished computation of all phi values
   fflush(output);
   }

   //Tell processors we are done
   for(i=1; i <= np-1; i++) {
      MPI_Send(vals, 5, MPI_DOUBLE, i, END_MSG, MPI_COMM_WORLD);
   }

return 0;
}

int slave_main(int argc, char **argv)
{

//holds name of inputfile
char inputFileName[80];
strcpy(inputFileName,argv[1]);

//read in parameters
OpenFile(input, inputFileName, ".input", "r");
ReadParameters (input);

gsl_set_error_handler_off ();

TabulateAlpha();

ReadInWF(A1, A2, wfTAG);
ReadInBFKL();

double myres;
double vals[5]; 

   for(;;){

      MPI_Recv(vals, 5, MPI_DOUBLE, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &Status);
      if (Status.MPI_TAG == END_MSG)  break;
      
      myres = d2N_BFKL( vals[0], vals[1], vals[2], vals[3], vals[4], rts );

      MPI_Send(&myres, 1, MPI_DOUBLE, 0, RET_MSG, MPI_COMM_WORLD);
   }

return 0;
}
