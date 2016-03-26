#define buffSize_RP 80

#include <stdio.h>
#include <string.h>

//using namespace std;

//these are used to see if we missed initializing anything
int global_variable_counter = 0;
int * global_variable_init ;

#define mydef_int(name) int name
#define mydef_double(name) double name
#define mydef_char(name) char * name
#define X(type, name, format) mydef_##type(name) ; 
   X_FIELDS
#undef X

#define setvar_int(name, tokValu) name = atol (tokValu)
#define setvar_double(name, tokValu) name = atof (tokValu)
#define setvar_char(name, tokValu)  name = (char *)malloc(buffSize_RP*sizeof(char)); strcpy(name,tokValu)

void fill_X_LIST(char *tokName, char *tokValu)
{
   int match = 0;
   int i = 0;
   #define X(type, name, format) \
      if (strcmp (tokName, #name) == 0) { \
         setvar_##type(name, tokValu) ;\
      global_variable_init[i] = 1; \
      match = 1; \
      } \
      i ++ ; 
   X_FIELDS
   #undef X
        
   if (match == 0) 
      printf("Warning: variable \"%s\" is present in input file but not used.\n",tokName); \
}

void print_X_LIST()
{
   int i = 0;
   #define X(type, name, format) \
      if ( global_variable_init[i] == 0 ) \
      printf("Warning: global variable \"%s\" has not been initialized\n",#name); \
      i ++ ; 
   X_FIELDS
   #undef X
   
   i = 0;
   printf("---Begin Global Parameter List---\n");
   #define X(type, name, format) \
      if ( global_variable_init[i] == 1 ) { \
            printf("%s: "format"\n",#name,name); \
         } \
   i ++ ; 
   X_FIELDS
   #undef X
   printf("---End Global Parameter List---\n");
}

int ReadParameters( FILE *fp )
{

   #define X(type, name, format) global_variable_counter ++ ; 
      X_FIELDS
   #undef X

   global_variable_init = (int *) malloc(global_variable_counter * sizeof (int));
   int i;
   for (i = 0; i < global_variable_counter; i++)
      global_variable_init[i] = 0;

   char * tokName, *tokValu;
   char buff[buffSize_RP];
   while (1) {
      if (fgets(buff, buffSize_RP, fp) == NULL)
      {
         return 0;
         exit(2);
      }

      if (feof (fp) ) break;

      tokName = strtok(buff, " \t\n");
      tokValu = strtok(NULL, " \t\n");
      if (! tokValu) break;

      fill_X_LIST(tokName, tokValu);

   }
  
   print_X_LIST();
   return 1;
}

char fn_filename[256];
#define OpenFile(fn, tag, name, mode)   strcpy(fn_filename, tag); strcat       (fn_filename, name); FILE * fn = fopen( fn_filename, mode); if (fn == NULL)     {printf("Error: file \"%s\" could not be opened.\n\n",fn_filename); exit(1); }

