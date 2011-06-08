#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int main(int argc, char **argv){
  
  if (argc < 4) {
    printf("Usage: %s assemblyFile.fasta runningWindow output\n(C) 2011 Scott Clark\nFinds the running window GC content in a fasta file\nJust sequence! No gene name/info!\n", argv[0]);
    return 0;
  }
  
  // attempt to open the first input file
  FILE *ins = fopen(argv[argc - 3], "r");
  if(ins == NULL){
      printf("Error! Could not open the first input file: %s\n", argv[argc - 3]);
  }
  
  // attempt to open the output file
  FILE *fout = fopen(argv[argc - 1], "w");
  if(fout == NULL){
      printf("Error! Could not open the output file: %s\n", argv[argc - 1]);
  }
  
  int runningWindow = atoi(argv[argc - 2]);
  int i;
  char base;
  char *backTrace = malloc(sizeof(char)*runningWindow);
  int gcTot = 0;
  int keepGoing;
  
  for(i = 0; i < runningWindow; i++){
    keepGoing = fscanf("%c", base);
//     if(base == '\0' || base == '\n'){
//       i--;
//     }
    backTrace[i] = base;
    if(base == 'G' || base == 'g' || base == 'C' || base == 'c'){
      gcTot++;
    }
  }
  
  printf("%d\n", (float)gcTot/((float)runningWindow));
  
  while(keepGoing > 0){
    keepGoing = fscanf("%c", base);
    if(base == '\0' || base == '\n'){
      i--;
    }
    if(base == 'G' || base == 'g' || base == 'C' || base == 'c'){
      gcTot++;
    }
  }
}
