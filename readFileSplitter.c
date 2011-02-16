#include <stdio.h>
#include <stdlib.h>

int main(int argc, char **argv){
    if (argc < 2) {
        printf("Usage: %s [options] reads\nOutputs two files part1_reads and part2_reads for bowtie consumption\nOptions: -ro for reads only (no quality info)\n", argv[0]);
        return 0;
    }
    printf("Input file to split: %s\n", argv[argc - 1]);
    
    int hasQualityInfo = 1;
    
    if (argc == 3) {
        if(strcmp(argv[1], "-ro")){
            hasQualityInfo = 0;
        }else{
            printf("Could not find option %s\n", argv[1]);
            return 0;
        }
    }
    
    // attempt to open the input file
    FILE *ins = fopen(argv[argc - 1], "r");
    if(ins == NULL){
        printf("Error! Could not open input file: %s\n", argv[argc - 1]);
    }
    
    
    
    // open up output file
    FILE *fo1, *fo2;
    char fileName1[100] = "part1_", fileName2[100] = "part2_";
    strcat(fileName1, argv[argc - 1]);
    strcat(fileName2, argv[argc - 1]);
    fo1 = fopen(fileName1, "w");
    fo2 = fopen(fileName2, "w");
    if(fo1 == NULL || fo2 == NULL){
        printf("Error! Could not open output files: %s or %s\n", fileName1, fileName2);
    }
    
    // read in and output the files
    int keepGoing = 1;
    char seqName[256];
    char seq[256];
    char qual[256];
    char temp[5];
    if(hasQualityInfo == 1){
        while(keepGoing > 0){
            // first read
            keepGoing = fscanf( ins, "%255s%255s%5s%255s", seqName, seq, temp, qual);
            if(keepGoing > 0){
                fprintf(fo1, "%s\n%s\n+\n%s\n", seqName, seq, qual);
                keepGoing = fscanf( ins, "%255s%255s%5s%255s", seqName, seq, temp, qual);
                fprintf(fo2, "%s\n%s\n+\n%s\n", seqName, seq, qual);
            }
        }
        close(fo1);
        close(fo2);
        close(ins);
    }else{
        while(keepGoing > 0){
            // first read
            keepGoing = fscanf( ins, "%255s%255s", seqName, seq);
            if(keepGoing > 0){
                fprintf(fo1, "%s\n%s\n", seqName, seq);
                keepGoing = fscanf( ins, "%255s%255s", seqName, seq);
                fprintf(fo2, "%s\n%s\n", seqName, seq);
            }
        }
        close(fo1);
        close(fo2);
        close(ins);
    }
    return 1;
}