// ALE.c

// (C) 2010 Scott Clark, JGI, LBNL

// Compiling instructions
// $ cc -g -O2 ALEv2.c -o ALE -lz -lm
// $ ./ALE <parameter file>
// requires: ALE.h, kseq.h, zlib.h, geneTree.h

// bowtie -t -I 0 -X 400 --fr -a -l 10 -v 3 -e 300 -S --threads 6 --sam-nohead Ecoli_complete -1 part1_EcoliReads2800k.fna -2 part2_EcoliReads2800k.fna Ecoli_complete.map.sam

#include "ALE.h"

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "ALEhelpers.h"
#include "ALElike.h"
#include "geneTree.h"

int main(int argc, char **argv){
  
    // Handles input variations
    if (argc < 2) {
        printf(USAGE, argv[0]);
        printf("%s", SHORT_OPTIONS);
        return 0;
    }
    if (argc < 4) {
        if(argc <= 1 || (argv[1][0] == '-' && argv[1][1] == 'h')){
            printf("%s", WELCOME_MSG);
            printf(USAGE, argv[0]);
            printf("\n%s", LONG_OPTIONS);
            return 0;
        }else{
            printf(USAGE, argv[0]);
            printf("%s", SHORT_OPTIONS);
            return 0;
        }
    }
    
    int options;
    int kmerLen = 4;
    char placementOut[256];
    samfile_t *placementBam = NULL;
    *placementOut = '\0';
    libraryParametersT *libParams = NULL; // TODO load/save this data structure
    double outlierFraction = 0.02;
    int qOff = -1;
    
    if(argc > 5) { // look for command line options
        for(options = 1; options < argc - 4; options++){ // search over all options
           if(strcmp(argv[options], "-kmer") == 0){
                kmerLen = atoi(argv[options+1]);
                if(kmerLen > 20){
                    printf("-kmer option of %i not in range [2,20], set to default [4].\n", kmerLen);
                    kmerLen = 4;
                }
                options++;
            }else if(strcmp(argv[options], "-qOff") == 0){
                qOff = atoi(argv[options+1]);
                if(qOff != 64 && qOff != 33 && qOff != 0){
                    printf("-qOff option of %i not in set [33,64], will be set to 33.\n", qOff);
                    qOff = 33;  // SAM/BAM specification is for ascii - 33.
                }
                options++;
            }else if(strcmp(argv[options], "-pl") == 0){
            	strcpy(placementOut, argv[options+1]);
                options++;
            } else{
                printf("Could not find option %s\n", argv[options]);
            }
        }
    }
    
    // input and output files
    printf("BAM file: %s\n", argv[argc - 3]);
    printf("Assembly fasta file: %s\n", argv[argc - 2]);
    printf("ALE Output file: %s\n", argv[argc - 1]);

    // attempt to open the bam input file
    samfile_t *ins = openSamOrBam(argv[argc - 3]);
    
    printf("Reading in assembly...\n");
    assemblyT *theAssembly = loadAssembly(argv[argc - 2]);
    if (!validateAssemblyIsSameAsAlignment(ins->header, theAssembly)) {
        printf("Error! Assembly fasta %s does not match alignment file %s!\n", argv[argc-2], argv[argc-3]);
        exit(1);
    }
    
    if (*placementOut != '\0') {
        printf("Placement file: %s\n", placementOut);
        char *mode = "wbu";
        if (placementOut[0] != '-' )
            mode[2] = '\0'; // compress bam if not part of a pipe
        placementBam = samopen(placementOut, mode, ins->header);
        if(placementBam == 0){
            printf("Error! Could not write to the placement BAM file: %s\n", placementOut);
            exit(1);
        }
    }

    printf("Reading in the map and computing statistics...\n");

    // calculate the insert mean/std if not given
    if(libParams == NULL){
         printf("Insert length and std not given, will be calculated from input map.\n");

        libParams = computeLibraryParameters(ins, outlierFraction, qOff);

        // close and re-open bam file
        samclose(ins);
        ins = openSamOrBam(argv[argc - 3]);
    }

    printf("Calculating GC content of reference over average read size\n");
    calculateGCcont(theAssembly, libParams->avgReadSize);

    // place reads and compute statistics on the assembly
    printf("Computing read placements and depths\n");
    computeReadPlacements(ins, theAssembly, libParams, placementBam);
    
    // compute statistics on assembly
    printf("Computing k-mer statistics...\n");
    computeKmerStats(theAssembly, kmerLen);
    printf("Done computing k-mer statistics.\n");
    
    printf("Computing depth statistics...\n");
    computeDepthStats(theAssembly);
    printf("Done computing depth statistics.\n");
    
    FILE *out = fopen(argv[argc - 1], "w");
    if(out == NULL){
        printf("Error! Could not open output file: %s\n", argv[argc - 1]);
    }
    writeToOutput(theAssembly, out);
    fclose(out);
    
    printf("Output is in file: %s\n", argv[argc - 1]);

    if (placementBam != NULL) {
        printf("Closing placement file\n");
        samclose(placementBam);
    }

    printf("Closing input file\n");
    samclose(ins);
    
    if (libParams != NULL)
        free(libParams);

    freeAssembly(theAssembly);
    return 0;
}
